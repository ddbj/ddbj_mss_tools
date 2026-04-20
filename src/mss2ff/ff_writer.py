"""Generate DDBJ flat file entries from parsed MSS annotation data."""

from __future__ import annotations

import re
import sys
import textwrap
from datetime import date
from io import StringIO
from typing import Optional

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature

from .ann_parser import CommonBlock, EntryBlock, Feature
from .location import expand_location, parse_mss_location
from .taxonomy import get_lineage

# Qualifiers that use integer format (no quotes)
_INT_QUALS = frozenset({"codon_start", "transl_table", "estimated_length"})

# Width constants
_LINE_WIDTH = 80
_FIELD_INDENT = "            "          # 12 spaces (LOCUS-level field continuation)
_QUAL_INDENT = "                     "  # 21 spaces (feature qualifier continuation)
_FEAT_PREFIX = 5                        # spaces before feature name
_FEAT_NAME_WIDTH = 16                   # feature name field width


# ── Text helpers ──────────────────────────────────────────────────────────────

def _wrap(text: str, indent: str, first_prefix: str = "") -> list[str]:
    """Wrap text to _LINE_WIDTH, using first_prefix for the first line."""
    first_avail = _LINE_WIDTH - len(first_prefix)
    cont_avail = _LINE_WIDTH - len(indent)

    if not text:
        return [first_prefix]

    words = text.split()
    lines: list[str] = []
    current: list[str] = []
    avail = first_avail

    for word in words:
        needed = len(word) if not current else len(word) + 1
        if current and len(" ".join(current)) + needed > avail:
            lines.append((first_prefix if not lines else indent) + " ".join(current))
            current = [word]
            avail = cont_avail
        else:
            current.append(word)

    if current:
        lines.append((first_prefix if not lines else indent) + " ".join(current))

    return lines


def _format_date(d: date) -> str:
    return d.strftime("%d-%b-%Y").upper()


def _authors_str(ab_names: list[str]) -> str:
    """Format author list: 'A, B and C' or 'A and B' or 'A'."""
    names = ab_names
    if len(names) == 0:
        return ""
    if len(names) == 1:
        return names[0]
    return ", ".join(names[:-1]) + " and " + names[-1]


# ── Qualifier formatting ──────────────────────────────────────────────────────

def _format_qualifier(key: str, value: str) -> list[str]:
    """Return formatted qualifier line(s) with 21-space indent."""
    indent = _QUAL_INDENT
    if not value:
        return [f"{indent}/{key}"]
    if key in _INT_QUALS or (value.lstrip("-").isdigit()):
        return [f"{indent}/{key}={value}"]

    prefix = f"{indent}/{key}=\""
    suffix = '"'
    # Characters available on the first line for value text
    first_avail = _LINE_WIDTH - len(prefix) - len(suffix)
    cont_avail = _LINE_WIDTH - len(indent) - len(suffix)

    if len(value) <= first_avail:
        return [f"{prefix}{value}{suffix}"]

    # For translation: wrap at character boundaries (no word split needed)
    if key == "translation":
        lines = []
        chunk = value[:first_avail]
        lines.append(f"{prefix}{chunk}")
        remaining = value[first_avail:]
        while remaining:
            chunk = remaining[:cont_avail]
            remaining = remaining[cont_avail:]
            if not remaining:
                lines.append(f"{indent}{chunk}{suffix}")
            else:
                lines.append(f"{indent}{chunk}")
        return lines

    # For other qualifiers: wrap at word boundaries
    lines = []
    remaining = value
    first = True
    while remaining:
        avail = first_avail if first else cont_avail
        if len(remaining) <= avail:
            tail = f"{prefix if first else indent}{remaining}{suffix}"
            lines.append(tail)
            break
        # Try to break at last space within avail chars
        chunk = remaining[:avail]
        sp = chunk.rfind(" ")
        if sp > 0:
            chunk = chunk[:sp]
        if first:
            lines.append(f"{prefix}{chunk}")
        else:
            lines.append(f"{indent}{chunk}")
        remaining = remaining[len(chunk):].lstrip(" ")
        first = False

    return lines


# ── Feature location & meta expansion ────────────────────────────────────────

def _resolve_meta(value: str, source_quals: dict[str, str], entry_id: str, seq_len: int) -> str:
    """Replace @@[key]@@ placeholders with actual values from source qualifiers."""
    def replacer(m: re.Match) -> str:
        k = m.group(1)
        if k == "entry":
            return entry_id
        return source_quals.get(k, m.group(0))

    result = re.sub(r"@@\[([^\]]+)\]@@", replacer, value)
    # Replace bare E used as sequence-length marker (word boundary)
    result = re.sub(r"\bE\b", str(seq_len), result)
    return result


# ── Feature section ───────────────────────────────────────────────────────────

def _format_feature_block(
    feat: Feature,
    source_quals: dict[str, str],
    entry_id: str,
    seq_len: int,
    parent_seq: Optional[Seq],
) -> list[str]:
    """Format a single feature block (feature line + qualifier lines)."""
    location = expand_location(feat.location, seq_len, entry_id)
    name_field = f"{'':>{_FEAT_PREFIX}}{feat.name:<{_FEAT_NAME_WIDTH}}"
    feat_line = f"{name_field}{location}"
    lines = [feat_line]

    # Pre-compute source qualifier map for meta expansion
    translation_needed = (
        feat.name == "CDS"
        and not feat.has_qualifier("pseudo")
        and not feat.has_qualifier("pseudogene")
        and parent_seq is not None
    )

    for key, val in feat.qualifiers:
        if key == "ff_definition":
            continue  # used for DEFINITION, not written as qualifier
        resolved = _resolve_meta(val.strip(), source_quals, entry_id, seq_len)

        # assembly_gap: expand estimated_length=known to actual gap length
        if (
            feat.name == "assembly_gap"
            and key == "estimated_length"
            and resolved.lower() == "known"
        ):
            resolved = str(_calc_location_length(location))

        lines.extend(_format_qualifier(key, resolved))

    if translation_needed:
        translation = _translate_cds(feat, location, parent_seq, seq_len)
        if translation:
            lines.extend(_format_qualifier("translation", translation))

    return lines


def _calc_location_length(loc_str: str) -> int:
    """Return total base count covered by a location string."""
    try:
        return len(parse_mss_location(loc_str))
    except Exception:
        return 0


def _translate_cds(feat: Feature, location: str, parent_seq: Seq, seq_len: int) -> str:
    """Translate CDS feature, returning protein string (no stop codon)."""
    try:
        from .translate_with_transl_except import translate_cds_with_transl_except

        bio_loc = parse_mss_location(location)
        # Build minimal SeqFeature for translation
        qualifiers: dict[str, list[str]] = {}
        for k, v in feat.qualifiers:
            qualifiers.setdefault(k, []).append(v.strip())

        sf = SeqFeature(location=bio_loc, type="CDS", qualifiers=qualifiers)
        protein = translate_cds_with_transl_except(sf, parent_seq)
        return str(protein).rstrip("*")
    except Exception as exc:  # noqa: BLE001
        print(f"[ff_writer] Warning: translation failed for {feat.location}: {exc}", file=sys.stderr)
        return ""


# ── LOCUS line ────────────────────────────────────────────────────────────────

def _locus_line(
    name: str,
    seq_len: int,
    topology: str,
    division: str,
    file_date: date,
) -> str:
    topo = "circular" if topology == "circular" else "linear  "
    d = _format_date(file_date)
    return f"LOCUS       {name:<21}{seq_len:>7} bp    DNA     {topo} {division:<3} {d}"


# ── Header fields ─────────────────────────────────────────────────────────────

def _definition_lines(definition: str) -> list[str]:
    prefix = "DEFINITION  "
    cont = _FIELD_INDENT
    return _wrap(definition, cont, prefix)


def _dblink_lines(common: CommonBlock) -> list[str]:
    lines = []
    prefix = "DBLINK      "
    cont = _FIELD_INDENT
    if common.dblink.project:
        lines.append(f"{prefix}BioProject:{common.dblink.project}")
        prefix = cont
    if common.dblink.biosample:
        lines.append(f"{prefix}BioSample:{common.dblink.biosample}")
        prefix = cont
    if common.dblink.sra:
        sra_text = "Sequence Read Archive: " + ", ".join(common.dblink.sra)
        lines.extend(_wrap(sra_text, cont, prefix))
        prefix = cont
    for label, val in common.dblink.extra:
        lines.extend(_wrap(f"{label}:{val}", cont, prefix))
        prefix = cont
    return lines


def _keywords_line(keywords: list[str]) -> str:
    if keywords:
        return "KEYWORDS    " + "; ".join(keywords) + "."
    return "KEYWORDS    ."


def _source_organism_lines(
    organism: str,
    lineage: str,
) -> list[str]:
    lines = []
    lines.extend(_wrap(organism, _FIELD_INDENT, "SOURCE      "))
    lines.extend(_wrap(organism, " " * 12, "  ORGANISM  "))
    if lineage:
        lines.extend(_wrap(lineage, " " * 12, " " * 12))
    return lines


# ── REFERENCE blocks ──────────────────────────────────────────────────────────

def _reference_1_lines(
    common: CommonBlock,
    seq_len: int,
    submission_date: date,
) -> list[str]:
    """Reference 1: Direct Submission (SUBMITTER information)."""
    sub = common.submitter
    cont = " " * 12
    lines = []
    lines.append(f"REFERENCE   1  (bases 1 to {seq_len})")

    # AUTHORS line (omit if no names but CONSRTM is present)
    if sub.ab_names:
        lines.extend(_wrap(_authors_str(sub.ab_names), cont, "  AUTHORS   "))
    if sub.consrtm:
        lines.extend(_wrap(sub.consrtm, cont, "  CONSRTM   "))

    lines.append("  TITLE     Direct Submission")

    # JOURNAL
    date_str = _format_date(submission_date)
    lines.append(f"  JOURNAL   Submitted ({date_str})")

    if sub.contact:
        lines.extend(_wrap(f"Contact:{sub.contact}", cont, cont))

    # Affiliation: [department, ]institute; (department first per DDBJ format)
    affil_parts = []
    if sub.department:
        affil_parts.append(sub.department)
    if sub.institute:
        affil_parts.append(sub.institute)
    if affil_parts:
        lines.extend(_wrap(", ".join(affil_parts) + ";", cont, cont))

    # Location: street, city, state zip, country
    loc_parts = []
    if sub.street:
        loc_parts.append(sub.street)
    if sub.city:
        loc_parts.append(sub.city)
    state_zip = ""
    if sub.state:
        state_zip = sub.state
    if sub.zip:
        state_zip = (state_zip + " " if state_zip else "") + sub.zip
    if state_zip:
        loc_parts.append(state_zip)
    if sub.country:
        loc_parts.append(sub.country)
    if loc_parts:
        lines.extend(_wrap(", ".join(loc_parts), cont, cont))

    # URL    : aligned with Contact: (both 7 chars before colon)
    if sub.url:
        lines.append(f"{cont}URL    :{sub.url}")

    return lines


def _reference_n_lines(ref_num: int, ref) -> list[str]:
    """Format a non-submission reference.

    JOURNAL format by status:
      Unpublished (no year)  : Unpublished.
      Unpublished (with year): Unpublished. (year)
      In press               : journal (year) In press
      Published              : journal volume, start_page[-end_page] (year)
    """
    cont = " " * 12
    lines = []
    lines.append(f"REFERENCE   {ref_num}  ")
    # AUTHORS line (omit if no names but CONSRTM is present)
    if ref.ab_names:
        lines.extend(_wrap(_authors_str(ref.ab_names), cont, "  AUTHORS   "))
    if ref.consrtm:
        lines.extend(_wrap(ref.consrtm, cont, "  CONSRTM   "))
    lines.extend(_wrap(ref.title, cont, "  TITLE     "))

    status = (ref.status or "").strip().lower()

    if status == "in press":
        journal_str = f"{ref.journal} ({ref.year}) In press"
    elif status == "published":
        pages = ref.from_page
        if ref.to_page:
            pages = f"{pages}-{ref.to_page}"
        journal_str = f"{ref.journal} {ref.volume}, {pages} ({ref.year})"
    else:
        # Unpublished (default)
        journal_str = f"Unpublished. ({ref.year})" if ref.year else "Unpublished."

    lines.extend(_wrap(journal_str, cont, "  JOURNAL   "))

    if ref.doi:
        lines.extend(_wrap(f"DOI:{ref.doi}", cont, "  REMARK    "))
    elif ref.pubmed:
        lines.extend(_wrap(f"PUBMED:{ref.pubmed}", cont, "  REMARK    "))

    return lines


# ── COMMENT section ───────────────────────────────────────────────────────────

def _comment_lines(common: CommonBlock) -> list[str]:
    """Build COMMENT block with regular comments then ST_COMMENTs."""
    all_regular: list[str] = []
    for block in common.comment_blocks:
        all_regular.extend(block)

    st_blocks = common.st_comments
    if not all_regular and not st_blocks:
        return []

    prefix = "COMMENT     "
    cont = " " * 12
    lines: list[str] = []

    for line_text in all_regular:
        if not lines:
            lines.extend(_wrap(line_text, cont, prefix))
        else:
            lines.extend(_wrap(line_text, cont, cont))

    for st in st_blocks:
        header = f"##{st.tagset_id}-START##"
        footer = f"##{st.tagset_id}-END##"

        if lines:
            lines.append(cont)  # blank separator line
        else:
            prefix = cont

        if not lines:
            lines.extend(_wrap(header, cont, "COMMENT     "))
        else:
            lines.append(f"{cont}{header}")

        for key, val in st.fields:
            field_prefix = f"{cont}{key:<22}:: "
            field_cont = " " * len(field_prefix)
            lines.extend(_wrap(val, field_cont, field_prefix))

        lines.append(f"{cont}{footer}")

    return lines


# ── Accession number utilities ────────────────────────────────────────────────

# Accession patterns — longest prefix match first to avoid ambiguity.
# Serial number is ≥6 digits (variable length).
_ACC_PATTERNS = [
    re.compile(r"^([A-Z]{6}\d{2})(\d{6,})$"),   # AAXJEM010000001 …
    re.compile(r"^([A-Z]{4}\d{2})(\d{6,})$"),   # AAXJ010000001 …
    re.compile(r"^([A-Z]{2})(\d{6,})$"),         # AA000001 …
]


def _parse_accession(acc: str) -> tuple[str, str]:
    """Split accession into (prefix, serial_str).

    AA000001        → ('AA',       '000001')
    AAXJ010000001   → ('AAXJ01',   '0000001')
    AAXJEM010000001 → ('AAXJEM01', '0000001')

    The serial portion can be any length ≥ 6 digits; the padding width is
    preserved when incrementing.
    """
    acc = acc.strip().upper()
    for pat in _ACC_PATTERNS:
        m = pat.match(acc)
        if m:
            return m.group(1), m.group(2)
    raise ValueError(
        f"Cannot parse accession {acc!r}. "
        "Expected format: AA000001, AAXJ010000001, or AAXJEM010000001 "
        "(2/4/6 letters [+ 2-digit version] + ≥6 digit serial)."
    )


def _accession_at(prefix: str, serial_str: str, offset: int) -> str:
    """Return the accession number at *offset* steps from (prefix, serial_str)."""
    new_serial = int(serial_str) + offset
    return f"{prefix}{new_serial:0{len(serial_str)}d}"


# ── BASE COUNT ────────────────────────────────────────────────────────────────

def _base_count_line(seq: str) -> str:
    s = seq.lower()
    a = s.count("a")
    c = s.count("c")
    g = s.count("g")
    t = s.count("t")
    # Format matches DDBJ reference: 7 spaces before each count
    return f"BASE COUNT       {a} a       {c} c       {g} g       {t} t"


# ── ORIGIN section ────────────────────────────────────────────────────────────

def _origin_lines(seq: str) -> list[str]:
    lines = ["ORIGIN      "]
    seq = seq.lower()
    for i in range(0, len(seq), 60):
        pos = i + 1
        chunk = seq[i:i + 60]
        groups = " ".join(chunk[j:j + 10] for j in range(0, len(chunk), 10))
        lines.append(f"{pos:>9} {groups}")
    return lines


# ── Main entry builder ────────────────────────────────────────────────────────

def build_entry(
    entry: EntryBlock,
    common: CommonBlock,
    seq: str,
    division: str = "UNK",
    submission_date: Optional[date] = None,
    file_date: Optional[date] = None,
    email: str = "mss2ff@ddbj.nig.ac.jp",
    no_taxonomy: bool = False,
    accession: Optional[str] = None,
) -> str:
    """Build one DDBJ flat file record and return it as a string.

    *accession* overrides the entry ID in LOCUS/ACCESSION/VERSION lines.
    DEFINITION always comes from ff_definition (no accession substitution).
    """
    today = date.today()
    if submission_date is None:
        submission_date = today
    if file_date is None:
        file_date = today

    seq_len = len(seq)
    entry_id = entry.entry_id
    locus_acc = accession or entry_id   # used for LOCUS name, ACCESSION, VERSION
    parent_seq = Seq(seq) if seq else None

    # Collect source qualifiers for meta expansion
    src = entry.source_feature
    source_quals: dict[str, str] = {}
    if src:
        for k, v in src.qualifiers:
            source_quals[k] = v.strip()

    organism = source_quals.get("organism", "").strip()

    # Fetch taxonomy
    taxon_id = ""
    lineage = ""
    if organism and not no_taxonomy:
        taxon_id, lineage = get_lineage(organism, email=email)

    # Build DEFINITION from ff_definition in source feature
    definition = ""
    if src:
        ff_def = src.get_qualifier("ff_definition") or ""
        definition = _resolve_meta(ff_def, source_quals, entry_id, seq_len).strip()
    if not definition:
        definition = f"{organism} DNA."

    # Keywords (entry-level overrides common-level)
    keywords = entry.keywords or common.keywords

    out = StringIO()
    w = out.write

    def wl(line: str = ""):
        out.write(line + "\n")

    # LOCUS  (uses accession if assigned, otherwise entry_id)
    wl(_locus_line(locus_acc, seq_len, entry.topology, division, file_date))

    # DEFINITION (period at end) — never contains the accession number
    def_str = definition if definition.endswith(".") else definition + "."
    for line in _definition_lines(def_str):
        wl(line)

    # ACCESSION / VERSION  (uses accession if assigned)
    wl(f"ACCESSION   {locus_acc}")
    wl(f"VERSION     {locus_acc}.1")

    # DBLINK
    for line in _dblink_lines(common):
        wl(line)

    # KEYWORDS
    wl(_keywords_line(keywords))

    # SOURCE / ORGANISM
    for line in _source_organism_lines(organism, lineage):
        wl(line)

    # REFERENCE 1 (Direct Submission)
    for line in _reference_1_lines(common, seq_len, submission_date):
        wl(line)

    # REFERENCE 2+ (literature references)
    for idx, ref in enumerate(common.references, start=2):
        for line in _reference_n_lines(idx, ref):
            wl(line)

    # COMMENT
    comment_block = _comment_lines(common)
    for line in comment_block:
        wl(line)

    # FEATURES
    wl("FEATURES             Location/Qualifiers")

    # Add taxon db_xref to source feature output
    src_taxon_injected = False

    for feat in entry.features:
        expanded_loc = expand_location(feat.location, seq_len, entry_id)
        feat_for_write = Feature(
            name=feat.name,
            location=expanded_loc,
            qualifiers=feat.qualifiers,
        )

        if feat.name == "source" and taxon_id and not src_taxon_injected:
            # Inject db_xref="taxon:XXXXX" into source feature
            modified_quals = list(feat.qualifiers)
            if not any(k == "db_xref" for k, _ in modified_quals):
                # Insert after mol_type if present, else append
                insert_pos = next(
                    (i for i, (k, _) in enumerate(modified_quals) if k == "mol_type"),
                    len(modified_quals) - 1,
                )
                modified_quals.insert(insert_pos + 1, ("db_xref", f"taxon:{taxon_id}"))
            feat_for_write = Feature(
                name=feat.name,
                location=expanded_loc,
                qualifiers=modified_quals,
            )
            src_taxon_injected = True

        for line in _format_feature_block(feat_for_write, source_quals, entry_id, seq_len, parent_seq):
            wl(line)

    # BASE COUNT
    wl(_base_count_line(seq))

    # ORIGIN
    for line in _origin_lines(seq):
        wl(line)

    wl("//")

    return out.getvalue()


# ── File writer ───────────────────────────────────────────────────────────────

def write_ff(
    common: CommonBlock,
    entries: list[EntryBlock],
    sequences: dict[str, str],
    output,
    division: str = "UNK",
    submission_date: Optional[date] = None,
    file_date: Optional[date] = None,
    email: str = "mss2ff@ddbj.nig.ac.jp",
    no_taxonomy: bool = False,
    start_accession: Optional[str] = None,
) -> None:
    """Write DDBJ flat file to *output* (file object or path string).

    If *start_accession* is given (e.g. 'AAXJ010000001'), each entry is
    assigned a sequential accession number starting from that value.
    """
    close_after = False
    if isinstance(output, (str, __import__("pathlib").Path)):
        output = open(output, "w", encoding="utf-8")
        close_after = True

    # Pre-compute accession list if a starting accession is given
    acc_prefix: Optional[str] = None
    acc_serial: Optional[str] = None
    if start_accession:
        try:
            acc_prefix, acc_serial = _parse_accession(start_accession)
        except ValueError as exc:
            print(f"[ff_writer] Warning: {exc}", file=sys.stderr)

    try:
        for idx, entry in enumerate(entries):
            seq = sequences.get(entry.entry_id, "")
            if not seq:
                # Try partial-match fallback
                for sid, s in sequences.items():
                    if entry.entry_id in sid or sid in entry.entry_id:
                        seq = s
                        break
            if not seq:
                print(
                    f"[ff_writer] Warning: no sequence found for entry {entry.entry_id!r}",
                    file=sys.stderr,
                )

            accession = (
                _accession_at(acc_prefix, acc_serial, idx)
                if acc_prefix is not None
                else None
            )

            record_str = build_entry(
                entry=entry,
                common=common,
                seq=seq,
                division=division,
                submission_date=submission_date,
                file_date=file_date,
                email=email,
                no_taxonomy=no_taxonomy,
                accession=accession,
            )
            output.write(record_str)

        if acc_prefix is not None and entries:
            first_acc = _accession_at(acc_prefix, acc_serial, 0)
            last_acc = _accession_at(acc_prefix, acc_serial, len(entries) - 1)
            if len(entries) == 1:
                print(f"Assigned accession number: {first_acc}", file=sys.stderr)
            else:
                print(f"Assigned accession numbers: {first_acc}-{last_acc}", file=sys.stderr)
    finally:
        if close_after:
            output.close()
