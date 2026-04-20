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
from .taxonomy import get_lineage, map_tagset_id

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
        lines.extend(_format_qualifier(key, resolved))

    if translation_needed:
        translation = _translate_cds(feat, location, parent_seq, seq_len)
        if translation:
            lines.extend(_format_qualifier("translation", translation))

    return lines


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
        accessions = ", ".join(common.dblink.sra)
        lines.append(f"{prefix}Sequence Read Archive:{accessions}")
        prefix = cont
    for label, val in common.dblink.extra:
        lines.append(f"{prefix}{label}:{val}")
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
    lines = []
    lines.append(f"REFERENCE   1  (bases 1 to {seq_len})")
    authors = _authors_str(sub.ab_names)
    lines.extend(_wrap(authors, " " * 12, "  AUTHORS   "))
    lines.append("  TITLE     Direct Submission")

    # Build JOURNAL lines
    date_str = _format_date(submission_date)
    journal_prefix = "  JOURNAL   "
    cont = " " * 12

    lines.append(f"{journal_prefix}Submitted ({date_str})")
    if sub.contact:
        lines.extend(_wrap(f"Contact:{sub.contact}", cont, cont))

    # Address: institute[, department][; street, city[, state] zip][, country]
    address_parts = []
    if sub.institute:
        address_parts.append(sub.institute)
    if sub.department:
        address_parts.append(sub.department)

    addr = ", ".join(address_parts)
    location_parts = []
    if sub.street:
        location_parts.append(sub.street)
    if sub.city:
        location_parts.append(sub.city)

    loc_str = ", ".join(location_parts)
    if sub.state:
        loc_str = (loc_str + ", " if loc_str else "") + sub.state
    if sub.zip:
        loc_str = (loc_str + " " if loc_str else "") + sub.zip
    if loc_str:
        addr = addr + "; " + loc_str if addr else loc_str
    if sub.country:
        addr = addr + ", " + sub.country if addr else sub.country

    if addr:
        lines.extend(_wrap(addr, cont, cont))
    if sub.url:
        lines.extend(_wrap(f"URL:{sub.url}", cont, cont))

    return lines


def _reference_n_lines(ref_num: int, ref) -> list[str]:
    """Format a non-submission reference.

    JOURNAL format by status:
      Unpublished (no year)  : Unpublished.
      Unpublished (with year): Unpublished. (year)
      In press               : journal (year) In press
      Published              : journal volume, start_page[-end_page] (year)
    """
    lines = []
    lines.append(f"REFERENCE   {ref_num}  ")
    authors = _authors_str(ref.ab_names)
    lines.extend(_wrap(authors, " " * 12, "  AUTHORS   "))
    lines.extend(_wrap(ref.title, " " * 12, "  TITLE     "))

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

    lines.extend(_wrap(journal_str, " " * 12, "  JOURNAL   "))

    if ref.doi:
        lines.extend(_wrap(f"DOI:{ref.doi}", " " * 12, "  REMARK    "))
    elif ref.pubmed:
        lines.extend(_wrap(f"PUBMED:{ref.pubmed}", " " * 12, "  REMARK    "))

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
        mapped = map_tagset_id(st.tagset_id)
        header = f"##{mapped}-START##"
        footer = f"##{mapped}-END##"

        if lines:
            lines.append(cont)  # blank separator line
        else:
            prefix = cont

        if not lines:
            lines.extend(_wrap(header, cont, "COMMENT     "))
        else:
            lines.append(f"{cont}{header}")

        for key, val in st.fields:
            lines.append(f"{cont}{key:<22}:: {val}")

        lines.append(f"{cont}{footer}")

    return lines


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
) -> str:
    """Build one DDBJ flat file record and return it as a string."""
    today = date.today()
    if submission_date is None:
        submission_date = today
    if file_date is None:
        file_date = today

    seq_len = len(seq)
    entry_id = entry.entry_id
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

    # LOCUS
    wl(_locus_line(entry_id, seq_len, entry.topology, division, file_date))

    # DEFINITION (period at end)
    def_str = definition if definition.endswith(".") else definition + "."
    for line in _definition_lines(def_str):
        wl(line)

    # ACCESSION / VERSION
    wl(f"ACCESSION   {entry_id}")
    wl(f"VERSION     {entry_id}.1")

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
) -> None:
    """Write DDBJ flat file to *output* (file object or path string)."""
    close_after = False
    if isinstance(output, (str, __import__("pathlib").Path)):
        output = open(output, "w", encoding="utf-8")
        close_after = True

    try:
        for entry in entries:
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

            record_str = build_entry(
                entry=entry,
                common=common,
                seq=seq,
                division=division,
                submission_date=submission_date,
                file_date=file_date,
                email=email,
                no_taxonomy=no_taxonomy,
            )
            output.write(record_str)
    finally:
        if close_after:
            output.close()
