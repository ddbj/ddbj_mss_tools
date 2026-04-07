"""
Parse NCBI feature table files (.tbl) and format locations for DDBJ MSS.

Feature table format reference:
    https://www.ncbi.nlm.nih.gov/WebSub/html/help/feature-table.html

A feature table file consists of:
    >Feature lcl|ENTRY_ID
    start  end  feature_type
    start  end                ← additional exon of the same feature
        ...
                qualifier_name  qualifier_value
                flag_qualifier           ← no value

Coordinates are 1-based.
start > end  ←→  feature is on the minus strand.
Partial features: '<' prefix on start means 5'-partial;
                  '>' prefix on end   means 3'-partial.
"""

from __future__ import annotations


# ── Qualifier ordering and renaming ───────────────────────────────────────────

# Qualifiers to emit for each feature type, in order.
# Qualifiers absent from this mapping are silently dropped.
_QUAL_ORDER: dict[str, list[str]] = {
    "gene":          ["locus_tag", "gene", "gene_desc", "pseudo", "db_xref"],
    "mRNA":          ["locus_tag", "gene", "product",   "note",   "pseudo", "db_xref"],
    "CDS":           ["locus_tag", "gene", "product",   "codon_start", "transl_table",
                      "note", "pseudo", "exception", "db_xref"],
    "ncRNA":         ["locus_tag", "gene", "product", "ncRNA_class", "note", "db_xref"],
    "tRNA":          ["locus_tag", "gene", "product",  "note", "db_xref"],
    "rRNA":          ["locus_tag", "gene", "product",  "note", "db_xref"],
    "misc_RNA":      ["locus_tag", "gene", "product",  "note", "db_xref"],
    "precursor_RNA": ["locus_tag", "gene", "product",  "note"],
}

# Rename NCBI-specific qualifier names to DDBJ-acceptable ones.
# A value of None means the qualifier is dropped entirely.
_QUAL_RENAME: dict[str, str | None] = {
    "gene_desc":     "note",   # gene description → note
    "transcript_id": None,     # NCBI internal → drop
    "protein_id":    None,     # NCBI internal → drop
    "experiment":    None,     # evidence tag  → drop
}

# Default qualifier values added to CDS features when absent from the source.
_CDS_DEFAULTS = {"codon_start": "1", "transl_table": "1"}


# ── TBL parser ────────────────────────────────────────────────────────────────

def parse_tbl(tbl_path: str) -> dict[str, list[dict]]:
    """
    Parse an NCBI feature table file.

    Returns
    -------
    dict mapping entry_id → list of feature dicts, each with:
        "type"       : str
        "intervals"  : list of (start_str, end_str)  raw coordinate strings
        "qualifiers" : dict[str, list[str | None]]   qualifier name → values
                       (None value = flag qualifier with no value)
    """
    entries: dict[str, list[dict]] = {}
    current_entry: str | None = None
    features: list[dict] = []

    cur_type: str | None = None
    cur_intervals: list[tuple[str, str]] = []
    cur_quals: dict[str, list] = {}

    def _save_current() -> None:
        if cur_type is not None:
            features.append({
                "type": cur_type,
                "intervals": list(cur_intervals),
                "qualifiers": dict(cur_quals),
            })

    with open(tbl_path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")

            # ── New entry header ──────────────────────────────────────────
            if line.startswith(">Feature "):
                _save_current()
                if current_entry is not None:
                    entries[current_entry] = list(features)

                current_entry = line[len(">Feature "):].strip()
                if current_entry.startswith("lcl|"):
                    current_entry = current_entry[4:]
                features = []
                cur_type = None
                cur_intervals = []
                cur_quals = {}
                continue

            if current_entry is None:
                continue

            # ── Qualifier line (3 leading tabs) ──────────────────────────
            if line.startswith("\t\t\t"):
                rest = line[3:]
                parts = rest.split("\t", 1)
                key = parts[0]
                val = parts[1] if len(parts) > 1 else None
                cur_quals.setdefault(key, []).append(val)
                continue

            # ── Coordinate / feature-type line ────────────────────────────
            cols = line.split("\t")
            if len(cols) >= 2 and cols[0].lstrip("<>-").isdigit():
                start_str = cols[0].strip()
                end_str = cols[1].strip()
                feat_type = cols[2].strip() if len(cols) >= 3 and cols[2].strip() else None

                if feat_type:
                    _save_current()
                    cur_type = feat_type
                    cur_intervals = [(start_str, end_str)]
                    cur_quals = {}
                elif cur_type is not None:
                    cur_intervals.append((start_str, end_str))

    _save_current()
    if current_entry is not None:
        entries[current_entry] = list(features)

    return entries


# ── Location formatting ───────────────────────────────────────────────────────

def _strip_partial(s: str) -> int:
    """Return integer coordinate from a string that may start with < or >."""
    return int(s.lstrip("<>"))


def _is_minus_strand(intervals: list[tuple[str, str]]) -> bool:
    """True if any interval has start > end (NCBI tbl minus-strand convention)."""
    return any(_strip_partial(s) > _strip_partial(e) for s, e in intervals)


def format_location(intervals: list[tuple[str, str]]) -> str:
    """
    Convert NCBI feature-table intervals to a DDBJ MSS location string.

    NCBI tbl convention:
        Plus strand  → start < end  (e.g. 100  200)
        Minus strand → start > end  (e.g. 200  100)
        Partial 5'   → '<' prefix on the start of the 5'-most exon
        Partial 3'   → '>' prefix on the end   of the 3'-most exon

    DDBJ MSS convention:
        join(s1..e1,s2..e2,...)
        complement(join(s1..e1,...))   ← exons in ascending genomic order
        <N  (5'-partial) and >N (3'-partial) inside the location
    """
    minus = _is_minus_strand(intervals)

    if minus:
        processed: list[tuple[int, str, str]] = []
        for s_raw, e_raw in intervals:
            s_int = _strip_partial(s_raw)
            e_int = _strip_partial(e_raw)
            lo, hi = min(s_int, e_int), max(s_int, e_int)
            # On minus strand: s_raw is the 3' end, e_raw is the 5' end
            lo_str = f"<{lo}" if e_raw.startswith("<") else str(lo)
            hi_str = f">{hi}" if s_raw.startswith(">") else str(hi)
            processed.append((lo, lo_str, hi_str))

        processed.sort(key=lambda x: x[0])
        exon_strs = [f"{lo_s}..{hi_s}" for _, lo_s, hi_s in processed]
    else:
        exon_strs = [f"{s}..{e}" for s, e in intervals]

    loc = exon_strs[0] if len(exon_strs) == 1 else f"join({','.join(exon_strs)})"
    return f"complement({loc})" if minus else loc


# ── Qualifier collection ──────────────────────────────────────────────────────

def collect_qualifiers(feat: dict) -> list[tuple[str, str | None]]:
    """
    Return an ordered list of (qualifier_name, value) pairs for a feature,
    after applying renaming and filtering rules.
    """
    ftype = feat["type"]
    raw_quals = feat["qualifiers"]
    allowed = _QUAL_ORDER.get(ftype, list(raw_quals.keys()))

    result: list[tuple[str, str | None]] = []
    seen: set[str] = set()

    for qname in allowed:
        dest_name = _QUAL_RENAME.get(qname, qname)
        if dest_name is None:
            continue

        vals = raw_quals.get(qname, [])
        if not vals:
            if ftype == "CDS" and qname in _CDS_DEFAULTS:
                result.append((qname, _CDS_DEFAULTS[qname]))
                seen.add(qname)
            continue

        for v in vals:
            key = f"{dest_name}:{v}"
            if key not in seen:
                result.append((dest_name, v))
                seen.add(key)

    if ftype == "CDS":
        for qname, default in _CDS_DEFAULTS.items():
            if not any(n == qname for n, _ in result):
                result.append((qname, default))

    return result
