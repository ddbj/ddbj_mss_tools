"""Fetch taxonomy lineage from NCBI for a given organism name."""

from __future__ import annotations

import time

_cache: dict[str, tuple[str, str]] = {}

# Tagset-id mapping from MSS to DDBJ flat file representation
TAGSET_ID_MAP: dict[str, str] = {
    "Genome-Assembly-Data": "Assembly-Data",
}


def get_lineage(organism: str, email: str = "mss2ff@ddbj.nig.ac.jp") -> tuple[str, str]:
    """Return (taxon_id, lineage_str) for *organism*.

    lineage_str has the form used in DDBJ flat files:
        "Bacteria; Bacillati; Bacillota; Bacilli; Lactobacillales;
         Lactobacillaceae; Paucilactobacillus."
    Returns ("", "") if not found.
    """
    organism = organism.strip()
    if organism in _cache:
        return _cache[organism]

    try:
        from Bio import Entrez
        Entrez.email = email

        # --- search taxon ID ---
        handle = Entrez.esearch(db="taxonomy", term=f'"{organism}"[Scientific Name]')
        result = Entrez.read(handle)
        handle.close()
        time.sleep(0.34)

        if not result["IdList"]:
            # fallback: broader search
            handle = Entrez.esearch(db="taxonomy", term=organism)
            result = Entrez.read(handle)
            handle.close()
            time.sleep(0.34)

        if not result["IdList"]:
            _cache[organism] = ("", "")
            return "", ""

        taxon_id = result["IdList"][0]

        # --- fetch lineage ---
        handle = Entrez.efetch(db="taxonomy", id=taxon_id, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        time.sleep(0.34)

        if not records:
            _cache[organism] = (taxon_id, "")
            return taxon_id, ""

        record = records[0]
        lineage_raw = record.get("Lineage", "")

        # Remove "cellular organisms" and "root" prefixes
        parts = [p.strip() for p in lineage_raw.split(";")]
        parts = [p for p in parts if p and p not in ("cellular organisms", "root")]
        lineage_str = "; ".join(parts) + "." if parts else ""

        _cache[organism] = (taxon_id, lineage_str)
        return taxon_id, lineage_str

    except Exception as exc:  # noqa: BLE001
        import sys
        print(f"[taxonomy] Warning: could not fetch lineage for {organism!r}: {exc}", file=sys.stderr)
        _cache[organism] = ("", "")
        return "", ""


def map_tagset_id(tagset_id: str) -> str:
    """Map MSS tagset_id to DDBJ flat file ST_COMMENT header name."""
    return TAGSET_ID_MAP.get(tagset_id, tagset_id)
