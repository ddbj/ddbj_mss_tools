"""Assemble a DDBJ MSS .ann: features from gff2mss, COMMON/source from common."""
from __future__ import annotations

from Bio import SeqIO

from ddbj_gff import parse
from ddbj_gff.io import open_text
from gff2mss.config import load_config
from gff2mss.convert import build_entry_features
from gff2mss.emit import feature_rows
from gff2mss.product_map import load_product_map

from common.models import load_common_json
from common.common_builder import create_common
from common.source_builder import load_sequence_roles, source_qualifier, ff_definition
from common.gap_annotator import GapAnnotator, annotate_gaps
from common.submission_category import inject_defaults, validate_and_fill


def build_ann_text(gff_path, fasta_path, mss_config_path, common_path,
                   sequence_roles_path, submission_category, locus_tag_start=None):
    cfg, _ = load_config(mss_config_path)
    if cfg.product_map_path:
        cfg.product_map = load_product_map(cfg.product_map_path)
    if locus_tag_start is not None:
        cfg.locus_tag_start = locus_tag_start  # continue numbering across compartments (genome-unique)

    with open_text(fasta_path) as fh:
        seqs = {rec.id: rec.seq for rec in SeqIO.parse(fh, "fasta")}
    with open_text(gff_path) as fh:
        doc = parse(fh.read())

    diagnostics: list = []
    per_entry = build_entry_features(doc, seqs, cfg, diagnostics)

    common = load_common_json(common_path)
    common_dict = common.model_dump(exclude_none=True)
    # category may be given on CLI or inside the JSON (_submission_category)
    category = submission_category or common.model_extra.get("_submission_category", "")
    if category:
        common_dict["_submission_category"] = category
        inject_defaults(common_dict, category)      # DATATYPE / DIVISION / KEYWORD defaults
        validate_and_fill(common_dict, category)    # required source/DBLINK/ST_COMMENT
    roles = load_sequence_roles(sequence_roles_path) if sequence_roles_path else {}

    base_source = dict(common.SOURCE or {})
    src_id_key = common.SOURCE_IDENTIFIER
    mol_type = base_source.get("mol_type", "")

    # every sequence in the input FASTA gets an entry (FASTA order), whether or not it
    # has annotation; unannotated sequences get a source-only entry.
    all_ids = list(seqs.keys())
    is_wgs = all((roles.get(e) is None or roles.get(e).type == "unplaced") for e in all_ids)
    segment_count = sum(1 for e in roles.values() if e.type == "segment")
    chromosome_count = sum(1 for e in roles.values() if e.type == "chromosome")

    gap_annotators: list = []
    gap_cfg = common.ASSEMBLY_GAP
    if gap_cfg:
        cfgs = gap_cfg if isinstance(gap_cfg, list) else [gap_cfg]
        gap_annotators = [GapAnnotator(linkage_evidence=c.linkage_evidence,
                                       min_gap_length=c.min_gap_length,
                                       max_gap_length=c.max_gap_length,
                                       gap_type=c.gap_type,
                                       estimated_length=c.estimated_length)
                          for c in cfgs if c.enabled]

    rows: list = list(create_common(common_dict))
    for entry_id in all_ids:
        seq = seqs[entry_id]
        length = len(seq)
        role = roles.get(entry_id)
        is_circular = role.is_circular if role is not None else False
        if is_circular:
            rows.append([entry_id, "TOPOLOGY", "", "circular", ""])
        src = dict(base_source)
        src.update(source_qualifier(role, entry_id, is_wgs, segment_count=segment_count))
        src["ff_definition"] = ff_definition(role, src_id_key, mol_type, is_wgs,
                                             chromosome_count=chromosome_count, segment_count=segment_count)
        items = list(src.items())
        first_col = "" if is_circular else entry_id
        rows.append([first_col, "source", f"1..{length}", items[0][0], str(items[0][1])])
        for k, v in items[1:]:
            rows.append(["", "", "", k, str(v)])
        for feat in per_entry.get(entry_id, []):   # features if annotated, else source only
            rows.extend(feature_rows(feat))
        if gap_annotators:
            rows.extend(annotate_gaps(gap_annotators, str(seq)))

    ann_text = "\n".join("\t".join(r) for r in rows) + "\n"
    # output every input sequence (annotation or not), in FASTA order
    out_seqs = dict(seqs)
    return ann_text, out_seqs
