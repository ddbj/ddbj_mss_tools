import gzip
from Bio import SeqIO
import re

def read_fasta(file_name):
    if file_name.endswith(".gz"):
        return list(SeqIO.parse(gzip.open(file_name, "rt"), "fasta"))
    else:
        return list(SeqIO.parse(file_name, "fasta"))


def check_number_of_seqs(seq_list,  dict_sequence):
    seq_names, seq_types, seq_topologies = dict_sequence.setdefault("seq_names", []), dict_sequence.setdefault("seq_types", []), dict_sequence.setdefault("seq_topologies", [])
    if not seq_names:
        for seq in seq_list:
            seq_names.append(seq.id)
        
    if not (len(seq_list) == len(seq_names) == len(seq_types) == len(seq_topologies)):
        raise AssertionError("The number of sequences is not consistent.")

# constant values for assembly_gap feature
# see https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/

MIN_GAP_LENGTH = 10  # Default 10, Recommendation from the curation staff
# LINKAGE_EVIDENCE = args.linkage_evidence.replace("_", "")  # default "paired-ends"  # paired-ends, align genus, proximity ligation, etc. 
LINKAGE_EVIDENCE = "paired-ends"  # default "paired-ends"  # paired-ends, align genus, proximity ligation, etc. 
GAP_TYPE = "within scaffolds"
GAP_LENGTH = "known"

# if gap_type == "auto":
#     if LINKAGE_EVIDENCE == "paired-ends":
#         GAP_TYPE = "within scaffolds"
#     elif LINKAGE_EVIDENCE == "proximity ligation":
#         GAP_TYPE = "within scaffolds"
#     else:
#         raise AssertionError("Please specify gap_type")
# else:
#     GAP_TYPE = args.gap_type.replace("_", " ")  # within scaffold, between scaffolds, telomere, centromere, short arm, heterochromatin, repeat_within_scaffold, repeat between_scaffolds, contamination, unknown
# GAP_TYPE = "within scaffold"  # see https://www.ddbj.nig.ac.jp/ddbj/qualifiers.html#gap_type
# if args.gap_length == "auto":
#     if LINKAGE_EVIDENCE == "paired-ends":
#         GAP_LENGTH = "known"
#     elif LINKAGE_EVIDENCE == "proximity ligation":
#         GAP_LENGTH = "unknown"
#     else:
#         raise AssertionError("Please specify gap_length")
# else:
#     GAP_LENGTH = args.gap_length # known or unknown 

gap_pattern = re.compile('n{%d,}' % MIN_GAP_LENGTH)

# logging.basicConfig(level=logging.DEBUG)
# logger = logging.getLogger(__name__)

def create_gap_feature(seq, seq_name=None):
    ret = []
    # length_value = "known" if GAP_LENGTH else "unknown"

    for match in gap_pattern.finditer(seq.lower()):
        start = match.start() + 1  # +1 because INSDC coordinate is 1-based
        end = match.end()
        location = f"{start}..{end}"
        ret.append(["", "assembly_gap", location, "estimated_length", GAP_LENGTH])
        ret.append(["", "", "", "gap_type", GAP_TYPE])
        ret.append(["", "", "", "linkage_evidence", LINKAGE_EVIDENCE])
    if ret and seq_name:
        ret[0][0] = seq_name
    return ret

def create_source_feature(_trad_submission_category, seq_name, seq_type, seq_topology, source_dict):

    submitter_seqid = None
    environmental_sample = False
    plasmid = False
    if _trad_submission_category == "GNM":
        mol_type = "genomic DNA"
        if seq_type in ["c", "complete"]:
            ff_definition = "@@[organism]@@ @@[strain]@@ DNA, complete genome"
        elif seq_type in ["n", "nearly complete", "nearly-complete"]:
            ff_definition = "@@[organism]@@ @@[strain]@@ DNA, nearly complete genome"
        elif seq_type in ["p", "plasmid"]:
            ff_definition = "@@[organism]@@ @@[strain]@@ plasmid @@[plasmid]@@ DNA, complete sequence"
            plasmid = True
        else:
            submitter_seqid = "@@[entry]@@"
            ff_definition = "@@[organism]@@ @@[strain]@@ DNA, @@[submitter_seqid]@@"
    elif _trad_submission_category == "MAG":
        mol_type = "genomic DNA"
        environmental_sample = True
        if seq_type in ["c", "complete"]:
            ff_definition = "@@[organism]@@ @@[isolate]@@ DNA, complete genome"
        elif seq_type in ["n", "nearly complete", "nearly-complete"]:
            ff_definition = "@@[organism]@@ @@[isolate]@@ DNA, nearly complete genome"
        elif seq_type in ["p", "plasmid"]:
            ff_definition = "@@[organism]@@ @@[isolate]@@ plasmid @@[plasmid]@@ DNA, complete sequence"
            plasmid = True
        else:
            submitter_seqid = "@@[entry]@@"
            ff_definition = "@@[organism]@@ @@[isolate]@@ DNA, @@[submitter_seqid]@@"
    elif _trad_submission_category == "WGS":
        mol_type = "genomic DNA"
        submitter_seqid = "@@[entry]@@"
        ff_definition = "@@[organism]@@ @@[strain]@@ DNA, @@[submitter_seqid]@@"
    elif _trad_submission_category == "MAG-WGS":
        mol_type = "genomic DNA"
        environmental_sample = True
        submitter_seqid = "@@[entry]@@"
        ff_definition = "@@[organism]@@ @@[isolate]@@ DNA, @@[submitter_seqid]@@"

    ret = []
    ret.append(["", "source", "1..E", "mol_type", mol_type])
    ret.append(["", "", "", "ff_definition", ff_definition])
    if submitter_seqid:
        ret.append(["", "", "", "submitter_seqid", submitter_seqid])
    if environmental_sample:
        ret.append(["", "", "", "environmental_sample", ""])
    if plasmid:
        ret.append(["", "", "", "plasmid", seq_name])
    for key, value in source_dict.items():
        ret.append(["", "", "", key, value])
    if _trad_submission_category in ["WGS", "MAG-WGS"]:
        # The source feature will be appended to COMMON. Nothing to do.
        pass
    elif _trad_submission_category in ["GNM", "MAG"]:
        if seq_topology in ["c", "circular"]:
            ret = [["", "TOPOLOGY", "", "circular", ""]] + ret
        ret[0][0] = seq_name
    return ret
