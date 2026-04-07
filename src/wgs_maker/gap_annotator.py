from dataclasses import dataclass
import re

# constant values for assembly_gap feature
# see https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/


@dataclass
class GapAnnotator:

    min_gap_length: int = 10  # Default 10, Recommendation from the curation staff
    linkage_evidence: str = "paired-ends"  # default "paired-ends"  # paired-ends, align genus, proximity ligation, etc. 
    gap_type: str = "within scaffolds"
    gap_length:str = "known"

    gap_pattern: re.Pattern|None = None # re.compile('n{%d,}' % min_gap_length)


    @staticmethod
    def initialize(args):

        min_gap_length = args.min_gap_length  # Default 10, Recommendation from the curation staff
        linkage_evidence = args.linkage_evidence.replace("_", " ")  # default "paired-ends"  # paired-ends, align genus, proximity ligation, etc. 

        if args.gap_type == "auto":
            if linkage_evidence == "paired-ends":
                gap_type = "within scaffolds"
            elif linkage_evidence == "proximity ligation":
                gap_type = "within scaffolds"
            elif linkage_evidence == "align genus":
                gap_type = "within scaffolds"
            else:
                raise AssertionError("Please specify gap_type")
        else:
            gap_type = args.gap_type.replace("_", " ")  # within scaffold, between scaffolds, telomere, centromere, short arm, heterochromatin, repeat_within_scaffold, repeat between_scaffolds, contamination, unknown

        if args.gap_length == "auto":
            if linkage_evidence == "paired-ends":
                gap_length = "known"
            elif linkage_evidence == "proximity ligation":
                gap_length = "unknown"
            elif linkage_evidence == "align genus":
                gap_length = "unknown"
            else:
                raise AssertionError("Please specify gap_length")
        else:
            gap_length = args.gap_length # known or unknown  
        
        gap_pattern = re.compile('n{%d,}' % min_gap_length)

        return GapAnnotator(min_gap_length, linkage_evidence, gap_type, gap_length, gap_pattern)

    
    def create_gap_feature(self, seq, seq_name=None):
        ret = []
        # length_value = "known" if gap_length else "unknown"

        for match in self.gap_pattern.finditer(seq.lower()):
            start = match.start() + 1  # +1 because INSDC coordinate is 1-based
            end = match.end()
            location = f"{start}..{end}"
            ret.append(["", "assembly_gap", location, "estimated_length", self.gap_length])
            ret.append(["", "", "", "gap_type", self.gap_type])
            ret.append(["", "", "", "linkage_evidence", self.linkage_evidence])
        if ret and seq_name:  # For draft genome, (souce feature is described in COMMON)
            ret[0][0] = seq_name
        return ret

if __name__ == "__main__":
    config = GapAnnotator()
    print(config)
