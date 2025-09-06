#!/usr/bin/env python
import optparse
import pandas as pd 

def main():
    parser = optparse.OptionParser("%prog transcript_info.csv canonical.tsv transcript_info.added_canonical.csv")
    
    _, args = parser.parse_args()
    infile1, infile2, outfile = args
    
    m = pd.read_csv(infile1)
    
    canonicals = set()
    d = pd.read_csv(infile2, sep="\t")
    canonicals = set([tid.split(".")[0] for tid in d["transcript"]])

    final_canonicals = []
    for gid, tmp in m.groupby(by="GeneID"):
        tmp = tmp.sort_values(by="Length")
        tids1 = tmp["TranscriptID"].values
        tids2 = list(filter(lambda tid: tid.split(".")[0] in canonicals, tids1))
        if len(tids2) == 0:
            final_canonicals.append(tids1[-1])
        else:
            final_canonicals.append(tids2[-1])
    final_canonicals = set(final_canonicals)
    m["Canonical"] = [tid in final_canonicals for tid in m["TranscriptID"]]
    m.to_csv(outfile, sep="\t" if outfile.endswith(".tsv") else ",", index=False)


if __name__ == "__main__":
    main()
    
