#!/usr/scripts/env python

import sys

script, path = sys.argv

# path = "/home/chenzonggui/cluster_data1/wangyafen/20190313_wangyafen_fC/data/snp150Common.txt"
# opath = "/home/chenzonggui/cluster_data1/wangyafen/20190313_wangyafen_fC/data/snp150Common.tbl"

BASES = list("ACTGN")
header = "chrom start end name score strand rbase tbase".split()
print("\t".join(header))
with open(path) as f:
    for line in f:
        row = line.strip("\n").split("\t")
        if row[11] != "single":
            continue

        chrom = row[1]
        start = int(row[2])
        end = start + 1
        strand = row[6]
        rbases = set(row[7:9])

        for rbase in rbases:
            if rbase not in BASES:
                continue
            for tbase in row[9].split("/"):
                if tbase in rbases:
                    continue
                if tbase not in BASES:
                    continue

                s = "\t".join(map(str, [chrom, start, end, row[4], ".", strand, rbase, tbase]))
                print(s)
