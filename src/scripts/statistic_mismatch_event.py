#!/usr/bin/env python
import sys
import os
import json
from optparse import OptionParser
from collections import defaultdict
import numpy as np
import pandas as pd
from pyBioInfo.IO.File import FamFile, VcfFile, FastaFile
from pyBioInfo.SpecialRange import MismatchEventFactory
from pyBioInfo.Utils import ShiftLoader

usage = """
%prog [options] input.fam genome.fasta outdir
"""

def load_snp(options):
    if options.snp:
        with VcfFile(options.snp) as vcf:
            for snp in vcf:
                yield snp


def load_fragment(fam, options):
    if options.fr:
        assert not options.rf
        for fragment in fam:
            fragment.strand = fragment.mate1.strand
            yield fragment
    elif options.rf:
        assert not options.fr
        for fragment in fam:
            fragment.strand = fragment.mate2.strand
            yield fragment
    else:
        for fragment in fam:
            fragment.strand = "+"
            yield fragment


def main():
    parser = OptionParser(usage=usage)
    parser.add_option("--rf", action="store_true", dest="rf", default=False)
    parser.add_option("--fr", action="store_true", dest="fr", default=False)
    parser.add_option("-s", "--snp", dest="snp", default=None)
    parser.add_option("-q", "--min-quality",
                      dest="quality", type="int", default=0)
    parser.add_option("-d", "--min-distance",
                      dest="distance", type="int", default=0)
    (options, args) = parser.parse_args()
    print(options)
    print(args)

    path1, path2, outdir = args
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    references = defaultdict(int)
    bases = list("ACGTN")
    counter = dict()
    for b1 in bases:
        for b2 in bases:
            if b1 != b2:
                mtype = "%s-%s" % (b1, b2)
                counter[mtype] = dict()

    with FamFile(path1) as fam, FastaFile(path2) as fasta:
        sloader = ShiftLoader(load_snp(options))
        for i, fragment in enumerate(load_fragment(fam, options)):
            # if i > 100000:
            #     break
            # Counting reference base
            seqeunces = fasta.fetch(obj=fragment).upper()
            for base in seqeunces:
                references[base] += 1

            # Fetching mismatch event
            events = [e for e in MismatchEventFactory.from_fragment(fragment)]
            events = list(sorted(events))
            if fragment.strand == "-":
                for e in events:
                    e.reverse_complement()
                    e = fragment.strand

            # Mark SNPs
            loader1 = ShiftLoader(sloader.fetch(obj=fragment))
            for e in events:
                for s in loader1.fetch(obj=e):
                    if s.induce(e.rbase, e.tbase):
                        e.is_snp = True
                        break

            # Counting
            for e in events:
                mtype = "%s-%s" % (e.rbase, e.tbase)
                if mtype == "N-N":
                    continue
                dist = e.distance
                score = e.score
                dict1 = counter[mtype]
                if dist not in dict1.keys():
                    dict1[dist] = dict()
                dict2 = dict1[dist]
                if score not in dict2.keys():
                    dict2[score] = {"snp": 0, "clean": 0}
                if e.is_snp:
                    dict2[score]["snp"] += 1
                else:
                    dict2[score]["clean"] += 1

    bpath = os.path.join(outdir, "references.tsv")
    spath = os.path.join(outdir, "snp.stat.tsv")
    qpath = os.path.join(outdir, "score.stat.tsv")
    dpath = os.path.join(outdir, "distance.stat.tsv")
    mpath = os.path.join(outdir, "matrix.tsv")
    epath = os.path.join(outdir, "events.json")
    rpath1 = os.path.join(outdir, "ratio.raw.tsv")
    rpath2 = os.path.join(outdir, "ratio.filtered.tsv")

    with open(bpath, "w+") as fw:
        total = sum(references.values())
        fw.write("Base\tCount\tPercent\n")
        for base in sorted(references.keys()):
            count = references[base]
            percent = np.divide(count, total)
            fw.write("\t".join(map(str, [base, count, percent])))

    with open(spath, "w+") as fw:
        snp = 0
        clean = 0
        for mtype, d1 in counter.items():
            for distance, d2 in d1.items():
                for score, d3 in d2.items():
                    snp += d3["snp"]
                    clean += d3["clean"]
        total = snp + clean
        fw.write("Type\tCount\tPercent\n")
        fw.write("SNP\t%d\t%.6f\n" % (snp, snp / total))
        fw.write("Clean\t%d\t%.6f\n" % (clean, clean / total))

    counter1 = defaultdict(int)  # score
    counter2 = defaultdict(int)  # distance
    counter3 = defaultdict(int)  # score and distance
    counter4 = defaultdict(int)  # event (raw)
    counter5 = defaultdict(int)  # event (filtered)
    for mtype, d1 in counter.items():
        for distance, d2 in d1.items():
            for score, d3 in d2.items():
                count = d3["clean"]
                counter1[score] += count
                counter2[distance] += count
                counter3[(score, distance)] += count
                counter4[mtype] += count
                if distance >= options.distance and score >= options.quality:
                    counter5[mtype] += count

    with open(qpath, "w+") as fw:
        fw.write("Score\tCount\tPercent\n")
        total = sum(counter1.values())
        for score in sorted(counter1.keys()):
            count = counter1[score]
            percent = count / total
            fw.write("\t".join(map(str, [score, count, percent])) + "\n")

    with open(dpath, "w+") as fw:
        fw.write("Distance\tCount\tPercent\n")
        total = sum(counter2.values())
        for distance in sorted(counter2.keys()):
            count = counter2[distance]
            percent = count / total
            fw.write("\t".join(map(str, [distance, count, percent])) + "\n")

    scores = list(sorted(set([key[0] for key in counter3.keys()])))
    distances = list(sorted(set([key[1] for key in counter3.keys()])))
    matrix = np.zeros(len(scores) * len(distances),
                      dtype=np.int).reshape(len(distances), -1)
    for key, value in counter3.items():
        # print(key, value)
        cidx = scores.index(key[0])
        ridx = distances.index(key[1])
        matrix[ridx][cidx] = value
    dat = pd.DataFrame(index=pd.Index(
        distances, name="Distance"), columns=scores, data=matrix)
    dat.to_csv(mpath, sep="\t")

    with open(rpath1, "w+") as fw:
        fw.write("Type\tCount\tRef\tRatio\n")
        for mtype in sorted(counter4.keys()):
            count = counter4[mtype]
            ref = references[mtype[0]]
            ratio = count / ref
            fw.write("\t".join(map(str, [mtype, count, ref, ratio])) + "\n")

    with open(rpath2, "w+") as fw:
        fw.write("Type\tCount\tRef\tRatio\n")
        for mtype in sorted(counter5.keys()):
            count = counter5[mtype]
            ref = references[mtype[0]]
            ratio = count / ref
            fw.write("\t".join(map(str, [mtype, count, ref, ratio])) + "\n")


if __name__ == "__main__":
    main()
