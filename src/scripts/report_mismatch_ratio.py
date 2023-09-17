#!/usr/bin/env python
import sys, optparse
from collections import defaultdict
import pandas as pd
from pyBioInfo.IO.File import BamFile, FamFile


MAPPER = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "N": "N",
    "-": "-"}


def make_table(d1, d2, name):
    rows = []
    for k, v in d1.items():
        r = d2[k[0]]
        p = v / r
        rows.append([k, v, r, p])
    dat = pd.DataFrame(rows)
    dat.index = dat[0]
    dat = dat[dat.columns[1:]]
    dat.columns = ["Count", "Reference", "Ratio"]
    dat.index.name = "Type"
    dat = dat[["-" not in x for x in dat.index]]
    dat = dat[["N" not in x for x in dat.index]]
    dat = dat.sort_index()
    s = dat["Ratio"]
    s.name = name
    return s


def process_se(infiles):
    raise NotImplementedError()


def process_pe(infiles):
    counter1 = defaultdict(int)
    counter2 = defaultdict(int)
    counter3 = defaultdict(int)
    total_fragments = 0
    data = dict()
    for a in "ACGT":
        for b in "ACGT":
            if a != b:
                mtype = "%s%s" % (a, b)
                data[mtype] = [defaultdict(int), defaultdict(int)]
    for infile in infiles:
        with FamFile(infile) as f:
            for n, frag in enumerate(f):
                # if n > 50000:
                #     break
                strand = frag.mate1.segment.get_tag("ST")
                try:
                    pe = frag.mate1.segment.get_tag("CE")
                except KeyError:
                    pe = frag.mate1.segment.get_tag("PE")
                events = []
                for event in pe.split(";"):
                    if event != "":
                        event = event.split(",")
                        ref_base = event[1]
                        alt_base = event[2]
                        quality = int(event[3])
                        if quality < 20:
                            continue
                        if strand == "-":
                            ref_base = MAPPER[ref_base]
                            alt_base = MAPPER[alt_base]
                        mtype = "%s%s" % (ref_base, alt_base)
                        events.append(mtype)

                for event in events:
                    counter1[event] += 1

                pc = frag.mate1.segment.get_tag("PC")
                counts = []
                for item in pc.split(";"):
                    if item != "":
                        base, count = item.split(",")
                        count = int(count)
                        if strand == "-":
                            base = MAPPER[base]
                        counts.append([base, count])
                for base, count in counts:
                    counter2[base] += count

                total_fragments += 1
                for mtype in data.keys():
                    if mtype in events:
                        counter3[mtype] += 1
                        d1, d2 = data[mtype]
                        for event in events:
                            d1[event] += 1
                        for base, count in counts:
                            d2[base] += count

    array = [make_table(counter1, counter2, "Total")]
    for mtype in data.keys():
        c1, c2 = data[mtype]
        array.append(make_table(c1, c2, mtype))

    dat = pd.concat(array, axis=1, sort=False)
    dat = dat.fillna(0)
    dat.index.name = "Type"
    columns = list(dat.columns)
    for i in range(len(columns)):
        if columns[i] == "Total":
            columns[i] = "%s (%d)" % (columns[i], total_fragments)
        else:
            columns[i] = "%s (%d)" % (columns[i], counter3[columns[i]])
    dat.columns = columns
    return dat



def main():
    infiles = sys.argv[1:1]
    outfile = sys.argv[-1]
    assert len(infiles) >= 1
    
    if infiles[0].endswith(".bam"):
        dat = process_se(infiles)
    elif infiles[0].endswith(".fam"):
        dat = process_pe(infiles)
    else:
        raise RuntimeError()

    dat.to_csv(outfile, sep="\t")