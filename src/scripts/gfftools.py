#!/usr/scripts/env python
import sys
import os
import pandas as pd
from pyBioInfo.Range import GRange


class Table(object):
    def __init__(self, path_in, path_out):
        self.path_in = path_in
        self.path_out = path_out

    def execute(self):
        rows = []
        with open(self.path_in) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                columns = line.strip("\n").split("\t")
                columns[3] = str(int(columns[3]) - 1)
                attris = dict()
                for item in columns[8].split(";"):
                    item = item.strip()
                    if item == "":
                        continue
                    i = item.find("=")
                    assert i != -1
                    key = item[:i]
                    value = item[i + 1:]
                    attris[key] = value
                columns[8] = attris
                rows.append(columns)

        headers1 = ["chrom", "source", "feature", "start", "end", "score", "strand", "frame"]
        headers2 = set()
        for row in rows:
            for key in row[8].keys():
                headers2.add(key)
        headers2 = list(headers2)
        headers = headers1 + headers2
        with open(self.path_out, "w+") as fw:
            fw.write("\t".join(headers) + "\n")
            for row in rows:
                values1 = row[:8]
                values2 = [row[8].get(key, "") for key in headers2]
                values = values1 + values2
                fw.write("\t".join(values) + "\n")

class Annotation(object):
    def __init__(self, path_in, path_out):
        self.path_in = path_in
        self.path_out = path_out

    def execute(self):
        tpath = self.path_out + ".table.temp.tsv"
        obj = Table(self.path_in, tpath)
        obj.execute()
        table = pd.read_csv(tpath, sep="\t", low_memory=False)
        table = table[(table["feature"] == "exon") | (table["feature"] == "CDS")]
        assert False
        if os.path.exists(tpath):
            os.remove()

class Transcript(object):
    def __init__(self, path_in, path_out):
        self.path_in = path_in
        self.path_out = path_out

    def execute(self):
        tpath = self.path_out + ".table.temp.tsv"
        obj = Table(self.path_in, tpath)
        obj.execute()
        table = pd.read_csv(tpath, sep="\t")
        table = table[(table["feature"] == "exon") | (table["feature"] == "CDS")]
        transcripts = []
        for tid, dat in table.groupby(by="Parent"):
            exons = dat[dat["feature"] == "exon"]
            cdss = dat[dat["feature"] == "CDS"]
            if len(exons) == 0:
                continue
            blocks = []
            thick = None
            for x, y in exons[["start", "end"]].values:
                if x >= y:
                    print(x, y)
                    print(tid)
                blocks.append([x, y])
            if len(cdss) > 0:
                thicks = []
                for x, y in cdss[["start", "end"]].values:
                    thicks.append([x, y])
                thicks = list(sorted(thicks, key=lambda item: item[0]))
                tx = thicks[0][0]
                ty = thicks[-1][1]
                thick = [tx, ty]
                if len(blocks) == 0:
                    for x, y in cdss[["start", "end"]].values:
                        blocks.append([x, y])
            blocks = list(sorted(blocks, key=lambda item: item[0]))
            chroms = dat["chrom"].values
            strands = dat["strand"].values
            assert len(set(chroms)) == 1
            assert len(set(strands)) == 1
            chrom = chroms[0]
            strand = strands[0]
            obj = GRange(chrom=chrom, strand=strand, name=tid, blocks=blocks, thick=thick)
            transcripts.append(obj)
        with open(self.path_out, "w+") as fw:
            for transcript in sorted(transcripts):
                fw.write(transcript.format("BED") + "\n")
        if os.path.exists(tpath):
            os.remove(tpath)


class GffTools(object):
    @classmethod
    def table(cls, args):
        obj = Table(args[0], args[1])
        obj.execute()

    @classmethod
    def annotation(cls, args):
        obj = Annotation(args[0], args[1])
        obj.execute()

    @classmethod
    def transcript(cls, args):
        obj = Transcript(args[0], args[1])
        obj.execute()


def main():
    usage = """Usage:
    gfftools.py command [args]

Commands:
    table       Information table.
    transcript  Extract transcript.bed (sorted).
    annotation  Extract annotation table.
"""

    if len(sys.argv) < 2:
        print(usage)
        exit(1)

    command = sys.argv[1]
    args = sys.argv[2:]

    if command == "table":
        GffTools.table(args)
    elif command == "transcript":
        GffTools.transcript(args)
    elif command == "annotation":
        GffTools.annotation(args)
    else:
        print(usage)
        exit(1)



if __name__ == "__main__":
    main()
