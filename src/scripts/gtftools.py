#!/usr/scripts/env python
import sys


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
                    i = item.find(" ")
                    assert i != -1
                    key = item[:i]
                    value = item[i + 2:-1]
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


class TranscriptTable(object):
    pass


class TranscriptBed(object):
    pass


class GtfTools(object):
    @classmethod
    def table(cls, args):
        obj = Table(args[0], args[1])
        obj.execute()


def main():
    usage = """Usage:
    gtftools.py command [args]

Commands:
    table       Information table.
"""

    if len(sys.argv) < 2:
        print(usage)
        exit(1)

    command = sys.argv[1]
    args = sys.argv[2:]

    if command == "table":
        GtfTools.table(args)
    else:
        print(usage)
        exit(1)


if __name__ == "__main__":
    main()
