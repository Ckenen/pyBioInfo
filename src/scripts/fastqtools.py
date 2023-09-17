#!/usr/scripts/env python

import sys
from optparse import OptionParser

class FastqUniq(object):
    def __init__(self, args):
        pass

class SubCommands(object):
    @classmethod
    def uniq(cls, args):
        usage = "\n".join([
        "fastqtools.py fpkm input.fastq [input2.fastq] output.fastq [output.fastq]"
        ])
        parser = OptionParser(usage=usage)
        options, args = parser.parse_args(args)

        if len(args) != 2:
            parser.parse_args(["-h"])
            exit(1)

    @classmethod
    def uniq_single_end(cls, infile, outfile):
        raise NotImplementedError()

    @classmethod
    def uniq_paired_end(cls, infile1, infile2, outfile1, outfile2):
        pass


def main():
    lines = [
        "fastqtools.py command [args]",
        "",
        "Available command:",
        "\nuniq",
        ""
    ]
    usage = "\n".join(lines)

    if len(sys.argv) < 2:
        print(usage)
        exit(1)

    command = sys.argv[1]
    args = sys.argv[2:]
    # print(args)

    if command == "uniq":
        SubCommands.uniq(args)
    else:
        print(usage)
        exit(1)


if __name__ == "__main__":
    main()