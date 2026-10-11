#!/usr/bin/env python
import argparse
from typing import cast
from pyBioInfo.Tools.fastqtools import FastqTools

def main():
    parser = argparse.ArgumentParser(description="A toolkit to process FAM file")

    subparsers = parser.add_subparsers(title="Available subcommands",dest="command")

    parser2 = subparsers.add_parser("fastq-to-txt", help="Convert FASTQ file to TXT file", description="Convert FASTQ file to TXT file")
    _ = parser2.add_argument("-1", "--input1", metavar="PATH", help="Input FASTQ file. (stdin, text or gzip)")
    _ = parser2.add_argument("-2", "--input2", metavar="PATH", help="Input FASTQ file of read 2 (text or gzip)")
    _ = parser2.add_argument("-o", "--output", metavar="PATH", help="Output TXT file (stdout, text or gzip)")

    parser2 = subparsers.add_parser("dedup-fastq", help="Dedup FASTQ file", description="Dedup FASTQ file")
    _ = parser2.add_argument("-1", "--input1", metavar="PATH", help="Input FASTQ file. (text or gzip)")
    _ = parser2.add_argument("-2", "--input2", metavar="PATH", help="Input FASTQ file of read 2. (text or gzip)")
    _ = parser2.add_argument("-o", "--output1", metavar="PATH", help="Output FASTQ file. (text or gzip)")
    _ = parser2.add_argument("-p", "--output2", metavar="PATH", help="Output FASTQ file of read2. (text or gzip)")
    _ = parser2.add_argument("-t", "--threads", type=int, default=4, help="Number of threads")

    parser2 = subparsers.add_parser("txt-to-fastq", help="Convert TXT file to FASTQ file", description="Convert TXT file to FASTQ file")
    _ = parser2.add_argument("-i", "--input", metavar="PATH", help="Input TXT file. (stdin, text or gzip)")
    _ = parser2.add_argument("-o", "--output1", metavar="PATH", help="Output FASTQ file. (text or gzip)")
    _ = parser2.add_argument("-p", "--output2", metavar="PATH", help="Output FASTQ file of read2. (text or gzip)")
    _ = parser2.add_argument("-t", "--threads", type=int, default=4, help="Number of threads")

    args = parser.parse_args()
    command = cast(str | None, args.command)
    if command is None:
        parser.print_help()
        return
    dispatch = {
        "fastq-to-txt": FastqTools.fastq_to_txt,
        "dedup-fastq": FastqTools.dedup_fastq,
        "txt-to-fastq": FastqTools.txt_to_fastq,
    }
    dispatch[command](args)

if __name__ == "__main__":
    main()