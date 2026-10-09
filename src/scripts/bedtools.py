#!/usr/bin/env python
import argparse
from typing import cast
from pyBioInfo.Tools.bedtools import BedTools
        
def main():
    parser = argparse.ArgumentParser(description="A toolkit to process BED file")
    subparsers = parser.add_subparsers(title="Available subcommands", dest="command")

    bed2bw_parser = subparsers.add_parser("bed-to-bigwig", help="Convert BED to BIGWIG", description="Convert BED to BIGWIG")
    _ = bed2bw_parser.add_argument("-i", "--input", required=True, help="Input BED file")
    _ = bed2bw_parser.add_argument("-o", "--output", required=True, help="Output BIGWIG file")
    _ = bed2bw_parser.add_argument("-g", "--genome-size", required=True, help="Genome size file")
    _ = bed2bw_parser.add_argument("-p", "--output2", help="Output BIGWIG file (For stranded)")
    _ = bed2bw_parser.add_argument("-s", "--scale", type=float, default=1, help="Scale factor")
    _ = bed2bw_parser.add_argument("-S", "--stranded", type=bool, default=False, help="Seperate by strand (true or false). Default is false.")

    args = parser.parse_args()
    command = cast(str | None, args.command)
    if command is None:
        parser.print_help()
        return None
    dispatch = {
        "bed-to-bigwig": BedTools.bed_to_bigwig,
    }
    dispatch[command](args)
    
    

if __name__ == '__main__':
    main()