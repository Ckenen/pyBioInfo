import argparse
import os
import shutil
import subprocess
from typing import cast
import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def run_cmd(cmd: str) -> int:
    logging.info("RUN CMD: %s" % cmd)
    return subprocess.check_call(cmd, shell=True)

class BedTools(object):
    @classmethod
    def bed_to_bigwig(cls, args: argparse.Namespace):
        bedfile = cast(str, args.input)
        outfile1 = cast(str, args.output)
        outfile2 = cast(str, args.output2)
        scale = cast(float, args.scale)
        stranded = cast(bool, args.stranded)
        genome_size = cast(str, args.genome_size)
        
        tmpdir = outfile1 + ".bedtools_bed_to_bigwig_tmp"
        os.makedirs(tmpdir, exist_ok=True)

        if stranded:
            bg1 = os.path.join(tmpdir, "converted.+.bg")
            bg2 = os.path.join(tmpdir, "converted.-.bg")
            _ = run_cmd(f"bedtools genomecov -bg -i {bedfile} -g {genome_size} -scale {scale} -strand + > {bg1}")
            _ = run_cmd(f"bedGraphToBigWig {bg1} {genome_size} {outfile1}")
            _ = run_cmd(f"bedtools genomecov -bg -i {bedfile} -g {genome_size} -scale {scale} -strand - > {bg2}")
            _ = run_cmd(f"bedGraphToBigWig {bg2} {genome_size} {outfile2}")

        else:
            bg1 = os.path.join(tmpdir, "converted.bg")
            _ = run_cmd(f"bedtools genomecov -bg -i {bedfile} -g {genome_size} -scale {scale} > {bg1}")
            _ = run_cmd(f"bedGraphToBigWig {bg1} {genome_size} {outfile1}")

        shutil.rmtree(tmpdir)

        logging.info("Done")
        