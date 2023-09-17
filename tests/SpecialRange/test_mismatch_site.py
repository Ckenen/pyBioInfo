#!/usr/scripts/env python

import unittest
from pyBioInfo.IO.File import BamFile, FamFile
# from pyBioInfo.SpecialRange import MismatchEvent, MismatchEventFactory


# class TestMismatchSite(unittest.TestCase):
#     def test_mismatch_site(self):
#         infile = "/home/chenzonggui/1_cwork/AMUC_seq/results/mapping/coincided/20191029_PDP0_Rep1/chr19.fam"
#         with FamFile(infile) as f:
#             for frag in f:
#                 for site in MismatchEventFactory.from_fragment(frag):
#                     print(site.chrom, site.start, site.end, site.ref, site.alt, site.strand)
#                 break
#                 # print(align.format("BED"))


# if __name__ == '__main__':
#     unittest.main()
