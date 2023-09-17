#!/usr/scripts/env python
import os
import filecmp
import unittest
from pyBioInfo.IO.File import VcfFile

DIR = os.path.dirname(__file__)


class TestVcfFile(unittest.TestCase):
    def test_vcf_file(self):
        path1 = os.path.join(DIR, "data/test.vcf.gz")
        path2 = os.path.join(DIR, "data/test.temp.vcf")
        path3 = os.path.join(DIR, "data/test.vcf.ref")
        # with VcfFile(path1) as f, VcfFile(path2, "w", template=f) as fw:
        #     for record in f:
        #         fw.write(record)
        # self.assertTrue(filecmp.cmp(path2, path3))


if __name__ == '__main__':
    unittest.main()
