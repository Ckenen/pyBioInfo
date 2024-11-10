#!/usr/scripts/env python
import os
import filecmp
import unittest
from pyBioInfo.IO.File import BamFile, FamFile

DIR = os.path.dirname(__file__)

class TestFamFile(unittest.TestCase):
    def test_fam_file(self):    
        path1 = os.path.join(DIR, "data/test.bam")
        path2 = os.path.join(DIR, "data/test.fam")
        with FamFile(path1, random=True) as f, FamFile(path2, "wb", template=f) as fw:
            for record in f:
                fw.write(record)

        path3 = os.path.join(DIR, "data/test.fam.bed")
        with FamFile(path2) as f, open(path3, "w+") as fw:
            for record in f:
                fw.write(record.format("BED") + "\n")

        path4 = os.path.join(DIR, "data/test.fam.bed.ref")
        filecmp.cmp(path3, path4)


if __name__ == '__main__':
    unittest.main()
