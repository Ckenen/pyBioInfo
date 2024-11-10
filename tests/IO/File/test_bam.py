#!/usr/scripts/env python

import os
import filecmp
import unittest
from pyBioInfo.IO.File import BamFile

DIR = os.path.dirname(__file__)

class TestBamFile(unittest.TestCase):
    def test_bam_file(self):
        # Read bam
        records = []
        path1 = os.path.join(DIR, "data/test.bam")
        with BamFile(path1) as f:
            for record in f:
                # print(record.format("BED"))
                records.append(record)
    
        # Write to BED
        path2 = os.path.join(DIR, "data/test.bam.bed")
        with open(path2, "w+") as fw:
            for record in records:
                fw.write(record.format("BED") + "\n")
    
        path3 = os.path.join(DIR, "data/test.bam.bed.ref")
        self.assertTrue(filecmp.cmp(path2, path3))
        
        # Write to BAM
        path4 = os.path.join(DIR, "data/test.bam.out.bam")
        with BamFile(path1) as f, BamFile(path4, "wb", template=f) as fw:
            for obj in f:
                fw.write(obj)
                
        path5 = os.path.join(DIR, "data/test.bam.out.bam.bed")
        with BamFile(path4) as f, open(path5, "w+") as fw:
            for obj in f:
                fw.write(obj.format("BED") + "\n")
        self.assertTrue(filecmp.cmp(path3, path5))

        """
        Random access
        """
        

if __name__ == '__main__':
    unittest.main()
