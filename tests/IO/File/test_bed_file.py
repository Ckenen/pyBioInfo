#!/usr/scripts/env python
import os
import unittest
import filecmp
from pyBioInfo.IO.File import BedFile

DIR = os.path.dirname(__file__)

class TestBedFile(unittest.TestCase):
    def test_bed_file(self):
        # Read
        records = []
        path1 = os.path.join(DIR, "data/test.bed")
        with BedFile(path1) as f:
            for record in f:
                records.append(record)

        # Write BED12
        path2 = os.path.join(DIR, "data/test.bed.out")
        with BedFile(path2, "w") as fw:
            for record in records:
                fw.write(record)
        self.assertTrue(filecmp.cmp(path1, path2))

        # Write BED6
        path3 = os.path.join(DIR, "data/test.bed6.out")
        with BedFile(path3, "w", ncol=6) as fw:
            for record in records:
                fw.write(record)
        path4 = os.path.join(DIR, "data/test.bed6.ref")
        self.assertTrue(filecmp.cmp(path3, path4))

        # Read and write BED12
        path1 = os.path.join(DIR, "data/test.bed12")
        path2 = os.path.join(DIR, "data/test.out.bed12")
        with BedFile(path1) as f, open(path2, "w+") as fw:
            for record in f:
                fw.write(record.format("bed", ncol=12) + "\n")
        self.assertTrue(filecmp.cmp(path1, path2))

        # Read and write BED3
        path1 = os.path.join(DIR, "data/test.bed3")
        path2 = os.path.join(DIR, "data/test.out.bed3")
        with BedFile(path1) as f, open(path2, "w+") as fw:
            for record in f:
                fw.write(record.format("bed", ncol=3) + "\n")
        self.assertTrue(filecmp.cmp(path1, path2))

        """
        Random access
        """
        path1 = os.path.join(DIR, "data/test.bed.gz")
        path2 = os.path.join(DIR, "data/test.bed.gz.out.gz")
        with BedFile(path1) as bed, BedFile(path2, "w") as fw:
            for item in bed.fetch("chr1", 14403, 29570):
                fw.write(item)
        names = [
            "ENST00000456328", "ENST00000515242", "ENST00000518655", "ENST00000541675", 
            "ENST00000423562", "ENST00000438504", "ENST00000488147"
        ]
        records = []
        with BedFile(path2) as bed:
            for record in bed:
                records.append(record)
        self.assertEqual(len(names), len(records))
        for name, record in zip(names, records):
            self.assertEqual(name, record.name)
        



if __name__ == '__main__':
    unittest.main()
