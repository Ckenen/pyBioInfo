#!/usr/scripts/env python
import os
import unittest
import filecmp
from pyBioInfo.Range import IRange, GRange
from pyBioInfo.IO.File import FastaFile

DIR = os.path.dirname(__file__)

class TestFastaFile(unittest.TestCase):
    def test_fasta_file(self):
        # Read and write
        path1 = os.path.join(DIR, "data/test.fasta.gz")
        path2 = os.path.join(DIR, "data/test.out.fasta")
        path3 = os.path.join(DIR, "data/test.out.fasta.ref")
        with FastaFile(path1, "r") as f, FastaFile(path2, "w") as fw:
            for record in f:
                fw.write(record)
        self.assertTrue(filecmp.cmp(path2, path3))

        # Random access
        with FastaFile(path2, random=True) as f:
            seq = list(f.fetch("seq1", 40, 50))[0].sequence
            self.assertEqual("CAATACTTTC", seq)

            seq = list(f.fetch("seq1", 50, 60))[0].sequence
            self.assertEqual("TACCAGCTAT", seq)

            seq = list(f.fetch("seq1", 45, 55))[0].sequence
            self.assertEqual("CTTTCTACCA", seq)

            obj = GRange(chrom="seq1", blocks=[(40, 45), (55, 60)], strand="-")
            seq = list(f.fetch(obj=obj, strandness=False))[0].sequence
            self.assertEqual("CAATAGCTAT", seq)
            seq = list(f.fetch(obj=obj))[0].sequence
            self.assertEqual("ATAGCTATTG", seq)
            
        # Fast random access
        with FastaFile(path2, random=True) as f:
            seq = f.fast_fetch("seq1", 40, 50)
            self.assertEqual("CAATACTTTC", seq)

            seq = f.fast_fetch("seq1", 50, 60)
            self.assertEqual("TACCAGCTAT", seq)

            seq = f.fast_fetch("seq1", 45, 55)
            self.assertEqual("CTTTCTACCA", seq)

            obj = GRange(chrom="seq1", blocks=[(40, 45), (55, 60)], strand="-")
            seq = f.fast_fetch(obj=obj, strandness=False)
            self.assertEqual("CAATAGCTAT", seq)
            seq = f.fast_fetch(obj=obj)
            self.assertEqual("ATAGCTATTG", seq)


if __name__ == '__main__':
    unittest.main()
