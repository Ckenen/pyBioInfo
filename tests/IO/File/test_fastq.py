#!/usr/scripts/env python
import os
import unittest
import filecmp
from pyBioInfo.IO.File import FastqFile

DIR = os.path.dirname(__file__)

class TestFastqFile(unittest.TestCase):
    def test_fastq_file(self):
        path1 = os.path.join(DIR, "data/test.fastq")
        records = []
        with FastqFile(path1) as f:
            for record in f:
                records.append(record)

        names = [
            "E00477:539:H5Y7TCCX2:5:1101:6928:1608",
            "E00477:539:H5Y7TCCX2:5:1101:10155:1608",
            "E00477:539:H5Y7TCCX2:5:1101:10236:1608",
            "E00477:539:H5Y7TCCX2:5:1101:10784:1608"
        ]

        sequences = [
            "NGCACAATCTCCGGGGGCAGATGAAGGTAATCACGGAGATACTGGATACCCTCATTGGTAAGGTACCAGTAGAAAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTGAAAAAAGTGAA",
            "NTTTGATCGTGGTGATTTAGAGGGTGAACTCACTGGAACGGGGATGCTTGCATGTGTAATCTTACTAAGAGCTAAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTGAAAAAAGTGTT",
            "NTCCCCGTTCCTCTGGCTGCTGGCCCCGGACACGAGACCCTCCGTCTCCTGCCGAGAGAGGCTCCCCCGGCCGCCTCGCCACGCCGCAGCCCCAGATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTT",
            "NGTCACAGGGTATCATCTTGGTGCCTTGGTGGATGAGCTGGTAAAAAGCATGCTGGCCATTGGTCCCTGGCTCCCCCAGATCGGAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTGAAAAAACAGT",
        ]

        self.assertEqual(len(names), len(records))

        for name, sequence, record in zip(names, sequences, records):
            self.assertEqual(name, record.name)
            self.assertEqual(sequence, record.sequence)

        # Write
        path2 = os.path.join(DIR, "data/test.out.fastq")
        path3 = os.path.join(DIR, "data/test.out.fastq.ref")
        with FastqFile(path2, "w") as fw:
            for record in records:
                fw.write(record)
        self.assertTrue(filecmp.cmp(path2, path3))


if __name__ == '__main__':
    unittest.main()
