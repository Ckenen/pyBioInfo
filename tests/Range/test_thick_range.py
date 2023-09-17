import unittest
from pyBioInfo.Range import TRange


class TestThickRange(unittest.TestCase):
    def test_thick_range(self):
        blocks = [(10, 20), (30, 40)]
        
        obj = TRange(chrom="chr1", blocks=blocks)
        self.assertEqual((0, None, None, 20), obj.indexes())

        obj = TRange(chrom="chr1", blocks=blocks, thick=(15, 35))
        self.assertEqual((0, 5, 15, 20), obj.indexes())

        obj = TRange(chrom="chr1", blocks=blocks, thick=(12, 35), strand="-")
        self.assertEqual((0, 5, 18, 20), obj.indexes())
        self.assertEqual((0, 2, 15, 20), obj.indexes(strandness=False))

if __name__ == '__main__':
    unittest.main()
