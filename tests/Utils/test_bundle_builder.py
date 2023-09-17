#!/usr/scripts/env python

import unittest
from pyBioInfo.Range import CRange
from pyBioInfo.Utils import BundleBuilder


class TestBundleBuilder(unittest.TestCase):
    def test_bundle_builder(self):
        array = [
            CRange(chrom="chr1", start=10, end=20),
            CRange(chrom="chr1", start=15, end=25),
            CRange(chrom="chr1", start=20, end=30),
            CRange(chrom="chr1", start=25, end=35),
            CRange(chrom="chr1", start=50, end=60),
            CRange(chrom="chr1", start=80, end=90),
            CRange(chrom="chr2", start=10, end=20),
            CRange(chrom="chr2", start=50, end=60),
            CRange(chrom="chr2", start=70, end=90),
            CRange(chrom="chr2", start=80, end=90),
            CRange(chrom="chr2", start=85, end=90),
            CRange(chrom="chr2", start=90, end=100),
            CRange(chrom="chr3", start=10, end=20),
            CRange(chrom="chr3", start=10, end=20),
            CRange(chrom="chr3", start=10, end=20),
            CRange(chrom="chr3", start=10, end=20),
        ]

        results = [
            ["chr1", 10, 11, 1],
            ["chr1", 15, 16, 1],
            ["chr1", 20, 21, 1],
            ["chr1", 25, 26, 1],
            ["chr1", 50, 51, 1],
            ["chr1", 80, 81, 1],
            ["chr2", 10, 11, 1],
            ["chr2", 50, 51, 1],
            ["chr2", 70, 71, 1],
            ["chr2", 80, 81, 1],
            ["chr2", 85, 86, 1],
            ["chr2", 90, 91, 1],
            ["chr3", 10, 11, 4]
        ]
        bundles = list(BundleBuilder(array, min_capacity=1, mode=BundleBuilder.MODE_START_ONLY))
        self.assertEqual(len(results), len(bundles))
        for result, bundle in zip(results, bundles):
            self.assertEqual(result, [bundle.chrom, bundle.start_min, bundle.start_max + 1, bundle.count])

        results = [
            ["chr1", 10, 16, 2],
            ["chr1", 20, 26, 2],
            ["chr1", 50, 81, 2],
            ["chr2", 10, 51, 2],
            ["chr2", 70, 81, 2],
            ["chr2", 85, 91, 2],
            ["chr3", 10, 11, 4]
        ]
        bundles = list(BundleBuilder(array, min_capacity=2, mode=BundleBuilder.MODE_START_ONLY))
        self.assertEqual(len(results), len(bundles))
        for result, bundle in zip(results, bundles):
            self.assertEqual(result, [bundle.chrom, bundle.start_min, bundle.start_max + 1, bundle.count])

        results = [
            ['chr1', 10, 35, 4],
            ['chr1', 50, 60, 1],
            ['chr1', 80, 90, 1],
            ['chr2', 10, 20, 1],
            ['chr2', 50, 60, 1],
            ['chr2', 70, 90, 3],
            ['chr2', 90, 100, 1],
            ['chr3', 10, 20, 4],
        ]
        bundles = list(BundleBuilder(array, min_capacity=1, mode=BundleBuilder.MODE_START_AND_END))
        self.assertEqual(len(results), len(bundles))
        for result, bundle in zip(results, bundles):
            self.assertEqual(result, [bundle.chrom, bundle.start_min, bundle.end_max, bundle.count])

        results = [
            ['chr1', 10, 35, 4],
            ['chr1', 50, 90, 2],
            ['chr2', 10, 60, 2],
            ['chr2', 70, 90, 3],
            ['chr2', 90, 100, 1],
            ['chr3', 10, 20, 4],
        ]
        bundles = list(BundleBuilder(array, min_capacity=2, mode=BundleBuilder.MODE_START_AND_END))
        self.assertEqual(len(results), len(bundles))
        for result, bundle in zip(results, bundles):
            self.assertEqual(result, [bundle.chrom, bundle.start_min, bundle.end_max, bundle.count])
    

if __name__ == '__main__':
    unittest.main()
