#!/usr/bin/env python
import os
import unittest
import pysam
# from pyBioInfo.Common import BuildFragment

DIR = os.path.dirname(__file__)


class TestBuildFragment(unittest.TestCase):
    def check_fragment_file(self, path, segments):
        pass
        # segments1 = []
        # with pysam.AlignmentFile(path) as handle:
        #     for segment in handle:
        #         segments1.append(segment)
        # self.assertEqual(len(segments), len(segments1))
        # pairs = []
        # for i in range(0, len(segments1), 2):
        #     mate1 = segments1[i]
        #     mate2 = segments1[i + 1]
        #     self.assertEqual(mate1.reference_name, mate2.reference_name)
        #     self.assertEqual(mate1.query_name, mate2.query_name)
        #     self.assertTrue(mate1.is_read1)
        #     self.assertTrue(mate2.is_read2)
        #     chrom = mate1.reference_name
        #     start = min(mate1.reference_start, mate2.reference_start)
        #     end = max(mate1.reference_end, mate2.reference_end)
        #     pair = [chrom, start, end]
        #     pairs.append(pair)
        # last = None
        # for pair in pairs:
        #     if last:
        #         self.assertGreaterEqual(pair[0], last[0])
        #         if pair[0] == last[0]:
        #             self.assertGreaterEqual(pair[1], last[1])
        #             if pair[1] == last[1]:
        #                 self.assertGreaterEqual(pair[2], last[2])
        #     last = pair

    def test_build_fragment(self):
        pass
        # path1 = os.path.join(DIR, "data/example.bam")
        # segments = []
        # with pysam.AlignmentFile(path1) as handle:
        #     for segment in handle:
        #         segments.append(segment)

        # path2 = os.path.join(DIR, "data/example.build.single.fam")
        # BuildFragment.build_fragment(bam=path1, fam=path2, processor=1)
        # self.check_fragment_file(path2, segments)

        # path3 = os.path.join(DIR, "data/example.build.multiple.fam")
        # BuildFragment.build_fragment(bam=path1, fam=path3, processor=2, capacity=200)
        # self.check_fragment_file(path3, segments)


if __name__ == '__main__':
    unittest.main()
