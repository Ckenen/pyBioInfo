import os
from pyBioInfo.Range import MRange
from pyBioInfo.Utils import SegmentPairBuilder
from .bam import BamFile, Alignment

class Fragment(MRange):
    def __init__(self, segment1, segment2):
        # assert segment1.query_name == segment2.query_name
        mate1 = Alignment(segment1)
        mate2 = Alignment(segment2)
        super(Fragment, self).__init__(chrom=mate1.chrom,
                                       name=mate1.name,
                                       strand=mate1.strand,
                                       blocks_array=[mate1.blocks,
                                                     mate2.blocks])
        self._mate1 = mate1
        self._mate2 = mate2

    @property
    def mate1(self):
        return self._mate1

    @property
    def mate2(self):
        return self._mate2


class FamFile(BamFile):
    def __init__(self, path, mode="rb", template=None, header=None):
        assert mode == "rb" or mode == "wb"
        super(FamFile, self).__init__(path, mode, template, header)

    def fetch(self, chrom=None, start=None, end=None):
        if chrom is None:
            last = None
            for i, segment in enumerate(self._handle.fetch(until_eof=True)):
                if i % 2 == 1:
                    assert last.query_name == segment.query_name
                    yield Fragment(last, segment)
                last = segment
        else:
            raise RuntimeError("Random access is not supported! Please try FamFileRandom.")

    def write(self, obj):
        if isinstance(obj, Fragment):
            self._handle.write(obj.mate1.segment)
            self._handle.write(obj.mate2.segment)
        else:
            super(FamFile, self).write(obj)


class FamFileRandom(FamFile):
    def __init__(self, path):
        assert path.endswith(".bam")
        assert os.path.exists(path)
        assert os.path.exists(path + ".bai")
        super(FamFileRandom, self).__init__(path, "rb")

    @property
    def handle(self):
        return self._handle
    
    @property
    def references(self):
        return self._references

    def fetch(self, chrom=None, start=None, end=None):
        if chrom is None:
            for c in sorted(self._references.keys()):
                segments = self._handle.fetch(contig=c)
                for pair in SegmentPairBuilder(segments):
                    yield Fragment(pair.mate1, pair.mate2)
        else:
            segments = self._handle.fetch(contig=chrom, start=start, stop=end)
            for pair in SegmentPairBuilder(segments):
                yield Fragment(pair.mate1, pair.mate2)
