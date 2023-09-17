from pyBioInfo.Range import MRange
from pyBioInfo.Utils import SegmentPairBuilder
from .bam import BamFile, Alignment


# class SegmentPair(object):
#     """
#     An object that contains segment1 and segment2.
#     """

#     def __init__(self, mate1, mate2):
#         # assert mate1.reference_name == mate2.reference_name
#         # assert mate1.is_read1
#         # assert mate2.is_read2
#         self.chrom = mate1.reference_name
#         self.start = min(mate1.reference_start, mate2.reference_start)
#         self.end = max(mate1.reference_end, mate2.reference_end)
#         self.mate1 = mate1
#         self.mate2 = mate2

#     def __lt__(self, other):
#         if self.chrom < other.chrom:
#             return True
#         elif self.chrom == other.chrom:
#             if self.start < other.start:
#                 return True
#             elif self.start == other.start:
#                 if self.end < other.end:
#                     return True
#         return False


# class SegmentPairBuilder(object):
#     """
#     Create SegmentPaired objects.
#     """

#     def __init__(self, segments):
#         self.segments = segments

#     @classmethod
#     def _build_pair(cls, segments):
#         # [s1, s2, s3, s4, ...]
#         array1 = []
#         # [[s1, s2, s3, s4, ...], [s5, s6, s7, s8, ...]]
#         array2 = []
#         # [pair1, pair2, ...]
#         array3 = []

#         # step 1: grouped by query name
#         for segment in sorted(segments, key=lambda obj: obj.query_name):
#             if len(array1) == 0:
#                 array1.append(segment)
#             else:
#                 if segment.query_name == array1[0].query_name:
#                     array1.append(segment)
#                 else:
#                     array2.append(array1)
#                     array1 = [segment]
#         if len(array1) > 0:
#             array2.append(array1)

#         # step 2: build pair
#         for array1 in array2:
#             array1 = list(sorted(array1, key=lambda obj: obj.reference_start))
#             while len(array1) >= 2:
#                 s1 = array1.pop(0)  # segment 1
#                 s2 = None  # segment 2
#                 idx = 0
#                 for i, tmp in enumerate(array1):
#                     if tmp.reference_start < s1.next_reference_start:
#                         continue
#                     elif tmp.reference_start == s1.next_reference_start:
#                         if tmp.next_reference_start == s1.reference_start:
#                             idx = i
#                             s2 = tmp
#                             break
#                     else:
#                         break
#                 if s2:
#                     array1.pop(idx)
#                     if s1.is_read2:
#                         s1, s2 = s2, s1
#                     pair = SegmentPair(s1, s2)
#                     array3.append(pair)

#         # step 3: sorting pairs
#         array3 = list(sorted(array3))
#         for item in array3:
#             yield item

#     def __iter__(self):
#         chrom, require_start = None, None
#         count = None
#         array = None
#         last = None
#         for segment in self.segments:
#             assert segment.is_proper_pair
#             chrom1 = segment.reference_name
#             start1 = segment.reference_start
#             start2 = segment.next_reference_start
#             if chrom is None:
#                 chrom = chrom1
#                 require_start = max(start1, start2)
#                 array = [segment]
#                 count = 1
#             else:
#                 if chrom1 == chrom:
#                     assert start1 >= last.reference_start
#                     if count >= 32 and start1 > require_start:
#                         for pair in self._build_pair(array):
#                             yield pair
#                         require_start = max(start1, start2)
#                         array = [segment]
#                         count = 1
#                     else:
#                         require_start = max(require_start, start1, start2)
#                         array.append(segment)
#                         count += 1
#                 elif chrom1 > chrom:
#                     for pair in self._build_pair(array):
#                         yield pair
#                     chrom = chrom1
#                     require_start = max(start1, start2)
#                     array = [segment]
#                     count = 1
#                 else:
#                     raise RuntimeError()
#             last = segment
#         if chrom:
#             for pair in self._build_pair(array):
#                 yield pair


class Fragment(MRange):
    def __init__(self, segment1, segment2):
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
        if path.endswith(".fam"):
            require_index = False
        else:
            require_index = True
        super(FamFile, self).__init__(path=path,
                                      mode=mode,
                                      template=template,
                                      header=header,
                                      require_index=require_index)

    def fetch(self, chrom=None, start=None, end=None):
        assert self._mode == "rb"
        if chrom is None:
            if self._require_index:  # build
                for chrom in sorted(self._references.keys()):
                    for pair in SegmentPairBuilder(self._handle.fetch(contig=chrom)):
                        yield Fragment(pair.mate1, pair.mate2)
            else:  # not build
                last = None
                for i, segment in enumerate(self._handle.fetch(until_eof=True)):
                    if i % 2 == 1:
                        yield Fragment(last, segment)
                    last = segment
        else:
            if self._require_index:  # build
                for pair in SegmentPairBuilder(self._handle.fetch(contig=chrom,
                                                                  start=start,
                                                                  stop=end)):
                    yield Fragment(pair.mate1, pair.mate2)
            else:
                raise RuntimeError("Can not random access!")

    def write(self, obj):
        if isinstance(obj, Fragment):
            self._handle.write(obj.mate1.segment)
            self._handle.write(obj.mate2.segment)
        else:
            super().write(obj)
