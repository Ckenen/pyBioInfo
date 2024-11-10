import os
import logging
import pysam
from .base import BaseFile
from pyBioInfo.Range import GRange
from pyBioInfo.Utils import SegmentTools


class Alignment(GRange):
    def __init__(self, segment):
        chrom = segment.reference_name
        name = segment.query_name
        strand = "-" if segment.is_reverse else "+"
        blocks = SegmentTools.get_block_from_segment(segment)
        super(Alignment, self).__init__(chrom=chrom,
                                        name=name,
                                        blocks=blocks,
                                        strand=strand)
        self._segment = segment
        self._aligned_sequences = None
        self._aligned_qualities = None
        self._parsed_cigar = None
        self._parsed_md_tag = None

    @property
    def segment(self):
        return self._segment

    def _init_aligned_sequence_and_qualities(self):
        self._aligned_sequences, self._aligned_qualities = SegmentTools.get_aligned_sequence(
            self.segment)

    def get_aligned_sequences(self):
        if self._aligned_sequences is None:
            self._init_aligned_sequence_and_qualities()
        return self._aligned_sequences, self._aligned_qualities
    
    def get_parsed_cigar(self):
        if self._parsed_cigar is None:
            self._parsed_cigar = SegmentTools.parse_cigar(self.segment)
        return self._parsed_cigar
    
    def get_pared_md_tag(self):
        if self._parsed_md_tag is None:
            self._parsed_md_tag = SegmentTools.parse_md_tag(self.segment)
        return self._parsed_md_tag

    def get_query_base(self, position):
        return SegmentTools.get_query_base(self.segment, position, self.get_parsed_cigar())
    

class BamFile(BaseFile):
    def __init__(self, path, mode="rb", template=None, header=None, random=None):
        if mode[0] == "r" and random is None:
            if path.endswith(".bam"):
                if os.path.exists(path + ".bai"):
                    logging.warning("Index file exists, automatic set random=True")
                    random = True
                else:
                    logging.warning("Index file does not exists, automatic set random=False")
                    random = False
            else:
                random = False
        
        if random:
            assert path.endswith(".bam")
            assert os.path.exists(path)
            assert os.path.exists(path + ".bai")
        assert mode == "rb" or mode == "wb"
        self._template = template
        self._header = header
        self._references = None
        self._random = random
        super(BamFile, self).__init__(path, mode)
        self.open()

    @property
    def header(self):
        return self._header

    @property
    def references(self):
        return self._references

    def open(self):
        if self._handle is None:
            if self.mode == "rb":
                self._handle = pysam.AlignmentFile(self.path)
                self._header = self._handle.header.as_dict().copy()
            elif self.mode == "wb":
                header = self._header
                if header is None:
                    if self._template is None:
                        raise RuntimeError()
                    else:
                        if isinstance(self._template, BamFile):
                            header = self._template.handle.header.as_dict().copy()
                        elif isinstance(self._template, pysam.AlignmentFile):
                            header = self._template.header.as_dict().copy()
                        else:
                            raise RuntimeError()
                assert header
                self._handle = pysam.AlignmentFile(self._path, "wb", header=header)
                self._header = header
            else:
                raise RuntimeError()
            self._references = dict()
            for item in self._handle.header["SQ"]:
                self._references[item["SN"]] = item["LN"]

    def close(self):
        if self._handle:
            self._handle.close()
            self._handle = None

    def fetch(self, chrom=None, start=None, end=None):
        if self._random:
            if chrom is None:
                for c in list(sorted(self._references)):
                    for s in self._handle.fetch(contig=c):
                        if s.is_unmapped:
                            continue
                        yield Alignment(s)
            else:
                for s in self._handle.fetch(contig=chrom, start=start, stop=end):
                    yield Alignment(s)
        else:
            if chrom is None:
                for segment in self._handle.fetch(until_eof=True):
                    if segment.is_unmapped:
                        continue
                    yield Alignment(segment)
            else:
                raise RuntimeError("Random access is not supported! Please set random=True.")

    def write(self, obj):
        if isinstance(obj, Alignment):
            self._handle.write(obj.segment)
        elif isinstance(obj, pysam.AlignedSegment):
            self._handle.write(obj)
        else:
            raise TypeError("Unsupported record type: %s!" % type(obj))

# class BamFileRandom(BamFile):
#     def __init__(self, path):
#         assert path.endswith(".bam")
#         assert os.path.exists(path)
#         assert os.path.exists(path + ".bai")
#         super(BamFileRandom, self).__init__(path, "rb")

#     def fetch(self, chrom=None, start=None, end=None):
#         if chrom is None:
#             for c in list(sorted(self._references)):
#                 for s in self._handle.fetch(contig=c):
#                     if s.is_unmapped:
#                         continue
#                     yield Alignment(s)
#         else:
#             for s in self._handle.fetch(contig=chrom, start=start, stop=end):
#                 yield Alignment(s)