import os
import pysam
from .file import BaseFile
from pyBioInfo.Range import GRange
from pyBioInfo.Utils import SegmentTools


# class SegmentTools(object):
#     """
#     SegmentTools contains many functions to manipulate Segment instance (from pysam package).
#     """
#     @classmethod
#     def get_block_from_segment(cls, segment, fill_gap=True):
#         blocks = []
#         bs = None  # block start
#         be = None  # block end
#         le = segment.reference_start  # last end
#         for m, n in segment.cigartuples:
#             if m == pysam.CMATCH:  # M
#                 bs, be = le, le + n
#                 le = be
#                 if len(blocks) > 0 and bs == blocks[-1][1]:  # concat
#                     blocks[-1][1] = be
#                 else:
#                     blocks.append([bs, be])
#             elif m == pysam.CDEL:  # D
#                 bs, be = le, le + n
#                 le = be
#                 if fill_gap:
#                     if len(blocks) > 0 and bs == blocks[-1][1]:
#                         blocks[-1][1] = be
#                     else:
#                         blocks.append([bs, be])
#             elif m == pysam.CREF_SKIP:  # N
#                 le = le + n
#             elif m == pysam.CINS:  # I
#                 continue
#             elif m == pysam.CSOFT_CLIP:  # S
#                 continue
#             elif m == pysam.CHARD_CLIP:  # H
#                 continue
#             else:
#                 raise RuntimeError("Unknown cigar type: %s" % m)
#         return blocks

#     @classmethod
#     def get_aligned_sequence(cls, segment):
#         # 把deletion的地方用横杠和-1的质量补上，insertion的地方去除掉
#         query_sequence = segment.query_sequence
#         query_qualities = segment.query_qualities
#         aligned_query_sequence = []
#         aligned_query_qualities = []

#         offset = 0
#         for flag, count in segment.cigartuples:
#             if flag == pysam.CHARD_CLIP:
#                 continue
#             elif flag == pysam.CSOFT_CLIP:
#                 offset = offset + count
#             elif flag == pysam.CREF_SKIP:
#                 continue
#             elif flag == pysam.CDEL:
#                 aligned_query_sequence.append("-" * count)
#                 aligned_query_qualities.extend([-1] * count)
#             elif flag == pysam.CINS:
#                 offset += count
#             elif flag == pysam.CMATCH:
#                 tmp = offset + count
#                 aligned_query_sequence.append(query_sequence[offset:tmp])
#                 aligned_query_qualities.extend(query_qualities[offset:tmp])
#                 offset = tmp
#             else:
#                 raise RuntimeError(
#                     "Unknown cigar type: %s [%d]" % (flag, count))
#         aligned_query_sequence = "".join(aligned_query_sequence)
#         return aligned_query_sequence, aligned_query_qualities

#     @classmethod
#     def merge_deletion_gap(cls, items):
#         items = items.copy()
#         while True:
#             found = False
#             i = 0
#             for item in items:
#                 if item[0] == pysam.CDEL:
#                     found = True
#                     break
#                 i += 1
#             if found:
#                 if i == 0:
#                     if i == len(items) - 1:
#                         raise ValueError("Only one DEL")
#                     else:
#                         item = items[i]
#                         items.pop(i)
#                         if items[i][0] == pysam.CMATCH and items[i][1] == item[2]:
#                             items[i][1] = item[1]
#                 else:
#                     if i == len(items) - 1:
#                         item = items[i]
#                         items.pop(i)
#                         if items[-1][0] == pysam.CMATCH and items[-1][2] == item[1]:
#                             items[-1][2] = item[2]
#                     else:
#                         item = items[i]
#                         items.pop(i)
#                         if items[i - 1][0] == pysam.CMATCH:
#                             if items[i][0] == pysam.CMATCH:
#                                 if items[i - 1][2] == item[1]:
#                                     if items[i][1] == item[2]:
#                                         item1 = items[i]
#                                         items.pop(i)
#                                         items[i - 1][2] = item1[2]
#                                     else:
#                                         items[i - 1][2] = item[2]
#                                 else:
#                                     if items[i][1] == item[2]:
#                                         items[i][1] = item[1]
#                                     else:
#                                         continue
#                             else:
#                                 if items[i - 1][2] == item[1]:
#                                     items[i - 1][2] = item[2]
#                         else:
#                             if items[i][0] == pysam.CMATCH:
#                                 if items[i][1] == item[2]:
#                                     items[i][2] = item[1]
#                             else:
#                                 continue
#             else:
#                 break
#         return items

#     @classmethod
#     def merge_insertion_gap(cls, items):
#         items = items.copy()
#         return items

#     @classmethod
#     def parse_cigar(cls, segment):
#         """
#         return: [
#             (
#                 operation, 
#                 (mapped_index_start, mapped_index_end), 
#                 (read_index_start, read_index_end), 
#                 (ref_position_start, ref_position_end)
#             ),
#             ...
#         ]
#         """
#         results = []
#         rps = segment.reference_start  # ref_position_start
#         rpe = 0  # ref_position_end
#         mis = 0  # mapped_index_start
#         mie = 0  # mapped_index_end
#         ris = 0  # read_index_start
#         rie = 0  # read_index_end
#         for operation, length in segment.cigartuples:
#             if operation == pysam.CMATCH:
#                 rpe = rps + length
#                 mie = mis + length
#                 rie = ris + length
#                 results.append(("M", (mis, mie), (ris, rie), (rps, rpe)))
#                 rps = rpe
#                 mis = mie
#                 ris = rie
#             elif operation == pysam.CINS:
#                 rie = ris + length
#                 results.append(("I", (mis, mis), (ris, rie), (rps, rps)))
#                 ris = rie
#             elif operation == pysam.CDEL:
#                 rpe = rps + length
#                 mie = mis + length
#                 results.append(("D", (mis, mie), (ris, ris), (rps, rpe)))
#                 rps = rpe
#                 mis = mie
#             elif operation == pysam.CREF_SKIP:
#                 rpe = rps + length
#                 results.append(("N", (mis, mis), (ris, ris), (rps, rps)))
#                 rps = rpe
#             elif operation == pysam.CSOFT_CLIP:
#                 rie = ris + length
#                 results.append(("S", (mis, mis), (ris, rie), (rps, rps)))
#                 ris = rie
#             elif operation == pysam.CHARD_CLIP:
#                 results.append(("H", (mis, mis), (ris, ris), (rps, rps)))
#                 # raise RuntimeError()
#             elif operation == pysam.CPAD:
#                 raise RuntimeError()
#             elif operation == pysam.CEQUAL:
#                 raise RuntimeError()
#             elif operation == pysam.CDIFF:
#                 raise RuntimeError()
#             elif operation == pysam.CBACK:
#                 raise RuntimeError()
#             else:
#                 raise RuntimeError()
#         return results

#     @classmethod
#     def parse_md_tag(cls, segment):
#         """_summary_

#         Args:
#             segment (_type_): _description_

#         Raises:
#             RuntimeError: _description_

#         Returns:
#             _type_: _description_
#         """
#         results = []
#         mis = 0  # mapped_index_start
#         mie = 0  # mapped_index_end
#         md = segment.get_tag("MD")
#         i = 0
#         while i < len(md):
#             if "0" <= md[i] <= "9":
#                 j = i + 1
#                 while j < len(md) and "0" <= md[j] <= "9":
#                     j += 1
#                 mie = mis + int(md[i:j])
#                 results.append(("=", (mis, mie), None))
#                 mis = mie
#                 i = j
#             elif "A" <= md[i] <= "Z":
#                 j = i + 1
#                 while j < len(md) and "A" <= md[j] <= "Z":
#                     j += 1
#                 mie = mis + j - i
#                 results.append(("X", (mis, mie), md[i:j]))
#                 mis = mie
#                 i = j
#             elif md[i] == "^":
#                 j = i + 1
#                 while j < len(md) and "A" <= md[j] <= "Z":
#                     j += 1
#                 mie = mis + j - i - 1
#                 results.append(("D", (mis, mie), md[i + 1:j]))
#                 mis = mie
#                 i = j
#             else:
#                 raise RuntimeError()
#         return results

#     @classmethod
#     def fetch_snvs(cls, segment=None, parsed_cigar=None, parsed_md_tag=None):
#         """[summary]

#         Args:
#             segment ([type], optional): [description]. Defaults to None.
#             parsed_cigar ([type], optional): [description]. Defaults to None.
#             parsed_md_tag ([type], optional): [description]. Defaults to None.

#         Raises:
#             RuntimeError: [description]
#             RuntimeError: [description]
#             RuntimeError: [description]
#             RuntimeError: [description]
#             RuntimeError: [description]
#         """
#         if parsed_cigar is None:
#             parsed_cigar = cls.parse_cigar(segment)
#         if parsed_md_tag is None:
#             parsed_md_tag = cls.parse_md_tag(segment)
#         array1 = parsed_cigar
#         array2 = parsed_md_tag
#         sequence = segment.query_sequence
#         qualities = segment.query_qualities
#         results = []
#         for item1 in array1:
#             if item1[0] == "I":
#                 alt_seq = sequence[item1[2][0]:item1[2][1]]
#                 scores = qualities[item1[2][0]:item1[2][1]]
#                 results.append(("I",
#                                item1[1], item1[2], item1[3],
#                                None, alt_seq, scores))
#         sequence = segment.query_sequence
#         qualities = segment.query_qualities
#         i1 = 0
#         i2 = 0
#         while i2 < len(array2):
#             item2 = array2[i2]
#             if item2[0] == "=":
#                 i2 += 1
#                 continue
#             elif item2[0] == "X":
#                 assert item2[1][1] - item2[1][0] == 1
#                 found = False
#                 while i1 < len(array1):
#                     item1 = array1[i1]
#                     if item1[1][1] <= item2[1][0]:
#                         i1 += 1
#                         continue
#                     elif item1[1][0] >= item2[1][1]:
#                         raise RuntimeError()
#                     else:
#                         if item1[0] == "M":
#                             offset = item2[1][0] - item1[1][0]
#                             length = item2[1][1] - item2[1][0]
#                             ris = item1[2][0] + offset
#                             rie = ris + length
#                             rps = item1[3][0] + offset
#                             rpe = rps + length
#                             alt_seq = sequence[ris:rie]
#                             scores = qualities[ris:rie]
#                             results.append(("X",
#                                            item2[1], (ris, rie), (rps, rpe),
#                                            item2[2], alt_seq, scores))
#                             found = True
#                             break
#                         elif item1[0] == "I":
#                             i1 += 1
#                             continue
#                         else:
#                             raise RuntimeError()
#                 assert found
#                 i2 += 1
#             elif item2[0] == "D":
#                 found = False
#                 while i1 < len(array1):
#                     item1 = array1[i1]
#                     if item1[1][1] <= item2[1][0]:
#                         i1 += 1
#                         continue
#                     elif item1[1][0] >= item2[1][1]:
#                         raise RuntimeError()
#                     else:
#                         if item1[0] == "D":
#                             assert item1[1][0] == item2[1][0]
#                             assert item1[1][1] == item2[1][1]
#                             results.append(("D",
#                                            item2[1], item1[2], item1[3],
#                                            item2[2], None, None))
#                             i1 += 1
#                             found = True
#                             break
#                         else:
#                             raise RuntimeError()
#                 assert found
#                 i2 += 1
#             else:
#                 raise RuntimeError()
#         results = list(sorted(results, key=lambda item: item[1]))
#         return results


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
        parsed_cigar = self.get_parsed_cigar()
        query_sequence = self.segment.query_sequence
        base = None
        for block in parsed_cigar:
            if block[3][0] <= position < block[3][1]:
                if block[0] == "M":
                    offset = position - block[3][0]
                    read_idx = block[2][0] + offset
                    base = query_sequence[read_idx]
                elif block[0] == "D":
                    base = "-"
                else:
                    assert False
                break
        return base
    

class BamFile(BaseFile):
    def __init__(self, path, mode="rb", template=None, header=None, require_index=True):
        if mode == "rb":
            if require_index:
                index_path = path + ".bai"
                if not os.path.exists(index_path):
                    raise IOError("Index file %s is not exists!" % index_path)
        elif mode == "wb":
            if template is None and header is None:
                raise RuntimeError(
                    "Arguments template and header can not both be None.")
        else:
            raise RuntimeError(
                "Unsupported mode: %s, only rb and wb are supported!" % mode)

        self._template = template
        self._header = header
        self._require_index = require_index
        self._references = None
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
                self._handle = pysam.AlignmentFile(
                    self._path, self._mode, header=header)
                self._header = header
            else:
                raise RuntimeError()

            references = dict()
            for item in self.handle.header["SQ"]:
                references[item["SN"]] = item["LN"]
            self._references = references

    def close(self):
        if self._handle:
            self._handle.close()
            self._handle = None

    def __iter__(self):
        for obj in self.fetch():
            yield obj

    def fetch(self, chrom=None, start=None, end=None):
        assert self._mode == "rb"
        if chrom is None:
            if self._require_index:
                for chrom in sorted(self._references.keys()):
                    for segment in self._handle.fetch(contig=chrom):
                        yield Alignment(segment)
            else:
                for segment in self._handle.fetch(until_eof=True):
                    if segment.is_unmapped:
                        continue
                    yield Alignment(segment)
        else:
            if self._require_index:
                for segment in self._handle.fetch(contig=chrom, start=start, stop=end):
                    yield Alignment(segment)
            else:
                raise RuntimeError("Can not random access!")

    def write(self, obj):
        assert self._mode == "wb"
        if isinstance(obj, Alignment):
            self._handle.write(obj.segment)
        elif isinstance(obj, pysam.AlignedSegment):
            self._handle.write(obj)
        else:
            raise TypeError("Unsupported record type: %s!" % type(obj))
