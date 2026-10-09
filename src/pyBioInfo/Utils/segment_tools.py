# Class: SegmentTools

from collections import defaultdict
from collections.abc import Generator, Iterable, Sequence
from functools import cmp_to_key
from typing import cast
import pysam


class SegmentTools(object):
    """Static utilities operating on pysam.AlignedSegment objects."""

    @staticmethod
    def cmp(s1: pysam.AlignedSegment, s2: pysam.AlignedSegment) -> int:
        """
        Compare two aligned segments by chrom, start, end and query name.

        Both segments must be mapped. Returns -1, 0 or 1.
        """
        assert not s1.is_unmapped and not s2.is_unmapped
        chrom1, chrom2 = s1.reference_name, s2.reference_name
        assert chrom1 is not None and chrom2 is not None
        if chrom1 < chrom2:
            return -1
        elif chrom1 == chrom2:
            start1, start2 = s1.reference_start, s2.reference_start
            if start1 < start2:
                return -1
            elif start1 == start2:
                end1, end2 = s1.reference_end, s2.reference_end
                assert end1 is not None and end2 is not None
                if end1 < end2:
                    return -1
                elif end1 == end2:
                    name1, name2 = s1.query_name, s2.query_name
                    assert name1 is not None and name2 is not None
                    if name1 < name2:
                        return -1
                    elif name1 == name2:
                        return 0
        return 1

    @staticmethod
    def sort(segments: "Iterable[pysam.AlignedSegment]", by_coordinates: bool = True,
             by_name: bool = False) -> "list[pysam.AlignedSegment]":
        """
        Sort aligned segments either by coordinates (chrom, start, end,
        name) or by query name.

        :param segments: Iterable of aligned segments.
        :param by_coordinates: Sort by genomic coordinates.
        :param by_name: Sort by query name.
        :raises RuntimeError: If both or neither option is requested.
        """
        if by_coordinates:
            if by_name:
                raise RuntimeError("Cannot sort by both coordinates and name")
            else:
                return sorted(segments, key=cmp_to_key(SegmentTools.cmp))
        else:
            if by_name:
                return sorted(segments, key=lambda x: x.query_name or "")
            else:
                raise RuntimeError("Must sort by either coordinates or name")

    @staticmethod
    def get_blocks(segment: pysam.AlignedSegment, fill_deletion: bool = True) -> "list[list[int]]":
        """
        Extract the alignment blocks of a segment as [start, end] lists.

        Consecutive M (and optionally D) operations are concatenated;
        N (reference skip) breaks the blocks; insertions and clips are
        ignored.

        :param segment: An aligned segment.
        :param fill_deletion: If True, deletions (D) are treated as
                              part of the blocks.
        :raises RuntimeError: If an unknown cigar operation is found.
        """
        blocks: "list[list[int]]" = []
        cigartuples = segment.cigartuples
        assert cigartuples is not None
        y0: int = segment.reference_start
        for m, n in cigartuples:
            if m == pysam.CMATCH:  # M
                x, y = y0, y0 + n
                y0 = y
                if len(blocks) > 0 and x == blocks[-1][1]:  # concat
                    blocks[-1][1] = y
                else:
                    blocks.append([x, y])
            elif m == pysam.CDEL:  # D
                x, y = y0, y0 + n
                y0 = y
                if fill_deletion:
                    if len(blocks) > 0 and x == blocks[-1][1]:
                        blocks[-1][1] = y
                    else:
                        blocks.append([x, y])
            elif m == pysam.CREF_SKIP:  # N
                y0 = y0 + n
            elif m == pysam.CINS:  # I
                continue
            elif m == pysam.CSOFT_CLIP:  # S
                continue
            elif m == pysam.CHARD_CLIP:  # H
                continue
            else:
                raise RuntimeError("Unknown cigar type: %s" % m)
        return blocks

    @staticmethod
    def get_mapped_length(segment: pysam.AlignedSegment, fill_deletion: bool = True) -> int:
        """
        Return the length of the segment on the reference.

        M operations always count; D operations count only when
        ``fill_deletion`` is True.
        """
        length: int = 0
        cigartuples = segment.cigartuples
        assert cigartuples is not None
        for flag, count in cigartuples:
            if flag == pysam.CMATCH:  # Cigar = M
                length += count
            elif flag == pysam.CDEL:  # Cigar = D
                if fill_deletion:
                    length += count
        return length

    @staticmethod
    def get_reference_sequence(segment: pysam.AlignedSegment) -> str:
        """Return the reference sequence covered by the segment."""
        return segment.get_reference_sequence()

    @staticmethod
    def get_reference_sequence_base(segment: pysam.AlignedSegment) -> "defaultdict[str, int]":
        """Count the occurrence of each reference base (upper case)."""
        counter: "defaultdict[str, int]" = defaultdict(int)
        for base in segment.get_reference_sequence():
            counter[base.upper()] += 1
        return counter

    @staticmethod
    def get_aligned_sequence(segment: pysam.AlignedSegment) -> "tuple[str, list[int]]":
        """
        Build the read sequence aligned to the reference.

        Deletion positions are filled with "-" characters and quality
        -1; insertions and clips are removed.

        :param segment: An aligned segment.
        :return: (aligned_sequence, aligned_qualities).
        :raises RuntimeError: If an unknown cigar operation is found.
        """
        # Pad deletions with dashes and quality -1; remove insertions.
        query_sequence = segment.query_sequence
        query_qualities = cast("Sequence[int] | None", segment.query_qualities)
        aligned_query_sequence: "list[str]" = []
        aligned_query_qualities: "list[int]" = []
        if query_sequence is None:
            raise RuntimeError("Query sequence is None")
        if query_qualities is None:
            raise RuntimeError("Query qualities are None")

        offset = 0
        cigartuples = segment.cigartuples
        if cigartuples is None:
            raise RuntimeError("Cigar is None")
        for flag, count in cigartuples:
            if flag == pysam.CHARD_CLIP:
                continue
            elif flag == pysam.CSOFT_CLIP:
                offset = offset + count
            elif flag == pysam.CREF_SKIP:
                continue
            elif flag == pysam.CDEL:
                aligned_query_sequence.append("-" * count)
                aligned_query_qualities.extend([-1] * count)
            elif flag == pysam.CINS:
                offset += count
            elif flag == pysam.CMATCH:
                tmp = offset + count
                aligned_query_sequence.append(query_sequence[offset:tmp])
                aligned_query_qualities.extend(query_qualities[offset:tmp])
                offset = tmp
            else:
                raise RuntimeError(
                    "Unknown cigar type: %s [%d]" % (flag, count))
        aligned_sequence = "".join(aligned_query_sequence)
        return aligned_sequence, aligned_query_qualities


    @staticmethod
    def parse_cigar(segment: "pysam.AlignedSegment | None" = None,
                    start: "int | None" = None, cigarstring: "str | None" = None) -> "list[tuple[str, tuple[int, int], tuple[int, int], tuple[int, int]]]":
        """
        Parse a CIGAR into a list of per-operation index tuples.

        Provide exactly one of the following input forms:

        - ``segment``: parse the CIGAR of an aligned segment; or
        - ``start`` and ``cigarstring``: parse a raw CIGAR string (this
          requires the optional third-party ``cigar`` library).

        :param segment: An aligned segment carrying a ``cigartuples``.
        :param start: Reference start coordinate (required together with
                      ``cigarstring``).
        :param cigarstring: A CIGAR string to parse.
        :return: A list of ``(op, mapped, read, genomic)`` tuples, one
                 per CIGAR operation, where

                 - ``op`` is the operation character (``M``/``I``/``D``/``N``/``S``/``H``);
                 - ``mapped`` is the mapped-index interval ``(start, end)``
                   (accumulated by ``M`` and ``D``);
                 - ``read`` is the read-index interval ``(start, end)``
                   (accumulated by ``M``/``I``/``S``/``H``);
                 - ``genomic`` is the reference-index interval ``(start, end)``
                   (accumulated by ``M``/``D``/``N``).

        :raises RuntimeError: If the arguments are missing or ambiguous,
                 if the segment has no CIGAR, or if a CIGAR operation is
                 unsupported.
        """

        items: "list[tuple[str, int]]"= []
        op: str = ""
        operation: int = 0
        length: int = 0

        if segment is None:
            if cigarstring is None:
                raise RuntimeError("Must provide one of segment or cigarstring")
            else:
                if start is None:
                    raise RuntimeError("Must provide start with cigarstring")
                from cigar import Cigar  # pyright: ignore[reportMissingTypeStubs]
                for item in Cigar(cigarstring).items():
                    length, op = cast("tuple[int, str]", item)
                    items.append((op, length))
        else:
            if cigarstring is None:
                cigartuples: list[tuple[int, int]] | None = segment.cigartuples
                if cigartuples is None:
                    raise RuntimeError("Can not get cigartuples from segment")
                start = segment.reference_start
                for operation, length in cigartuples:
                    if operation == pysam.CMATCH:
                        op = "M"
                    elif operation == pysam.CINS:
                        op = "I"
                    elif operation == pysam.CDEL:
                        op = "D"
                    elif operation == pysam.CREF_SKIP:
                        op = "N"
                    elif operation == pysam.CSOFT_CLIP:
                        op = "S"
                    elif operation == pysam.CHARD_CLIP:
                        op = "H"
                    else:
                        raise RuntimeError("Unsupported cigar operation: %s" % operation)
                    items.append((op, length))
            else:
                raise RuntimeError("Can not provide both of segment and cigarstring")

        parsed_cigars: list[tuple[str, tuple[int, int], tuple[int, int], tuple[int, int]]] = []
        mapped_start: int = 0
        mapped_end: int = 0
        read_start: int = 0
        read_end: int = 0
        genomic_start: int = start
        genomic_end: int = genomic_start

        for op, length in items:
            mapped_start = mapped_end
            read_start = read_end
            genomic_start = genomic_end
            if op == "H" or op == "S" or op == "I":
                read_end = read_end + length
            elif op == "M":
                read_end = read_end + length
                mapped_end = mapped_end + length
                genomic_end = genomic_end + length
            elif op == "D":
                genomic_end = genomic_end + length
                mapped_end = mapped_end + length
            elif op == "N":
                genomic_end = genomic_end + length
            else:
                raise RuntimeError("Unsupported cigar operation: %s" % op)

            parsed_cigars.append((
                op, 
                (mapped_start, mapped_end), 
                (read_start, read_end), 
                (genomic_start, genomic_end)
            ))
        
        return parsed_cigars


    @staticmethod
    def parse_md_tag(segment: pysam.AlignedSegment) -> list[tuple[str, tuple[int, int], str | None]]:
        """
        Parse the MD tag of a segment.

        :param segment: An aligned segment carrying an MD tag.
        :return: List of (kind, mapped_indexes, bases) tuples, where
                 kind is one of "=" (match), "X" (mismatch) or "D"
                 (deletion); ``bases`` is None for matches and holds
                 the reference bases otherwise.
        :raises RuntimeError: If the MD tag contains an invalid token.
        """
        results: list[tuple[str, tuple[int, int], str | None]] = []
        mis: int = 0  # mapped_index_start
        mie: int = 0  # mapped_index_end
        md: str = cast("str", segment.get_tag("MD"))  # pyright: ignore[reportUnknownMemberType]
        i: int = 0
        while i < len(md):
            if "0" <= md[i] <= "9":
                j = i + 1
                while j < len(md) and "0" <= md[j] <= "9":
                    j += 1
                mie = mis + int(md[i:j])
                results.append(("=", (mis, mie), None))
                mis = mie
                i = j
            elif "A" <= md[i] <= "Z":
                j = i + 1
                while j < len(md) and "A" <= md[j] <= "Z":
                    j += 1
                mie = mis + j - i
                results.append(("X", (mis, mie), md[i:j]))
                mis = mie
                i = j
            elif md[i] == "^":
                j = i + 1
                while j < len(md) and "A" <= md[j] <= "Z":
                    j += 1
                mie = mis + j - i - 1
                results.append(("D", (mis, mie), md[i + 1:j]))
                mis = mie
                i = j
            else:
                raise RuntimeError("Invalid MD tag token at %d: %s" % (i, md[i]))
        return results

    @staticmethod
    def get_clipped(segment: pysam.AlignedSegment, soft_clip: bool = True,
                    hard_clip: bool = True, stranded: bool = False) -> "tuple[int, int]":
        """
        Return the number of clipped bases at both ends of the segment.

        :param segment: An aligned segment.
        :param soft_clip: Whether soft clips count.
        :param hard_clip: Whether hard clips count.
        :param stranded: If True, swap the two values on the reverse
                         strand so that they are reported in read
                         orientation (5' clip first).
        """
        c1: int = 0
        c2: int = 0
        cigars: list[tuple[int, int]] | None = segment.cigartuples
        assert cigars is not None
        if cigars[0][0] == pysam.CHARD_CLIP:
            if hard_clip:
                c1 = cigars[0][1]
        elif cigars[0][0] == pysam.CSOFT_CLIP:
            if soft_clip:
                c1 = cigars[0][1]
        if cigars[-1][0] == pysam.CHARD_CLIP:
            if hard_clip:
                c2 = cigars[-1][1]
        elif cigars[-1][0] == pysam.CSOFT_CLIP:
            if soft_clip:
                c2 = cigars[-1][1]
        if stranded and segment.is_reverse:
            c1, c2 = c2, c1
        return c1, c2

    @classmethod
    def get_query_base(cls, segment: pysam.AlignedSegment, position: int,
                       parsed_cigar: "list[tuple[str, tuple[int, int], tuple[int, int], tuple[int, int]]] | None" = None) -> "str | None":
        """
        Return the read base aligned at a reference position.

        :param segment: An aligned segment.
        :param position: Reference position (0-based).
        :param parsed_cigar: Pre-parsed cigar (parsed on demand if None).
        :return: The read base ("-" for deletions, None for reference
                 skips or uncovered positions).
        :raises RuntimeError: If the cigar contains an unexpected
                              operation at the position.
        """
        if parsed_cigar is None:
            parsed_cigar = cls.parse_cigar(segment)
        base: "str | None" = None
        for k, _, read_idxs, contig_idxs in parsed_cigar:
            if contig_idxs[1] <= position:
                continue
            elif contig_idxs[0] > position:
                break
            else:
                if k == "M":
                    read_idx = position - contig_idxs[0] + read_idxs[0]
                    if parsed_cigar[0][0] == "H":
                        i = read_idx - parsed_cigar[0][2][1]
                    else:
                        i = read_idx
                    query_sequence = segment.query_sequence
                    assert query_sequence is not None
                    base = query_sequence[i]
                elif k == "D":
                    base = "-"
                elif k == "N":
                    base = None
                else:
                    raise RuntimeError("Unexpected cigar operation: %s" % k)
                break

        return base

    @staticmethod
    def load_same_name_segments(path_or_loader: "str | Iterable[pysam.AlignedSegment]") -> "Generator[list[pysam.AlignedSegment], None, None]":
        """
        Group an aligned-segment stream by query name.

        :param path_or_loader: Either a path to a BAM/SAM file or an
                               iterable of segments sorted by query
                               name.
        :yield: Lists of segments sharing the same query name.

        The input must be sorted by query name; a decrease raises an
        AssertionError.
        """
        segments: "list[pysam.AlignedSegment]" = []

        f: "pysam.AlignmentFile | None" = None
        loader: "Iterable[pysam.AlignedSegment]"
        if isinstance(path_or_loader, str):
            f = pysam.AlignmentFile(path_or_loader)
            loader = f.fetch(until_eof=True)
        else:
            loader = path_or_loader

        for s in loader:
            if len(segments) == 0:
                segments = [s]
            else:
                assert s.query_name is not None and segments[0].query_name is not None
                if s.query_name == segments[0].query_name:
                    segments.append(s)
                else:
                    assert s.query_name > segments[0].query_name
                    yield segments
                    segments = [s]

        if len(segments) > 0:
            yield segments

        if f is not None:
            f.close()

    @staticmethod
    def is_mate(segment1: pysam.AlignedSegment, segment2: pysam.AlignedSegment) -> bool:
        """
        Test whether two segments are mates of the same pair.

        Checks query name, mate positions, read1/read2 flags,
        secondary flags and the MC (mate cigar) tags.

        :raises RuntimeError: If the read1/read2 or MC tag information
                              is inconsistent.
        """
        # read name
        if segment1.query_name != segment2.query_name:
            return False

        # chrom
        if segment1.next_reference_name != segment2.reference_name:
            return False
        if segment1.reference_name != segment2.next_reference_name:
            return False

        # start
        if segment1.next_reference_start != segment2.reference_start:
            return False
        if segment1.reference_start != segment2.next_reference_start:
            return False

        # read 1 and read 2
        if segment1.is_read1:
            if segment2.is_read1:
                return False
            elif segment2.is_read2:
                pass
            else:
                raise RuntimeError(
                    "segment1 is read1 but segment2 is neither read1 nor read2")
        elif segment1.is_read2:
            if segment2.is_read1:
                pass
            elif segment2.is_read2:
                return False
            else:
                raise RuntimeError(
                    "segment1 is read2 but segment2 is neither read1 nor read2")
        else:
            raise RuntimeError(
                "segment1 is neither read1 nor read2")

        # primary and secondary
        if segment1.is_secondary:
            if segment2.is_secondary:
                pass
            else:
                return False
        else:
            if segment2.is_secondary:
                return False
            else:
                pass

        # mate cigar
        if segment1.has_tag("MC"):
            if segment2.has_tag("MC"):
                mc1 = cast("str", segment1.get_tag("MC"))  # pyright: ignore[reportUnknownMemberType]
                mc2 = cast("str", segment2.get_tag("MC"))  # pyright: ignore[reportUnknownMemberType]
                if mc1 != segment2.cigarstring:
                    return False
                if segment1.cigarstring != mc2:
                    return False
            else:
                raise RuntimeError(
                    "segment1 has an MC tag but segment2 does not")
        else:
            if segment2.has_tag("MC"):
                raise RuntimeError(
                    "segment2 has an MC tag but segment1 does not")
            else:
                pass

        return True

