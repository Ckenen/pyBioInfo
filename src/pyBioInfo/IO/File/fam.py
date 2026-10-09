# Class: FamFile (reader and writer of FAM format files)

import sys
from collections.abc import Iterator

if sys.version_info >= (3, 12):
    from typing import override
else:
    from typing import Callable  # pyright: ignore[reportUnreachable]

    def override(func: Callable) -> Callable:
        """No-op fallback for Python < 3.12 compatibility."""
        return func

import pysam

from .bam import Alignment, BamFile
from ...Range import MRange
from ...Utils import SegmentPairBuilder


class Fragment(MRange):
    """
    Fragment is a class for a pair of alignments of BAM files.

    It is a :py:class:`MRange` (multiple range) grouping the two mates
    of a paired-end fragment: the blocks of the first and second mates
    are kept as the two block arrays, and the mates themselves are
    exposed by the :py:attr:`mate1` and :py:attr:`mate2` properties.
    """

    _mate1: Alignment
    _mate2: Alignment

    def __init__(self, segment1: pysam.AlignedSegment,
                 segment2: pysam.AlignedSegment) -> None:
        """
        Build a Fragment from two paired pysam segments.

        :param segment1: The first mate segment of the fragment.
        :param segment2: The second mate segment of the fragment.
        :raises AssertionError: If the mates are not from the same
                                template on the same chromosome.
        """
        if segment1.query_name != segment2.query_name:
            raise AssertionError(f"Mates are not from the same template: {segment1.query_name} != {segment2.query_name}")
        if segment1.reference_name != segment2.reference_name:
            raise AssertionError(f"Mates are not on the same chromosome: {segment1.reference_name} != {segment2.reference_name}")
        mate1 = Alignment(segment1)
        mate2 = Alignment(segment2)
        super(Fragment, self).__init__(
            chrom=mate1.chrom, 
            name=mate1.name, 
            strand=mate1.strand,
            blocks_array=[mate1.blocks, mate2.blocks])
        self._mate1 = mate1
        self._mate2 = mate2

    @property
    def mate1(self) -> Alignment:
        """The first mate alignment of the fragment."""
        return self._mate1

    @property
    def mate2(self) -> Alignment:
        """The second mate alignment of the fragment."""
        return self._mate2


class FamFile(BamFile):
    """
    FamFile is a class for reading and writing paired-end BAM files.

    It behaves like a :py:class:`BamFile`, except that records are
    yielded as :py:class:`Fragment` objects pairing the two mates of
    each template. Random-access files pair the mates with
    SegmentPairBuilder; sequential files pair records in file order,
    two by two. Fragments are written as their two mate segments.
    """

    @override
    def fetch(self, chrom: str | None = None,  # pyright: ignore[reportIncompatibleMethodOverride]
              start: int | None = None,
              end: int | None = None) -> Iterator[Fragment]:
        """
        Fetch the fragments of the given chromosome (region).

        For random-access files fragments of the given region are
        fetched via the index, and chrom=None iterates over all
        references; for sequential files only chrom=None (the whole
        file) is allowed, and records are paired two by two in file
        order.

        :param chrom: Chromosome name; None iterates over all references.
        :param start: Region start (inclusive); None for the whole chromosome.
        :param end: Region end (exclusive); None for the whole chromosome.
        :raises RuntimeError: If random access is not enabled or the
                              file is not opened.
        :return: Iterator of Fragment objects.
        """
        handle = self._handle
        if handle is None:
            raise RuntimeError("File is not opened!")
        if self._random:
            references = self._references
            if references is None:
                raise RuntimeError("File is not opened!")   
            if chrom is None:
                for contig in sorted(references):
                    for pair in SegmentPairBuilder(handle.fetch(contig=contig), is_sorted_by_coords=True):
                        yield Fragment(pair.mate1, pair.mate2)
            else:
                for pair in SegmentPairBuilder(handle.fetch(contig=chrom, start=start, stop=end), is_sorted_by_coords=True):
                    yield Fragment(pair.mate1, pair.mate2)
        else:
            if chrom is None:
                segment_last = None
                for i, segment in enumerate(handle.fetch(until_eof=True)):  
                    if i % 2 == 1:
                        if segment_last is None:
                            raise RuntimeError("The file is not paired-end!")
                        if segment_last.query_name != segment.query_name:
                            raise RuntimeError(f"Mates are not from the same template: {segment_last.query_name} != {segment.query_name}")
                        yield Fragment(segment_last, segment)
                    segment_last = segment
            else:
                raise RuntimeError("Random access is not supported! Please set random=True.")

    @override
    # Note: Fragment is an MRange (not an Alignment), so the record
    # parameter is intentionally narrower than the base class one.
    def write(self, obj: "Fragment | Alignment | pysam.AlignedSegment") -> None:  # pyright: ignore[reportIncompatibleMethodOverride]
        """
        Write a fragment (or single alignment) to the file.

        Fragments are written as their two mate segments; other record
        types are delegated to the base class.

        :param obj: A Fragment, Alignment or raw pysam AlignedSegment.
        :raises RuntimeError: If the file is not opened.
        """
        handle = self._handle
        assert handle is not None, "File is not opened!"
        if isinstance(obj, Fragment):
            _ = handle.write(obj.mate1.segment)
            _ = handle.write(obj.mate2.segment)
        else:
            super(FamFile, self).write(obj)
