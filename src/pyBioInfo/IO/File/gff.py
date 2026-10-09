import gzip
import os
import sys
from collections import OrderedDict, defaultdict
from collections.abc import Generator, Iterable, Iterator, Sequence
from typing import TextIO, cast

if sys.version_info >= (3, 12):
    from typing import override
else:
    from typing import Callable  # pyright: ignore[reportUnreachable]

    def override(func: Callable) -> Callable:
        """No-op fallback for Python < 3.12 compatibility."""
        return func

import pysam

from .base import BaseFile
from ...Range import GRange, TRange

"""
GFF format:
1   chrom
2   source
3   feature
4   start (1-base)
5   end (1-base, not included)
6   score
7   strand
8   frame
9   attributes  ID=AAA;Parent=BBB
"""

# Record


class GffRecord(GRange):
    """
    GffRecord is a class for a single record (line) of GFF files.

    Besides the interval coordinates it keeps the annotation columns
    (source, feature, score, strand, frame) and the parsed attributes
    dictionary. The record name is taken from the ``ID`` attribute.
    """

    # Optional child records attached by GffTreeBuilder when building
    # the feature hierarchy (gene -> transcript -> exon/CDS).
    childs: "list[GffRecord] | None" = None

    def __init__(self, chrom: str, start: int, end: int, name: "str | None",
                 strand: "str | None", score: "int | str | None",
                 frame: "int | str | None", source: "str | None",
                 feature: "str | None",
                 attributes: "OrderedDict[str, str]") -> None:
        """
        Initialize a GFF record.

        :param chrom: Chromosome name.
        :param start: Start coordinate (0-base, inclusive).
        :param end: End coordinate (0-base, exclusive).
        :param name: Record name (the ID attribute).
        :param strand: Strand: None, ".", "+" or "-".
        :param score: Score column value.
        :param frame: Frame column value.
        :param source: Source column value.
        :param feature: Feature type, e.g. "gene", "exon" or "CDS".
        :param attributes: Parsed attributes dictionary.
        """
        super(GffRecord, self).__init__(chrom=chrom, start=start, end=end, name=name, strand=strand)
        self.score: "int | float | str | None" = score
        self.frame: "int | str | None" = frame
        self.source: "str | None" = source
        self.feature: "str | None" = feature
        self.attributes: "OrderedDict[str, str]" = attributes

    @override
    def __str__(self) -> str:
        """Return the GFF formatted line of this record."""
        return self.format(fmt="GFF")

    @override
    def format(self, fmt: str = "GFF") -> str:  # pyright: ignore[reportIncompatibleMethodOverride]
        """
        Format the record as a line of the given text format.

        :param fmt: Format name, case-insensitive; only "GFF" is supported.
        :return: Tab-separated GFF line.
        :raises NotImplementedError: If the format is unsupported.
        """
        fmt = fmt.upper()
        if fmt == "GFF":
            attri_str = ";".join(["%s=%s" % (k, v)
                                  for k, v in self.attributes.items()])
            values: "list[object]" = [self.chrom, self.source, self.feature,
                                      self.start + 1, self.end, self.score,
                                      self.strand, self.frame, attri_str]
            return "\t".join(map(str, values))
        else:
            raise NotImplementedError("Format %s is not supported!" % fmt)


class GffFile(BaseFile):
    """
    GffFile is a class for reading GFF files.

    Both plain text GFF and block-compressed GFF (``.gff.gz`` with a
    ``.tbi`` index) are supported. For indexed files random access via
    :py:meth:`fetch` is available; otherwise records can only be read
    sequentially (iteration).

    Each line is parsed into a :py:class:`GffRecord` by
    :py:meth:`parse_gff_string`.
    """

    # Attributes assigned by the untyped BaseFile base class; declared
    # here to keep their types known.
    _path: str
    _mode: str
    _handle: "TextIO | pysam.TabixFile | None"

    def __init__(self, path: str, mode: str = "r", random: "bool | None" = None) -> None:
        """
        Initialize a GFF file object.

        :param path: Path to the GFF (or block-compressed .gff.gz) file.
        :param mode: File mode; only reading ("r") is supported.
        :param random: Enable random access. If None (default), it is
                       enabled automatically for indexed .gff.gz files.
        :raises ValueError: If the mode is unsupported.
        """
        super(GffFile, self).__init__(path, mode)
        if mode != "r":
            raise ValueError("Mode %s is not supported!" % mode)
        if random is None:
            random = path.endswith(".gff.gz") and os.path.exists(path + ".tbi")
        self._random: bool = bool(random)

    @property
    def random(self) -> bool:
        """Return True if random access is enabled for this file."""
        return self._random

    @override
    def open(self) -> None:
        """
        Open the underlying file handle.

        Plain GFF files are opened as text streams; block-compressed
        files are opened with pysam (tabix) when random access is
        requested, otherwise as a gzip text stream. Calling open on an
        already-opened file is a no-op.
        """
        if self._handle is None:
            if self._path.endswith(".gz"):
                if self._random:
                    self._handle = pysam.TabixFile(self._path)
                else:
                    self._handle = cast("TextIO", gzip.open(self._path, "rt"))
            else:
                self._handle = open(self._path)

    @override
    def close(self) -> None:
        """
        Close the underlying file handle.

        Calling close on an already-closed file is a no-op, so that
        explicit close followed by ``__del__`` does not raise.
        """
        if self._handle is not None:
            self._handle.close()
            # Note: the narrowed handle type is intentionally incompatible
            # with the heterogeneous "IO[Any] | None" of the base class.
            self._handle = None  # pyright: ignore[reportIncompatibleVariableOverride]

    @classmethod
    def parse_gff_string(cls, string: str) -> GffRecord:
        """
        Parse a single GFF line into a GffRecord.

        :param string: One line of a GFF file.
        :return: The parsed GffRecord instance.
        :raises AssertionError: If the line is malformed or has no ID
                                attribute.
        """
        items = string.rstrip("\n").split("\t")
        assert len(items) == 9
        chrom, base, feature, start, end, score, strand, frame, attributes = items
        assert strand in ("+", "-", ".")
        attributes = OrderedDict((key, value.strip()) for key, value in
                                 (item.split("=", 1) for item in attributes.split(";")))
        assert "ID" in attributes
        return GffRecord(chrom=chrom,
                         start=int(start) - 1,
                         end=int(end),
                         name=attributes["ID"],
                         strand=strand,
                         score=score,
                         frame=frame,
                         feature=feature,
                         source=base,
                         attributes=attributes)

    @override
    def __iter__(self) -> "Generator[GffRecord, None, None]":
        """
        Iterate over all records of the file sequentially.

        :raises RuntimeError: If the file is not opened.
        :return: Iterator of GffRecord objects.
        """
        if self._handle is not None:
            for line in cast("TextIO", self._handle):
                if line.startswith("#") or line == "\n":
                    continue
                yield self.parse_gff_string(line)
        else:
            raise RuntimeError("File is not opened!")

    @override
    def fetch(self, chrom: "str | None" = None, start: "int | None" = None, end: "int | None" = None) -> "Iterator[GffRecord]":
        """
        Fetch records of the given chromosome (region).

        :param chrom: Chromosome name; None iterates over all contigs.
        :param start: Region start (inclusive); None for the whole chromosome.
        :param end: Region end (exclusive); None for the whole chromosome.
        :raises RuntimeError: If random access is not enabled or the
                              file is not opened.
        :return: Iterator of GffRecord objects.
        """
        if self._random:
            handle = self._handle
            assert handle is not None
            assert isinstance(handle, pysam.TabixFile)
            if chrom is None:
                for contig in sorted(handle.contigs):
                    for line in cast("Iterator[str]", handle.fetch(contig)):
                        yield self.parse_gff_string(line)
            else:
                for line in cast("Iterator[str]", handle.fetch(chrom, start, end)):
                    yield self.parse_gff_string(line)
        else:
            raise RuntimeError("Random access is not supported! Please use a block gzip file and create with random=True.")


# Transcript


class GffTranscript(TRange):
    """
    GffTranscript is a class for transcript of GFF files.

    A transcript is a thick range whose blocks are the exons and whose
    thick interval spans the CDS features. The original GFF records are
    kept in :py:attr:`records`, grouped by feature type.
    """

    def __init__(self, chrom: str, start: int, end: int, name: "str | None",
                 strand: "str | None", blocks: "Iterable[Iterable[int]] | None",
                 thick: "tuple[int, int] | None",
                 records: "dict[str, list[GffRecord]]") -> None:
        """
        Initialize a GFF transcript.

        :param chrom: Chromosome name.
        :param start: Overall start (inclusive).
        :param end: Overall end (exclusive).
        :param name: Transcript ID.
        :param strand: Strand: None, ".", "+" or "-".
        :param blocks: Exon intervals.
        :param thick: CDS interval, or None for non-coding transcripts.
        :param records: Original records grouped by feature type.
        """
        super(GffTranscript, self).__init__(chrom=chrom,
                                            start=start,
                                            end=end,
                                            name=name,
                                            strand=strand,
                                            blocks=blocks,
                                            thick=thick)
        # {"feature": [record1, record2], ...}
        self.records: "dict[str, list[GffRecord]]" = records

    @override
    def format(self, fmt: str = "GFF", ncol: int = 12,
               attris: "list[str] | None" = None,
               header: "Iterable[str] | None" = None) -> str:
        """
        Format the transcript as lines of the given text format.

        For "GFF" all records (gene, transcript, exons and CDS) are
        formatted and joined by newlines; other formats are delegated
        to the base class.

        :param fmt: Format name, case-insensitive.
        :param ncol: Number of columns for the BED format.
        :param attris: Extra attribute keys for the GFF format.
        :param header: Column keys for the TBL format.
        :return: Formatted text line(s).
        """
        fmt = fmt.upper()
        if fmt == "GFF":
            lines: "list[str]" = []
            for key in self.records.keys():
                for record in self.records[key]:
                    lines.append(record.format("GFF"))
            return "\n".join(lines)
        else:
            return super(GffTranscript, self).format(fmt=fmt, ncol=ncol,
                                                     attris=attris, header=header)


class GffTranscriptBuilder(object):
    """
    The structure of feature type:
    gene -> transcript -> exon and CDS
    """

    def __init__(self, records: "Iterable[GffRecord]") -> None:
        """
        Initialize the builder.

        :param records: GFF records sorted by chromosome and coordinates.
        """
        self.records: "Iterable[GffRecord]" = records

    @classmethod
    def construct_gff_transcript(cls, records: "Sequence[GffRecord]") -> GffTranscript:
        """
        Construct a transcript from its child records (exons and CDS).

        :param records: Records sharing the same Parent transcript ID.
        :return: The constructed GffTranscript.
        :raises AssertionError: If no exon record exists.
        """
        assert len(records) > 0
        transcript_id = records[0].attributes["Parent"]
        transcript_records: "dict[str, list[GffRecord]]" = dict()
        for record in records:
            assert record.attributes["Parent"] == transcript_id
            feature = record.feature
            assert feature is not None
            if feature not in transcript_records:
                transcript_records[feature] = list()
            transcript_records[feature].append(record)
        assert "exon" in transcript_records
        exons = transcript_records["exon"]
        chrom = exons[0].chrom
        strand = exons[0].strand
        blocks = [(exon.start, exon.end) for exon in exons]
        blocks.sort()
        start, end = blocks[0][0], blocks[-1][1]
        thick: "tuple[int, int] | None" = None
        if "CDS" in transcript_records:
            cdss = transcript_records["CDS"]
            blocks_cds = [(cds.start, cds.end) for cds in cdss]
            blocks_cds.sort()
            thick = (blocks_cds[0][0], blocks_cds[-1][1])
        transcript = GffTranscript(chrom=chrom,
                                   start=start,
                                   end=end,
                                   name=transcript_id,
                                   strand=strand,
                                   blocks=blocks,
                                   thick=thick,
                                   records=transcript_records)
        return transcript

    def process_one_chromosome(self, records: "Sequence[GffRecord]") -> "list[GffTranscript]":
        """
        Build all transcripts of one chromosome.

        :param records: Records of a single chromosome.
        :return: Sorted list of GffTranscript objects.
        :raises AssertionError: If transcript IDs are duplicated.
        """
        # {ID: record, ...}
        id_to_record_dict: "dict[str, GffRecord]" = dict()
        # {Parent ID: [Child 1, Child 2, Child 3, ...], ...}
        parent_to_childs_dict: "defaultdict[str, list[GffRecord]]" = defaultdict(list)
        # {Transcript ID 1, Transcript ID 2, Transcript ID 3, ...}
        transcript_id_set: "set[str]" = set()
        for record in records:
            id_to_record_dict[record.attributes["ID"]] = record
            if "Parent" in record.attributes:
                parent_to_childs_dict[record.attributes["Parent"]].append(
                    record)
            if record.feature == "transcript":
                # Make sure the transcript id is unique.
                assert record.attributes["ID"] not in transcript_id_set
                transcript_id_set.add(record.attributes["ID"])

        transcripts: "list[GffTranscript]" = []
        for transcript_id in transcript_id_set:
            transcript = self.construct_gff_transcript(
                parent_to_childs_dict[transcript_id])
            transcript_record = id_to_record_dict[transcript_id]
            assert transcript_record.feature == "transcript" or transcript_record.feature == "mRNA"
            feature = transcript_record.feature
            assert feature is not None
            transcript.records[feature] = [transcript_record]
            if "Parent" in transcript_record.attributes:
                parent_id = transcript_record.attributes["Parent"]
                if parent_id in id_to_record_dict:
                    parent = id_to_record_dict[parent_id]
                    assert parent.feature == "gene"
                    parent_feature = parent.feature
                    assert parent_feature is not None
                    transcript.records[parent_feature] = [parent]
            transcripts.append(transcript)
        transcripts.sort()
        return transcripts

    def __iter__(self) -> "Iterator[GffTranscript]":
        """
        Iterate over transcripts, chromosome by chromosome.

        :raises RuntimeError: If the records are not sorted by chromosome.
        :return: Iterator of GffTranscript objects.
        """
        buffer: "list[GffRecord] | None" = None  # Records of current chromosome
        for record in self.records:
            if buffer is not None:
                if record.chrom == buffer[0].chrom:
                    buffer.append(record)
                    continue
                if record.chrom < buffer[0].chrom:
                    raise RuntimeError("Records is not sorted by chrom!")
                for transcript in self.process_one_chromosome(buffer):
                    yield transcript
            buffer = [record]
        if buffer is not None:
            for transcript in self.process_one_chromosome(buffer):
                yield transcript

# Gene


class GffGene(TRange):
    """
    GffGene is a class for gene of GFF files.

    A gene groups the transcripts sharing the same gene ID. Its
    coordinates span all transcripts; the original records are kept in
    :py:attr:`records` grouped by feature type.
    """

    def __init__(self, chrom: str, start: int, end: int, name: "str | None",
                 strand: "str | None", blocks: "Iterable[Iterable[int]] | None",
                 thick: "tuple[int, int] | None",
                 transcripts: "list[GffTranscript]",
                 records: "defaultdict[str, list[GffRecord]]") -> None:
        """
        Initialize a GFF gene.

        :param chrom: Chromosome name.
        :param start: Gene start (inclusive).
        :param end: Gene end (exclusive).
        :param name: Gene ID.
        :param strand: Strand: None, ".", "+" or "-".
        :param blocks: Not used for genes (None).
        :param thick: Not used for genes (None).
        :param transcripts: Transcripts belonging to this gene.
        :param records: Original records grouped by feature type.
        """
        super(GffGene, self).__init__(chrom=chrom,
                                      start=start,
                                      end=end,
                                      strand=strand,
                                      name=name,
                                      blocks=blocks,
                                      thick=thick)
        # [Transcript 1, Transcript 2, Transcript 3, ...]
        self.transcripts: "list[GffTranscript]" = transcripts
        # {
        #   "gene": [record],
        #   "transcript": [record1, record2, ...],
        #   "exon": [record1, record2, ...],
        #   "CDS": [record1, record2, ...],
        #   ....
        # }
        self.records: "defaultdict[str, list[GffRecord]]" = records


class GffGeneBuilder(object):
    """
    GffGeneBuilder is a class for building genes from GFF records.

    The structure of feature type:
    gene -> transcript -> exon and CDS
    """

    def __init__(self, records: "Iterable[GffRecord]") -> None:
        """
        Initialize the builder.

        :param records: GFF records sorted by chromosome and coordinates.
        """
        self.records: "Iterable[GffRecord]" = records

    def process_one_chromosome(self, records: "Sequence[GffRecord]") -> "list[GffGene]":
        """
        Build all genes of one chromosome.

        :param records: Records of a single chromosome.
        :return: Sorted list of GffGene objects.
        :raises AssertionError: If the gene coordinates mismatch the
                                 transcripts.
        """
        id_to_record_dict: "dict[str, GffRecord]" = dict()
        parent_child_dict: "defaultdict[str, list[GffRecord]]" = defaultdict(list)
        gene_id_set: "set[str]" = set()
        for record in records:
            id_to_record_dict[record.attributes["ID"]] = record
            if record.feature == "gene":
                gene_id_set.add(record.attributes["ID"])
            if "Parent" in record.attributes:
                parent_child_dict[record.attributes["Parent"]].append(record)

        genes: "list[GffGene]" = []
        for gene_id in gene_id_set:
            gene_record = id_to_record_dict[gene_id]
            assert gene_record.feature == "gene"
            gene_records: "defaultdict[str, list[GffRecord]]" = defaultdict(list)
            gene_records[gene_record.feature] = [gene_record]

            transcripts: "list[GffTranscript]" = []
            for transcript_record in parent_child_dict[gene_id]:
                assert transcript_record.feature == "transcript" or transcript_record.feature == "mRNA"
                feature = transcript_record.feature
                assert feature is not None
                gene_records[feature].append(transcript_record)
                transcript_id = transcript_record.attributes["ID"]
                transcript_record_list = parent_child_dict[transcript_id]
                for record in transcript_record_list:
                    assert record.feature is not None
                    gene_records[record.feature].append(record)
                transcript = GffTranscriptBuilder.construct_gff_transcript(transcript_record_list)
                transcript_feature = transcript_record.feature
                assert transcript_feature is not None
                transcript.records[transcript_feature] = [
                    transcript_record]
                gene_feature = gene_record.feature
                assert gene_feature is not None
                transcript.records[gene_feature] = [gene_record]
                transcripts.append(transcript)
            transcripts.sort()

            chrom = transcripts[0].chrom
            start = min([transcript.start for transcript in transcripts])
            assert start == gene_record.start
            end = max([transcript.end for transcript in transcripts])
            assert end == gene_record.end
            strand = transcripts[0].strand
            blocks: "Iterable[Iterable[int]] | None" = None
            thick: "tuple[int, int] | None" = None
            gene = GffGene(chrom=chrom,
                           start=start,
                           end=end,
                           strand=strand,
                           name=gene_id,
                           records=gene_records,
                           transcripts=transcripts,
                           blocks=blocks, thick=thick)
            genes.append(gene)
        genes.sort()
        return genes

    def __iter__(self) -> "Iterator[GffGene]":
        """
        Iterate over genes, chromosome by chromosome.

        :raises RuntimeError: If the records are not sorted by chromosome.
        :return: Iterator of GffGene objects.
        """
        buffer: "list[GffRecord] | None" = None  # Records of current chromosome
        for record in self.records:
            if buffer is not None:
                if record.chrom == buffer[0].chrom:
                    buffer.append(record)
                    continue
                if record.chrom < buffer[0].chrom:
                    raise RuntimeError("Records is not sorted by chrom!")
                for gene in self.process_one_chromosome(buffer):
                    yield gene
            buffer = [record]
        if buffer is not None:
            for gene in self.process_one_chromosome(buffer):
                yield gene

# Tree


class GffTreeBuilder(object):
    """
    GffTreeBuilder is a class for building feature trees from GFF records.

    Records are linked through the ``Parent`` attribute; top level
    records (without Parent) are returned, each carrying its child
    records in the ``childs`` attribute.
    """

    def __init__(self, records: "Iterable[GffRecord]") -> None:
        """
        Initialize the builder.

        :param records: GFF records sorted by chromosome and coordinates.
        """
        self.records: "Iterable[GffRecord]" = records

    def process_one_chromosome(self, records: "Sequence[GffRecord]") -> "list[GffRecord]":
        """
        Build the feature tree of one chromosome.

        :param records: Records of a single chromosome.
        :return: Sorted list of top level (no Parent) records.
        """
        id_to_record_dict: "dict[str, GffRecord]" = dict()
        parent_child_dict: "defaultdict[str, list[GffRecord]]" = defaultdict(list)
        no_parent_records: "list[GffRecord]" = []
        for record in records:
            assert "ID" in record.attributes.keys()
            id_to_record_dict[record.attributes["ID"]] = record
            if "Parent" in record.attributes:
                parent_child_dict[record.attributes["Parent"]].append(record)
            else:
                no_parent_records.append(record)
        for pid, items in parent_child_dict.items():
            id_to_record_dict[pid].childs = items

        array: "list[GffRecord]" = []
        for record in no_parent_records:
            array.append(record)
        array.sort()
        return array

    def __iter__(self) -> "Iterator[GffRecord]":
        """
        Iterate over top level records, chromosome by chromosome.

        :raises RuntimeError: If the records are not sorted by chromosome.
        :return: Iterator of top level GffRecord objects.
        """
        buffer: "list[GffRecord] | None" = None  # Records of current chromosome
        for record in self.records:
            if buffer is not None:
                if record.chrom == buffer[0].chrom:
                    buffer.append(record)
                    continue
                if record.chrom < buffer[0].chrom:
                    raise RuntimeError("Records is not sorted by chrom!")
                for record1 in self.process_one_chromosome(buffer):
                    yield record1
            buffer = [record]
        if buffer is not None:
            for record1 in self.process_one_chromosome(buffer):
                yield record1
