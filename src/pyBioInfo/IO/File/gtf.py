# Class: GtfFile (reader of GTF format files)

import gzip
import os
import sys
from collections import OrderedDict, defaultdict
from collections.abc import Iterable, Iterator, Sequence
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
from ...Utils import BlockTools

"""
GTF format:
colum-number    content    values/format
1	chromosome name	chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,M} or GRC accession a
2	annotation source	{ENSEMBL,HAVANA}
3	feature type	{gene,transcript,exon,CDS,UTR,start_codon,stop_codon,Selenocysteine}
4	genomic start location	integer-value (1-based)
5	genomic end location	integer-value
6	score(not used)	.
7	genomic strand	{+,-}
8	genomic phase (for CDS features)	{0,1,2,.}
9	additional information as key-value pairs	see below (GENCODE)

reference: https://www.gencodegenes.org/pages/data_format.html
"""

# Record


class GtfRecord(GRange):
    """
    GtfRecord is a class for a single record (line) of GTF files.

    Besides the interval coordinates it keeps the annotation columns
    (source, feature, score, strand, frame) and the parsed attributes
    dictionary.
    """

    def __init__(self, chrom: str, start: int, end: int, strand: str,
                 score: str, frame: str, source: str, feature: str,
                 attributes: "OrderedDict[str, str]") -> None:
        """
        Initialize a GTF record.

        :param chrom: Chromosome name.
        :param start: Start coordinate (0-base, inclusive).
        :param end: End coordinate (0-base, exclusive).
        :param strand: Strand: "+", "-" or ".".
        :param score: Score column value (usually ".").
        :param frame: Frame column value ("0", "1", "2" or ".").
        :param source: Source column value, e.g. "ENSEMBL".
        :param feature: Feature type, e.g. "gene", "exon" or "CDS".
        :param attributes: Parsed attributes dictionary.
        """
        super(GtfRecord, self).__init__(chrom=chrom,
                                        start=start,
                                        end=end,
                                        strand=strand)
        self.score: "int | float | str | None" = score
        self.frame: "int | str | None" = frame
        self.source: "str | None" = source
        self.feature: "str | None" = feature
        self.attributes: "OrderedDict[str, str]" = attributes

    @override
    def format(self, fmt: str = "GTF") -> str:  # pyright: ignore[reportIncompatibleMethodOverride]
        """
        Format the record as a line of the given text format.

        :param fmt: Format name, case-insensitive; only "GTF" is supported.
        :return: Tab-separated GTF line.
        :raises NotImplementedError: If the format is unsupported.
        """
        fmt = fmt.upper()
        if fmt == "GTF":
            attri_str = "; ".join(["%s \"%s\"" % (k, v)
                                   for k, v in self.attributes.items()]) + ";"
            values: "list[object]" = [self.chrom, self.source, self.feature,
                                      self.start + 1, self.end, self.score,
                                      self.strand, self.frame, attri_str]
            return "\t".join(map(str, values))
        else:
            raise NotImplementedError("Format %s is not supported!" % fmt)


class GtfFile(BaseFile):
    """
    GtfFile is a class for reading and writing GTF files.

    Plain text GTF, gzip-compressed GTF (``.gtf.gz``) and block-
    compressed GTF with a tabix index (``.gtf.gz`` + ``.tbi``) are
    supported. For indexed files random access via :py:meth:`fetch`
    is available; otherwise records can only be read sequentially
    (iteration). Comment lines encountered while reading sequentially
    are collected into :py:attr:`comments`. Each line is parsed into
    a :py:class:`GtfRecord` by :py:meth:`parse_gtf_string`.
    """

    # Attributes assigned by the untyped BaseFile base class; declared
    # here to keep their types known.
    _path: str
    _mode: str
    _handle: "TextIO | pysam.TabixFile | None"

    def __init__(self, path: str, mode: str = "r",
                 random: "bool | None" = None) -> None:
        """
        Initialize and open a GTF file object.

        :param path: Path to the GTF (or compressed .gtf.gz) file.
        :param mode: File mode; "r" for reading and "w" for writing.
        :param random: Enable random access. If None (default), it is
                       enabled automatically for indexed .gtf.gz files.
        :raises AssertionError: If the random-access preconditions
                                are not met.
        """
        # Enable random access automatically for indexed files.
        if mode == "r" and random is None:
            random = path.endswith(".gz") and os.path.exists(path + ".tbi")

        if random:
            assert mode == "r", "Random access is only supported in read mode!"
            assert path.endswith(".gz"), \
                "Random access requires a block gzip (.gz) file!"
            assert os.path.exists(path), "File %s does not exist!" % path
            assert os.path.exists(path + ".tbi"), \
                "The tabix index (.tbi) file is missing!"
        super(GtfFile, self).__init__(path, mode)
        # Comment lines ("#...") collected while reading sequentially.
        self.comments: "list[str]" = []
        self._random: bool = bool(random)
        self.open()

    @property
    def random(self) -> bool:
        """Return True if random access is enabled for this file."""
        return self._random

    @override
    def open(self) -> None:
        """
        Open the underlying file handle.

        Random-access files are opened with pysam (tabix); plain and
        gzip-compressed files are opened as text streams. Calling open
        on an already-opened file is a no-op.

        :raises ValueError: If the mode is unsupported.
        """
        if self._handle is None:
            if self._random:
                self._handle = pysam.TabixFile(self._path)
            elif self._mode == "r":
                if self._path.endswith(".gz"):
                    self._handle = cast("TextIO", gzip.open(self._path, "rt"))
                else:
                    self._handle = open(self._path)
            elif self._mode == "w":
                if self._path.endswith(".gz"):
                    self._handle = cast("TextIO", gzip.open(self._path, "wt"))
                else:
                    self._handle = open(self._path, "w+")
            else:
                raise ValueError("Mode %s is not supported!" % self._mode)

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

    @override
    def fetch(self, chrom: "str | None" = None, start: "int | None" = None,
              end: "int | None" = None) -> "Iterator[GtfRecord]":
        """
        Fetch records of the given chromosome (region).

        For random-access files records of the given region are fetched
        via tabix, and chrom=None iterates over all contigs; for
        sequential files only chrom=None (the whole file) is allowed,
        and comment lines are collected into :py:attr:`comments`.

        :param chrom: Chromosome name; None iterates over all contigs.
        :param start: Region start (inclusive); None for the whole chromosome.
        :param end: Region end (exclusive); None for the whole chromosome.
        :raises RuntimeError: If random access is not enabled or the
                              file is not opened.
        :return: Iterator of GtfRecord objects.
        """
        if self._random:
            handle = self._handle
            assert handle is not None, "File is not opened!"
            assert isinstance(handle, pysam.TabixFile)
            if chrom is None:
                for contig in sorted(handle.contigs):
                    for line in cast("Iterator[str]", handle.fetch(contig)):
                        yield self.parse_gtf_string(line)
            else:
                for line in cast("Iterator[str]", handle.fetch(chrom, start, end)):
                    yield self.parse_gtf_string(line)
        else:
            if chrom is None:
                handle = self._handle
                assert handle is not None, "File is not opened!"
                for line in cast("TextIO", handle):
                    line = line.strip("\n")
                    if line == "":
                        continue
                    if line.startswith("#"):
                        self.comments.append(line)
                    else:
                        yield self.parse_gtf_string(line)
            else:
                raise RuntimeError(
                    "Random access is not supported! Please set random=True.")

    @classmethod
    def parse_gtf_attribute_string(cls, s: str) -> "OrderedDict[str, str]":
        """
        Parse the attribute column of a GTF line.

        Attributes are ``key "value";`` pairs separated by semicolons
        (GENCODE style).

        :param s: The attribute column string.
        :return: Ordered dictionary of attribute key-value pairs.
        """
        attributes: "OrderedDict[str, str]" = OrderedDict()
        key_x, key_y = None, None
        value_x, value_y = None, None
        value_has_quotes = False
        i = 0
        while i < len(s):
            c = s[i]
            if key_x is None:
                if c != " " and c != ";":
                    key_x = i
            else:
                if key_y is None:
                    if c == " ":
                        key_y = i
                else:
                    if value_x is None:
                        if c != " ":
                            if c == "\"":
                                value_x = i + 1
                                value_has_quotes = True
                            else:
                                value_x = i
                                value_has_quotes = False
                    else:
                        if value_has_quotes:
                            if c == "\"" and s[i - 1] != "\\" and s[i + 1] == ";":
                                value_y = i
                                attributes[s[key_x:key_y]] = s[value_x:value_y]
                                key_x, key_y = None, None
                                value_x, value_y = None, None
                        else:
                            if c == ";":
                                value_y = i
                                attributes[s[key_x:key_y]] = s[value_x:value_y]
                                key_x, key_y = None, None
                                value_x, value_y = None, None
            i += 1
        return attributes

    @classmethod
    def parse_gtf_string(cls, line: str) -> GtfRecord:
        """
        Parse a single GTF line into a GtfRecord.

        :param line: One line of a GTF file.
        :return: The parsed GtfRecord instance.
        :raises RuntimeError: If the line is malformed.
        """
        try:
            values = line.split("\t")
            assert len(values) == 9
            chrom = values[0]
            source = values[1]
            feature = values[2]
            start = int(values[3]) - 1
            end = int(values[4])
            score = values[5]
            strand = values[6]
            if strand == "?":
                strand = "."
            frame = values[7]
            attributes_str = values[8]
            assert strand == "+" or strand == "-" or strand == "."
            assert frame in ["0", "1", "2", "."]
            attributes = cls.parse_gtf_attribute_string(attributes_str)
            assert "gene_id" in attributes
            if feature != "gene":
                assert "transcript_id" in attributes
            record = GtfRecord(chrom=chrom,
                               start=start,
                               end=end,
                               strand=strand,
                               score=score,
                               frame=frame,
                               source=source,
                               feature=feature,
                               attributes=attributes)
            return record
        except Exception as e:
            raise RuntimeError("Error in parsing line: %s" % line) from e


# Transcript


class GtfTranscript(TRange):
    """
    GtfTranscript is a class for transcripts of GTF files.

    A transcript is a thick range whose blocks are the exons and whose
    thick interval spans the CDS features. The original GTF records are
    kept in :py:attr:`records`, grouped by feature type.
    """

    def __init__(self, chrom: str, start: int, end: int, name: "str | None",
                 strand: "str | None", blocks: "Iterable[Iterable[int]] | None",
                 thick: "tuple[int, int] | None",
                 records: "defaultdict[str, list[GtfRecord]]") -> None:
        """
        Initialize a GTF transcript.

        :param chrom: Chromosome name.
        :param start: Overall start (inclusive).
        :param end: Overall end (exclusive).
        :param name: Transcript ID.
        :param strand: Strand: None, ".", "+" or "-".
        :param blocks: Exon intervals.
        :param thick: CDS interval, or None for non-coding transcripts.
        :param records: Original records grouped by feature type.
        """
        super(GtfTranscript, self).__init__(chrom=chrom,
                                            start=start,
                                            end=end,
                                            name=name,
                                            strand=strand,
                                            blocks=blocks,
                                            thick=thick)
        # {"feature": [record1, record2], ...}
        self.records: "defaultdict[str, list[GtfRecord]]" = records


class GtfTranscriptBuilder(object):
    """
    GtfTranscriptBuilder is a class for building transcripts from GTF
    records.

    The structure of feature type:
    gene -> transcript -> exon and CDS
    """

    def __init__(self, records: "Iterable[GtfRecord]") -> None:
        """
        Initialize the builder.

        :param records: GTF records sorted by chromosome and coordinates.
        """
        self.records: "Iterable[GtfRecord]" = records

    @classmethod
    def construct_gtf_transcript(cls, records: "Sequence[GtfRecord]") -> GtfTranscript:
        """
        Construct a transcript from its records (exons and CDS).

        :param records: Records sharing the same transcript ID.
        :return: The constructed GtfTranscript.
        :raises AssertionError: If no exon record exists.
        """
        assert len(records) > 0
        transcript_id = records[0].attributes["transcript_id"]
        transcript_records: "defaultdict[str, list[GtfRecord]]" = defaultdict(list)
        for record in records:
            assert record.attributes["transcript_id"] == transcript_id
            feature = record.feature
            assert feature is not None
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
        transcript = GtfTranscript(chrom=chrom,
                                   start=start,
                                   end=end,
                                   name=transcript_id,
                                   strand=strand,
                                   blocks=blocks,
                                   thick=thick,
                                   records=transcript_records)
        return transcript

    def process_one_chromosome(self, records: "Sequence[GtfRecord]") -> "list[GtfTranscript]":
        """
        Build all transcripts of one chromosome.

        :param records: Records of a single chromosome.
        :return: Sorted list of GtfTranscript objects.
        """
        # {"gene_id": record, "gene_id": record, ...}
        genes: "dict[str, GtfRecord]" = dict()
        # {"transcript_id": [record 1, record 2, record 3, ...]}
        container: "defaultdict[str, list[GtfRecord]]" = defaultdict(list)
        for record in records:
            if record.feature == "gene":
                genes[record.attributes["gene_id"]] = record
            else:
                container[record.attributes["transcript_id"]].append(record)

        transcripts: "list[GtfTranscript]" = []
        for transcript_id in container.keys():
            transcript = self.construct_gtf_transcript(
                container[transcript_id])
            gene_id = container[transcript_id][0].attributes["gene_id"]
            if gene_id in genes.keys():
                assert "gene" not in transcript.records.keys()
                transcript.records["gene"] = [genes[gene_id]]
            transcripts.append(transcript)
        transcripts.sort()
        return transcripts

    def __iter__(self) -> "Iterator[GtfTranscript]":
        """
        Iterate over transcripts, chromosome by chromosome.

        :return: Iterator of GtfTranscript objects.
        """
        chrom: "str | None" = None
        array: "list[GtfRecord]" = []
        for record in self.records:
            if chrom is None:
                chrom = record.chrom
                array = [record]
            else:
                if record.chrom == chrom:
                    array.append(record)
                else:
                    for transcript in self.process_one_chromosome(array):
                        yield transcript
                    chrom = record.chrom
                    array = [record]
        if chrom is not None:
            for transcript in self.process_one_chromosome(array):
                yield transcript


# Gene


class GtfGene(TRange):
    """
    GtfGene is a class for genes of GTF files.

    A gene groups the transcripts sharing the same gene ID. Its
    coordinates span all transcripts; the original records are kept in
    :py:attr:`records` grouped by feature type.
    """

    def __init__(self, chrom: str, start: int, end: int, name: "str | None",
                 strand: "str | None", blocks: "Iterable[Iterable[int]] | None",
                 thick: "tuple[int, int] | None",
                 transcripts: "list[GtfTranscript]",
                 records: "defaultdict[str, list[GtfRecord]]") -> None:
        """
        Initialize a GTF gene.

        :param chrom: Chromosome name.
        :param start: Gene start (inclusive).
        :param end: Gene end (exclusive).
        :param name: Gene ID.
        :param strand: Strand: None, ".", "+" or "-".
        :param blocks: Exon intervals of all transcripts, sutured.
        :param thick: Not used for genes (None).
        :param transcripts: Transcripts belonging to this gene.
        :param records: Original records grouped by feature type.
        """
        super(GtfGene, self).__init__(chrom=chrom,
                                      start=start,
                                      end=end,
                                      strand=strand,
                                      name=name,
                                      blocks=blocks,
                                      thick=thick)
        # [Transcript 1, Transcript 2, Transcript 3, ...]
        self.transcripts: "list[GtfTranscript]" = transcripts
        # {
        #   "gene": [record],
        #   "transcript": [record 1, record 2, ...],
        #   "exon": [record 1, record 2, ...],
        #   ...
        # }
        self.records: "defaultdict[str, list[GtfRecord]]" = records


class GtfGeneBuilder(object):
    """
    GtfGeneBuilder is a class for building genes from GTF records.

    The structure of feature type:
    gene -> transcript -> exon and CDS
    """

    def __init__(self, records: "Iterable[GtfRecord]") -> None:
        """
        Initialize the builder.

        :param records: GTF records sorted by chromosome and coordinates.
        """
        self.records: "Iterable[GtfRecord]" = records

    def process_one_chromosome(self, records: "Sequence[GtfRecord]") -> "list[GtfGene]":
        """
        Build all genes of one chromosome.

        :param records: Records of a single chromosome.
        :return: Sorted list of GtfGene objects.
        """
        # {"gene_id": [record1, record2, record3]}
        container: "defaultdict[str, list[GtfRecord]]" = defaultdict(list)
        for record in records:
            container[record.attributes["gene_id"]].append(record)

        genes: "list[GtfGene]" = []
        for gene_id in container.keys():
            array = container[gene_id]
            # {
            #   "gene": [record],
            #   "transcript": [record 1, record 2, ...],
            #   "exon": [record 1, record 2, ...],
            #   ...
            # }
            gene_records: "defaultdict[str, list[GtfRecord]]" = defaultdict(list)
            for record in array:
                feature = record.feature
                assert feature is not None
                gene_records[feature].append(record)

            # {"transcript_id": [record1, record2, record3, ...]}
            transcript_records: "defaultdict[str, list[GtfRecord]]" = defaultdict(list)
            transcripts: "list[GtfTranscript]" = []
            for record in array:
                if record.feature == "gene":
                    continue
                transcript_records[record.attributes["transcript_id"]].append(record)

            for transcript_id in transcript_records:
                transcript = GtfTranscriptBuilder.construct_gtf_transcript(
                    transcript_records[transcript_id])
                if "gene" in gene_records:
                    transcript.records["gene"] = gene_records["gene"]
                transcripts.append(transcript)
            transcripts.sort()
            chrom = transcripts[0].chrom
            start = min([transcript.start for transcript in transcripts])
            end = max([transcript.end for transcript in transcripts])
            strand = transcripts[0].strand
            blocks: "list[tuple[int, int]]" = []
            for transcript in transcripts:
                for block in transcript.blocks:
                    blocks.append(block)
            blocks = list(sorted(blocks))
            blocks = BlockTools.suture(blocks)

            thick = None
            gene = GtfGene(chrom=chrom,
                           start=start,
                           end=end,
                           strand=strand,
                           name=gene_id,
                           records=gene_records,
                           transcripts=transcripts,
                           blocks=blocks,
                           thick=thick)
            genes.append(gene)
        genes.sort()
        return genes

    def __iter__(self) -> "Iterator[GtfGene]":
        """
        Iterate over genes, chromosome by chromosome.

        :return: Iterator of GtfGene objects.
        """
        chrom: "str | None" = None
        array: "list[GtfRecord]" = []
        for record in self.records:
            if chrom is None:
                chrom = record.chrom
                array = [record]
            else:
                if record.chrom == chrom:
                    array.append(record)
                else:
                    for gene in self.process_one_chromosome(array):
                        yield gene
                    chrom = record.chrom
                    array = [record]
        if chrom is not None:
            for gene in self.process_one_chromosome(array):
                yield gene
