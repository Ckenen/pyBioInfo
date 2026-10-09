# Class: FastaFile (reader and writer of FASTA format files)

import gzip
import os
import sys
from collections.abc import Callable, Iterable, Iterator
from typing import TextIO, cast

from .base import BaseFile
from ...Range import GRange

if sys.version_info >= (3, 12):
    from typing import override
else:
    from typing import Callable  # pyright: ignore[reportUnreachable]

    def override(func: Callable) -> Callable:
        """No-op fallback for Python < 3.12 compatibility."""
        return func


REVERSE_COMPLEMENT_MAPPER = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
    "N": "N",
    "M": "K",
    "K": "M",
    "R": "Y",
    "Y": "R",
    "W": "W",
    "S": "S",
    "V": "B",
    "B": "V",
    "H": "D",
    "D": "H",
}

for _k, _v in list(REVERSE_COMPLEMENT_MAPPER.items()):
    REVERSE_COMPLEMENT_MAPPER[_k.lower()] = _v.lower()


class FastaRecord(object):
    """
    FastaRecord is a class for a single record of FASTA files.

    It holds a sequence ``name`` (without the ">" prefix) and the
    nucleotide ``sequence`` string.
    """

    def __init__(self, name: str, sequence: str) -> None:
        """
        Initialize a FASTA record.

        :param name: Name of the sequence.
        :param sequence: Nucleotide sequence string.
        """
        self.name: str = name
        self.sequence: str = sequence

    @override
    def __str__(self) -> str:
        """Format the record as a two-line FASTA text."""
        return ">%s\n%s\n" % (self.name, self.sequence)

    @property
    def seq(self) -> str:
        """Alias of :py:attr:`sequence`."""
        return self.sequence

    def reverse_complement(self) -> "FastaRecord":
        """
        Return the reverse complement of the record.

        :return: A new FastaRecord with the reversed and complemented
                 sequence.
        """
        s = "".join([REVERSE_COMPLEMENT_MAPPER[b] for b in self.sequence[::-1]])
        return FastaRecord(self.name, s)

    # Deprecated misspelled alias kept for backward compatibility.
    reverse_comlement: "Callable[..., FastaRecord]" = reverse_complement


class FastaFile(BaseFile):
    """
    FastaFile is a class for reading and writing FASTA files.

    Plain text FASTA, gzip-compressed FASTA (``.fasta.gz``) and
    indexed FASTA (``.fai``) are supported. For indexed files random
    access via :py:meth:`fetch` (generator) or :py:meth:`fast_fetch`
    (plain string) is available; otherwise records can only be read
    sequentially (iteration). Sequences are fetched from the
    reference strand, or reverse-complemented for negative-strand
    requests when ``strandness`` is enabled.
    """

    # Attributes assigned by the untyped BaseFile base class; declared
    # here to keep their types known.
    _path: str
    _mode: str
    _handle: "TextIO | None"

    def __init__(self, path: str, mode: str = "r",
                 random: "bool | None" = None) -> None:
        """
        Initialize and open a FASTA file object.

        :param path: Path to the FASTA (or compressed .fasta.gz) file.
        :param mode: File mode; "r" for reading and "w" for writing.
        :param random: Enable random access. If None (default), it is
                       enabled automatically when a ``.fai`` index
                       exists.
        :raises ValueError: If the mode is unsupported.
        :raises AssertionError: If the random-access preconditions
                                are not met.
        """
        # Enable random access automatically for indexed files.
        if mode == "r" and random is None:
            random = os.path.exists(path + ".fai")

        if random:
            assert mode == "r", "Random access is only supported in read mode!"
            assert os.path.exists(path), "File %s does not exist!" % path
            assert os.path.exists(path + ".fai"), \
                "The FASTA index (.fai) file is missing!"

        super(FastaFile, self).__init__(path, mode)
        self._random: bool = bool(random)
        # Sequence index: name -> (length, offset, line_bases, line_width).
        self._index: "dict[str, list[int]] | None" = None
        self.open()

    @property
    def random(self) -> bool:
        """Return True if random access is enabled for this file."""
        return self._random

    def _load_index(self) -> None:
        """
        Load the .fai index file into memory.

        The index maps each sequence name to its (length, byte offset,
        line bases, line width) tuple.
        """
        index: "dict[str, list[int]]" = dict()
        with open(self._path + ".fai") as f:
            for line in f:
                row = line.strip("\n").split("\t")
                seqname = row[0]
                length = int(row[1])
                offset = int(row[2])
                line_bases = int(row[3])
                line_width = int(row[4])
                index[seqname] = [length, offset, line_bases, line_width]
        self._index = index

    @override
    def open(self) -> None:
        """
        Open the underlying file handle.

        Random-access files also load the ``.fai`` index. Calling open
        on an already-opened file is a no-op.

        :raises ValueError: If the mode is unsupported.
        """
        if self._handle is None:
            if self._random:
                self._handle = open(self._path)
                self._load_index()
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

    def get_chrom_length(self, chrom: str) -> int:
        """
        Return the length of a sequence.

        :param chrom: Name of the sequence.
        :return: Length of the sequence.
        """
        index = self._index
        assert index is not None, "Index is not loaded!"
        return index[chrom][0]

    def _get_sequence_by_coordinates(self, chrom: str,
                                      start: "int | None" = None,
                                      end: "int | None" = None) -> str:
        """
        Fetch a sequence region using the .fai index.

        :param chrom: Name of the sequence.
        :param start: Region start (0-based, inclusive); None for 0.
        :param end: Region end (0-based, exclusive); None for the end
                    of the sequence.
        :return: The sequence of the region (reference strand).
        """
        assert self._random
        index = self._index
        assert index is not None, "Index is not loaded!"
        length, offset, width1, width2 = index[chrom]
        if start is None:
            start = 0
        if end is None:
            end = length
        assert start >= 0 and end <= length
        n1 = int(start / width1) * width2 + start % width1
        n2 = int((end - 1) / width1) * width2 + (end - 1) % width1 + 1
        handle = self._handle
        assert handle is not None, "File is not opened!"
        _ = handle.seek(offset + n1, 0)
        return handle.read(n2 - n1).replace("\n", "")

    def _get_sequence_by_blocks(self, chrom: str,
                                blocks: "Iterable[Iterable[int]]") -> str:
        """
        Fetch and concatenate the sequences of several regions.

        :param chrom: Name of the sequence.
        :param blocks: Iterable of (start, end) block coordinates.
        :return: The concatenated sequence (reference strand).
        """
        assert self._random
        ss: "list[str]" = []
        for start, end in blocks:
            s = self._get_sequence_by_coordinates(chrom, start, end)
            ss.append(s)
        return "".join(ss)

    @override
    def fetch(self, chrom: "str | None" = None, start: "int | None" = None,
              end: "int | None" = None, strand: str = "+",
              obj: "object | None" = None,
              strandness: bool = True) -> "Iterator[FastaRecord]":
        """
        Fetch sequences of the given chromosome (region).

        For random-access files the ``obj`` argument (an object with
        ``chrom``, ``blocks`` and ``strand`` attributes, e.g. a
        :py:class:`GRange`) takes precedence over the chrom/start/end
        coordinates, and chrom=None iterates over all sequences; for
        sequential files only chrom=None (the whole file) is allowed.
        Negative-strand requests are reverse-complemented when
        ``strandness`` is enabled.

        :param chrom: Sequence name; None iterates over all sequences.
        :param start: Region start (inclusive); None for 0.
        :param end: Region end (exclusive); None for the sequence end.
        :param strand: Strand of the request ("+" or "-").
        :param obj: Range object with chrom/blocks/strand attributes.
        :param strandness: Reverse-complement negative-strand requests.
        :raises RuntimeError: If random access is not enabled for a
                              region request.
        :return: Iterator of FastaRecord objects.
        """
        if self._random:
            if obj is not None:
                # Duck-typed range object (GRange, Fragment, ...) with
                # chrom/blocks/strand attributes.
                obj2 = cast("GRange", obj)
                strand2 = obj2.strand
                s = self._get_sequence_by_blocks(obj2.chrom, obj2.blocks)
                if strand2 == "-" and strandness:
                    s = self.reverse_complement(s)
                yield FastaRecord("unamed", s)
            else:
                index = self._index
                assert index is not None, "Index is not loaded!"
                if chrom is None:
                    for chrom2 in sorted(index):
                        s = self._get_sequence_by_coordinates(chrom2)
                        if strand == "-" and strandness:
                            s = self.reverse_complement(s)
                        yield FastaRecord(chrom2, s)
                else:
                    if start is None:
                        start = 0
                    if end is None:
                        end = self.get_chrom_length(chrom)
                    s = self._get_sequence_by_coordinates(chrom, start, end)
                    if strand == "-" and strandness:
                        s = self.reverse_complement(s)
                    yield FastaRecord("%s:%d-%s" % (chrom, start, end), s)
        else:
            if chrom is None:
                handle = self._handle
                assert handle is not None, "File is not opened!"
                name: "str | None" = None
                rows: "list[str]" = []
                for line in handle:
                    line = line.strip("\n").strip()
                    if line.startswith(">"):
                        if name is not None:
                            seq = "".join(rows)
                            yield FastaRecord(name, seq)
                        name = line[1:]
                        rows = []
                    else:
                        if line != "":
                            rows.append(line)
                if name is not None:
                    seq = "".join(rows)
                    yield FastaRecord(name, seq)
            else:
                raise RuntimeError(
                    "Random access is not supported! Please set random=True.")

    def fast_fetch(self, chrom: "str | None" = None, start: "int | None" = None,
                   end: "int | None" = None, strand: str = "+",
                   obj: "object | None" = None,
                   strandness: bool = True) -> str:
        """
        Fetch a sequence region as a plain string.

        A non-generator counterpart of :py:meth:`fetch` returning the
        sequence directly. The ``obj`` argument takes precedence over
        the chrom/start/end coordinates; chrom is mandatory for
        coordinate-based requests.

        :param chrom: Sequence name.
        :param start: Region start (inclusive); None for 0.
        :param end: Region end (exclusive); None for the sequence end.
        :param strand: Strand of the request ("+" or "-").
        :param obj: Range object with chrom/blocks/strand attributes.
        :param strandness: Reverse-complement negative-strand requests.
        :raises RuntimeError: If random access is not enabled or chrom
                              is not provided.
        :return: The sequence string.
        """
        if self._random:
            if obj is not None:
                # Duck-typed range object (GRange, Fragment, ...) with
                # chrom/blocks/strand attributes.
                obj2 = cast("GRange", obj)
                strand2 = obj2.strand
                s = self._get_sequence_by_blocks(obj2.chrom, obj2.blocks)
                if strand2 == "-" and strandness:
                    s = self.reverse_complement(s)
                return s
            else:
                if chrom is None:
                    raise RuntimeError("Must provided chrom.")
                else:
                    if start is None:
                        start = 0
                    if end is None:
                        end = self.get_chrom_length(chrom)
                    s = self._get_sequence_by_coordinates(chrom, start, end)
                    if strand == "-" and strandness:
                        s = self.reverse_complement(s)
                    return s
        else:
            raise RuntimeError(
                "Random access is not supported! Please set random=True.")

    def write(self, record: "FastaRecord | None" = None, name: "str | None" = None,
              sequence: "str | None" = None, width: "int | None" = 50) -> None:
        """
        Write a FASTA record to the file.

        Either a :py:class:`FastaRecord` or a name/sequence pair is
        given; the sequence is wrapped into lines of ``width`` bases.

        :param record: A FastaRecord to write; overrides name/sequence.
        :param name: Name of the sequence (without the ">" prefix).
        :param sequence: Nucleotide sequence string.
        :param width: Line width for wrapping; None writes the whole
                      sequence on one line.
        :raises RuntimeError: If the file is not opened.
        """
        handle = self._handle
        assert handle is not None, "File is not opened!"
        if record is not None:
            name = record.name
            sequence = record.sequence
        assert name is not None and sequence is not None
        _ = handle.write(">%s\n" % name)
        if width is None:
            _ = handle.write("%s\n" % sequence)
        else:
            for i1 in range(0, len(sequence), width):
                i2 = min(i1 + width, len(sequence))
                seq2 = sequence[i1:i2]
                _ = handle.write(seq2 + "\n")

    @classmethod
    def reverse_complement(cls, seq: str) -> str:
        """
        Reverse-complement a nucleotide sequence.

        :param seq: Nucleotide sequence string.
        :return: The reverse complement of the sequence.
        """
        return "".join([REVERSE_COMPLEMENT_MAPPER[b] for b in seq[::-1]])
