# Class: BaseFile (base class of file objects)

from collections.abc import Iterator
from typing import IO, Any, TypeVar

_T = TypeVar("_T", bound="BaseFile")

class BaseFile(object):
    """
    BaseFile is the base class for all file objects of pyBioInfo.

    It stores the common file state (path, mode, and the underlying
    handle) and defines the interface every subclass implements:
    ``open()`` acquires the handle, ``close()`` releases it, and
    ``fetch()`` yields the records. Iterating the object yields the
    records of ``fetch()``; the context-manager protocol
    (``with ...``) opens and closes the file automatically.

    Examples:
        >>> with BedFile("demo.bed", "r") as file:
        ...     for record in file:
        ...         ...
    """

    # Path of the file.
    _path: str
    # File mode: "r" for reading, "w" for writing, etc.
    _mode: str
    # Underlying file handle; None while the file is closed. The
    # concrete type is heterogeneous across subclasses (TextIO,
    # pysam.TabixFile, pysam.AlignmentFile, ...), hence "Any".
    _handle: "IO[Any] | None"  # pyright: ignore[reportExplicitAny]

    def __init__(self, path: str, mode: str = "r") -> None:
        """
        Initialize a file object without opening it.

        :param path: Path of the file.
        :param mode: File mode, e.g. "r" or "w".
        """
        self._path = path
        self._mode = mode
        self._handle = None

    # Access to protected attributes.

    @property
    def path(self) -> str:
        """Path of the file."""
        return self._path

    @property
    def mode(self) -> str:
        """Mode of the file, e.g. "r" or "w"."""
        return self._mode

    @property
    def handle(self) -> "IO[Any] | None":  # pyright: ignore[reportExplicitAny]
        """Underlying file handle, or None while the file is closed."""
        return self._handle

    # The interface implemented by subclasses.

    def open(self) -> None:
        """
        Open the file and acquire the underlying handle.

        :raises NotImplementedError: Always, in the base class.
        """
        raise NotImplementedError(
            "This method should be overridden in subclasses!")

    def close(self) -> None:
        """
        Close the file and release the underlying handle.

        :raises NotImplementedError: Always, in the base class.
        """
        raise NotImplementedError(
            "This method should be overridden in subclasses!")

    def fetch(self) -> "Iterator[Any]":  # pyright: ignore[reportExplicitAny]
        """
        Yield the records of the opened file.

        :raises NotImplementedError: Always, in the base class.
        """
        raise NotImplementedError(
            "This method should be overridden in subclasses!")

    # Iteration and context-manager protocols.

    def __iter__(self) -> "Iterator[Any]":  # pyright: ignore[reportExplicitAny]
        """Yield the records of ``fetch()``."""
        yield from self.fetch()

    def __enter__(self: _T) -> _T:
        """Open the file and return self (entering a with-block)."""
        self.open()
        return self

    def __exit__(self, exc_type: object, exc_val: object, exc_tb: object) -> None:
        """Close the file (leaving a with-block)."""
        self.close()

    def __del__(self) -> None:
        """Close the file when the object is garbage collected."""
        self.close()
