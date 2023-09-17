from pyBioInfo.Range import GRange
from .file_stream import TextFile


class TableFile(TextFile):
    def __init__(self, path, mode, column_names):
        self.column_names = column_names
        super(TableFile, self).__init__(path, mode)

    def _read_record(self):
        handle = self.handle
        for line in handle:
            line = line.strip("\n").strip()
            if line.startswith("#"):
                self.comments.append(line)
            elif line == "":
                continue
            else:
                raise NotImplementedError()

    def _write_record(self, record):
        raise NotImplementedError()

    # def __init__(self, path, mode="r", comment="#", header=0, delimiter="\t", *args, **kwargs):
    #     super(TableFile, self).__init__(path=path, mode=mode, *args, **kwargs)
    #     if header is not None:
    #         if isinstance(header, int):
    #             assert header >= 0
    #         elif isinstance(header, str):
    #             header = header.split(delimiter)
    #             assert "chrom" in header and "start" in header
    #         elif isinstance(header, list):
    #             assert "chrom" in header and "start" in header
    #         else:
    #             raise RuntimeError()
    #
    #     self.header = header
    #     self.delimiter = delimiter
    #     self.comment = comment
    #
    # def records(self):
    #     handle = self.handle
    #     header = self.header
    #     delimiter = self.delimiter
    #
    #     for i, line in enumerate(handle):
    #         if header is None:
    #             raise RuntimeError()
    #         else:
    #             if isinstance(header, int):
    #                 if i < header:
    #                     continue
    #                 elif i == header:
    #                     header = line.strip("\t").split(delimiter)
    #                     if header[0].startswith("#"):
    #                         header[0] = header[0][1:]
    #                     self.header = header
    #                     assert "chrom" in header and "start" in header
    #                 else:
    #                     raise RuntimeError()
    #             elif isinstance(header, list):
    #                 cols = line.strip("\n").split(delimiter)
    #                 assert len(header) <= len(cols)
    #                 kwargs = dict()
    #                 for key, value in zip(self.header, cols):
    #                     if key == "thick":
    #                         if value == "":
    #                             value = None
    #                         else:
    #                             x, y = value.split("-")
    #                             value = (int(x), int(y))
    #                     elif key == "blocks":
    #                         blocks = []
    #                         for item in value.split(";"):
    #                             x, y = item.split("-")
    #                             blocks.append((int(x), int(y)))
    #                         value = blocks
    #                     else:
    #                         if "." in value:
    #                             if value == ".":
    #                                 value = None
    #                             else:
    #                                 try:
    #                                     value = float(value)
    #                                 except ValueError:
    #                                     pass
    #                         else:
    #                             try:
    #                                 value = int(value)
    #                             except ValueError:
    #                                 pass
    #                     kwargs[key] = value
    #
    #                 yield GRange(**kwargs)
    #
    # def write(self, obj, *args, **kwargs):
    #     cols = []
    #     for head in self.header:
    #         value = obj.__getattribute__(head)
    #         cols.append(value)
    #     self._handle.write(self.delimiter)
