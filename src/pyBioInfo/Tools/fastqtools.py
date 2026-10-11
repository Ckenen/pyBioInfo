from typing import Any, IO
from io import TextIOWrapper
import argparse
import subprocess
import sys
import re
import gzip
from typing import cast
import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_cmd(cmd: str) -> int:
    logging.info("RUN CMD: %s" % cmd)
    return subprocess.check_call(cmd, shell=True)

class FastqTools(object):
    @classmethod
    def fastq_to_txt(cls, args: argparse.Namespace) -> None:
        fq1 = cast(str | None, args.input1)
        fq2 = cast(str | None, args.input2)
        out = cast(str | None, args.output)

        if out is None:
            fw = sys.stdout
        elif out.endswith(".gz"):
            fw = gzip.open(out, "wt")
        else:
            fw = open(out, "w+")

        p1, p2 = None, None
        f1, f2 = None, None
        if fq1 is None:
            f1 = sys.stdin
        elif fq1.endswith(".gz"):
            p1 = subprocess.Popen("gzip -d -c %s" % fq1, shell=True, stdout=subprocess.PIPE, text=True)
            f1 = p1.stdout
            if f1 is None:
                raise RuntimeError("Failed to open process!")
        else:
            f1 = open(fq1, "r")

        if fq2 is None:
            # single-end
            pass
        elif fq2.endswith(".gz"):
            p2 = subprocess.Popen("gzip -d -c %s" % fq2, shell=True, stdout=subprocess.PIPE, text=True)
            f2 = p2.stdout
            if f2 is None:
                raise RuntimeError("Failed to open process!")
        else:
            f2 = open(fq2, "r")
        
        if f2 is None:
            name, seq, qua = "", "", ""
            for i, line in enumerate(f1):
                j = i % 4
                if j == 0:
                    name = line[:-1].split()[0]
                elif j == 1:
                    seq = line[:-1]
                elif j == 2:
                    pass
                else:
                    qua = line[:-1]
                    _ = fw.write("%s\t%s\t%s\n" % (seq, name, qua))
        else:
            name1, name2, seq1, seq2, qua1, qua2 = "", "", "", "", "", ""
            for i, (line1, line2) in enumerate(zip(f1, f2)):
                j = i % 4
                if j == 0:
                    name1 = line1[:-1].split()[0]
                    name2 = line2[:-1].split()[0]
                    if name1 != name2:
                        if name1.endswith("/1") and name2.endswith("/2"):
                            name1 = name1[:-2]
                            name2 = name2[:-2]
                        elif name1.endswith("/2") and name2.endswith("/1"):
                            name1 = name1[:-2]
                            name2 = name2[:-2]
                        elif name1.endswith(".1") and name2.endswith(".2"):
                            name1 = name1[:-2]
                            name2 = name2[:-2]
                        elif name1.endswith(".2") and name2.endswith(".1"):
                            name1 = name1[:-2]
                            name2 = name2[:-2]
                        else:
                            raise RuntimeError("Invalid name!")
                    if name1 != name2:
                        raise RuntimeError("Invalid name!")
                elif j == 1:
                    seq1, seq2 = line1[:-1], line2[:-1]
                elif j == 2:
                    pass
                else:
                    qua1, qua2 = line1[:-1], line2[:-1]
                    _ = fw.write("%s\t%s\t%s\t%s\t%s\n" % (seq1, seq2, name1, qua1, qua2))

        if f1:
            _ = f1.close()
        if f2:
            _ = f2.close()
        if fw:
            _ = fw.close()

    @classmethod
    def get_and_remove_read_size(cls, name: str, report_read_size: bool=True) -> tuple[str, int]:
        if report_read_size:
            ret = re.search(r"\[RS=(\d+)\]", name)
            if ret is not None:
                s = ret.group()
                return name.replace(s, ""), int(s[4:-1])
        return name, 1

    @classmethod
    def write_fq_line_to_pipe(cls, fields: list[str], fw1: TextIOWrapper | IO[Any], fw2: TextIOWrapper | IO[Any] | None=None, read_size: int=1, report_read_size: bool=True) -> None:
        if len(fields) == 3:
            seq1, name, qua1 = fields
            if report_read_size:
                name = "%s[RS=%d]" % (name, read_size)
            _ = fw1.write("%s\n%s\n%s\n%s\n" % (name, seq1, "+", qua1))
        elif len(fields) == 5:
            seq1, seq2, name, qua1, qua2 = fields
            if report_read_size:
                name = "%s[RS=%d]" % (name, read_size)
            _ = fw1.write("%s\n%s\n%s\n%s\n" % (name, seq1, "+", qua1))
            if fw2 is None:
                raise RuntimeError("fw2 is None!")
            _ = fw2.write("%s\n%s\n%s\n%s\n" % (name, seq2, "+", qua2))
        else:
            raise RuntimeError("Invalid fields!")


    @classmethod
    def txt_to_fastq(cls, args: argparse.Namespace) -> None:
        infile = cast(str | None, args.input)
        outfile1 = cast(str | None, args.output1)
        outfile2 = cast(str | None, args.output2)
        threads = 4
        report_read_size = True
        f = None
        fw1 = None
        fw2 = None
        p1 = None
        p2 = None

        if infile is None:
            f = sys.stdin
        elif infile.endswith(".gz"):
            p = subprocess.Popen("gzip -d -c %s" % infile, shell=True, stdout=subprocess.PIPE, text=True)
            f = p.stdout
            if f is None:
                raise RuntimeError("Failed to open process!")
        else:
            f = open(infile, "r")

        if outfile1 is None:
            fw1 = sys.stdout
        elif outfile1.endswith(".gz"):
            p1 = subprocess.Popen("pigz -p %d -c > %s" % (threads, outfile1), shell=True, stdin=subprocess.PIPE, text=True)
            fw1 = p1.stdin
            if fw1 is None:
                raise RuntimeError("Failed to open process!")
        else:
            fw1 = open(outfile1, "w+")
        
        if outfile2 is None:
            fw2 = None
        elif outfile2.endswith(".gz"):
            p2 = subprocess.Popen("pigz -p %d -c > %s" % (threads, outfile2), shell=True, stdin=subprocess.PIPE, text=True)
            fw2 = p2.stdin
            if fw2 is None:
                raise RuntimeError("Failed to open process!")
        else:
            fw2 = open(outfile2, "w+")

        last = None
        read_size = 0
        reads_in = 0
        reads_out = 0

        if fw2 is None:
            for i, line in enumerate(f):
                if i % 100000 == 0:
                    logging.info(f"Processing line {i}")
                reads_in += 1
                row = line.strip("\n").split("\t")
                new_name, rs = cls.get_and_remove_read_size(row[1], report_read_size)
                row[1] = new_name
                if last is None or row[0] != last[0]:
                    if last is not None:
                        cls.write_fq_line_to_pipe(fields=last, fw1=fw1, read_size=read_size, report_read_size=report_read_size)
                        reads_out += 1
                    last = row
                    read_size = rs
                else:
                    read_size += rs

            if last is not None:
                cls.write_fq_line_to_pipe(fields=last, fw1=fw1, read_size=read_size, report_read_size=report_read_size)
                reads_out += 1
        else:
            for i, line in enumerate(f):
                if i % 100000 == 0:
                    logging.info(f"Processing line {i}")
                reads_in += 1
                row = line.strip().split("\t")
                new_name, rs = cls.get_and_remove_read_size(row[2], report_read_size)
                row[2] = new_name
                if last is None or row[0] != last[0] or row[1] != last[1]:
                    if last is not None:
                        cls.write_fq_line_to_pipe(fields=last, fw1=fw1, fw2=fw2, read_size=read_size, report_read_size=report_read_size)
                        reads_out += 1
                    last = row
                    read_size = rs
                else:
                    read_size += rs
            if last is not None:
                cls.write_fq_line_to_pipe(fields=last, fw1=fw1, fw2=fw2, read_size=read_size, report_read_size=report_read_size)
                reads_out += 1

        if fw1:
            _ = fw1.close()
        if fw2:
            _ = fw2.close()
        if f:
            _ = f.close()       

        logging.info(f"Input reads: {reads_in}") 
        logging.info(f"Output reads: {reads_out} ({reads_out / reads_in * 100:.2f}%)")


    @classmethod
    def dedup_fastq(cls, args: argparse.Namespace) -> None:
        fq1 = cast(str, args.input1)
        fq2 = cast(str | None, args.input2)
        out1 = cast(str, args.output1)
        out2 = cast(str | None, args.output2)
        threads = cast(int, args.threads)
        buffer_size = "8G"
        report_read_size=True

        logging.info("Starting fastq dedup")

        logging.info(f"Input fastq 1: {fq1}")
        logging.info(f"Input fastq 2: {fq2}")
        logging.info(f"Output fastq 1: {out1}")
        logging.info(f"Output fastq 2: {out2}")
        logging.info(f"Threads: {threads}")
        logging.info(f"Buffer size: {buffer_size}")
        logging.info(f"Report read size: {report_read_size}")


        tmpdir = out1 + ".fastqtools_dedup_fastq_tmp"
        _ = run_cmd(f"mkdir -p {tmpdir}")

        if fq2 is None:
            cmd1 = f"fastqtools.py fastq-to-txt -1 {fq1}"
            # cmd2 = f"sort -k1,1 -k2,2 -T {tmpdir} -S {buffer_size} --parallel={threads}"
            cmd2 = f"sort -T {tmpdir} -S {buffer_size} --parallel={threads}"
            cmd3 = f"fastqtools.py txt-to-fastq -t {threads} -o {out1}"
        else:
            cmd1 = f"fastqtools.py fastq-to-txt -1 {fq1} -2 {fq2}"
            # cmd2 = f"sort -k1,1 -k2,2 -k3,3 -T {tmpdir} -S {buffer_size} --parallel={threads}"
            cmd2 = f"sort -T {tmpdir} -S {buffer_size} --parallel={threads}"
            cmd3 = f"fastqtools.py txt-to-fastq -t {threads} -o {out1} -p {out2}"
        cmd = "%s | %s | %s" % (cmd1, cmd2, cmd3)
        _ = run_cmd(cmd)
        _ = run_cmd(f"rm -r {tmpdir}")

        logging.info("Fastq dedup finished")
