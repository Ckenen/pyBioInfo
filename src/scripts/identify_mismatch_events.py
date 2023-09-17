#!/usr/bin/env python
import sys
import os
from collections import Counter
from pyBioInfo.IO.File import BamFile, FamFile, FastaFile, SegmentTools


def get_events(align):
    events = []
    for item in align.parsed_md:
        if item[0] == "X":
            mapped_start, mapped_end = item[1]
            assert mapped_end - mapped_start == 1
            distance = min(mapped_start, align.mapped_length - mapped_start - 1)
            ref_base = item[2]
            alt_base = None
            quality = None
            position = None
            found = False
            for d in align.parsed_cigar:
                if d[0] == "M" and d[1][0] <= mapped_start and mapped_end <= d[1][1]:
                    offset = mapped_start - d[1][0]
                    read_start, read_end = d[2]
                    chrom_start, chrom_end = d[3]
                    alt_base = align.sequence[read_start + offset]
                    quality = align.qualities[read_start + offset]
                    position = chrom_start + offset
                    events.append(
                        [position, ref_base, alt_base, quality, distance])
                    found = True
                    break
            assert found
        elif item[0] == "D":
            mapped_start, mapped_end = item[1]
            for d in align.parsed_cigar:
                if d[0] == "D" and d[1][0] <= mapped_start and mapped_end <= d[1][1]:
                    assert d[1][0] == mapped_start
                    assert d[1][1] == mapped_end
                    for i in range(mapped_end - mapped_start):
                        distance = min(
                            mapped_start + i, align.mapped_length - 1 - (mapped_start + i))
                        ref_base = item[2][i]
                        position = d[3][0] + i
                        events.append([position, ref_base, "-", 0, distance])
                    break
    return events


def check_events(events, align):
    array = []
    for event in events:
        position = event[0]
        base = None
        quality = None
        distance = None
        for item in align.parsed_cigar:
            if item[3][0] <= position < item[3][1]:
                if item[0] == "M":
                    offset = position - item[3][0]
                    read_index = item[2][0] + offset
                    base = align.sequence[read_index]
                    quality = align.qualities[read_index]
                    distance = min(
                        item[1][0] + offset, align.mapped_length - 1 - (item[1][0] + offset))
                    break
                elif item[0] == "D":
                    offset = position - item[3][0]
                    base = "-"
                    quality = 0
                    distance = min(
                        item[1][0] + offset, align.mapped_length - 1 - (item[1][0] + offset))
                    break
        event.append(base)
        event.append(quality)
        event.append(distance)
        array.append(event)
    return array


def merge_events(events1, events2):
    data = dict()
    for event in events1:
        data[event[0]] = event
    for event in events2:
        if event[0] in data:
            pass
        else:
            data[event[0]] = event
    return data


def process_alignment(align, fasta):
    align.sequence = align.segment.query_sequence
    align.qualities = align.segment.query_qualities
    align.parsed_md = SegmentTools.parse_md_tag(align.segment)
    align.parsed_cigar = SegmentTools.parse_cigar(align.segment)
    align.mapped_length = align.parsed_cigar[-1][1][1]
    # ME tag
    align.events = get_events(align)
    align.segment.set_tag("ME", ";".join([",".join(map(str, e)) for e in align.events]))
    # RC tag
    s = fasta.fetch(obj=align, strandness=False).upper()
    c = Counter(s)
    s = ";".join(["%s,%d" % (x[0], x[1]) for x in c.items()])
    align.segment.set_tag("RC", s)
    return align
    
    
def process_fragment(frag, fasta):
    align1 = process_alignment(frag.mate1, fasta)
    align2 = process_alignment(frag.mate2, fasta)
    events1 = align1.events
    events2 = align2.events
    events1 = check_events(events1, align2)
    events2 = check_events(events2, align1)
    events = merge_events(events1, events2)
    # PE tag
    pe = ""
    array = []
    for position in events.keys():
        event = events[position]
        ref_base = event[1]
        alt_base = event[2]
        quality = event[3]
        distance = event[4]
        valid = True
        if event[6] is not None:
            if event[6] > quality:
                alt_base = event[5]
                quality = event[6]
                distance = event[7]
            elif event[6] == quality:
                if event[5] != alt_base:
                    valid = False
        if valid and alt_base != ref_base:
            array.append(
                ",".join(map(str, [position, ref_base, alt_base, quality, distance])))
    pe = ";".join(array)
    align1.segment.set_tag("PE", pe)
    # PC tag
    s = fasta.fetch(obj=frag, strandness=False).upper()
    c = Counter(s)
    s = ";".join(["%s,%d" % (x[0], x[1]) for x in c.items()])
    frag.mate1.segment.set_tag("PC", s)
    

def identify_mismatch_events_se(infile1, infile2, outfile):
    with BamFile(infile1) as f, FastaFile(infile2) as fasta, BamFile(outfile, "wb", f) as fw:
        for align in f:
            align = process_alignment(align, fasta)
            fw.write(align)
             

def identify_mismatch_events_pe(infile1, infile2, outfile):
    with FamFile(infile1) as f, FastaFile(infile2) as fasta, FamFile(outfile, "wb", f) as fw:
        for frag in f:
            process_fragment(frag, fasta)
            fw.write(frag)

                
def main():
    infile1 = None # input.bam or input.fam
    infile2 = None # genome.fasta
    outfile = None # output.bam or output.fam
    
    try:
        infile1, infile2, outfile = sys.argv[1:]
    except Exception:
        print("Usage: %s input.bam genome.fasta output.bam" % os.path.basename(__file__))
        exit(1)
        
    if infile1.endswith(".bam"):
        identify_mismatch_events_se(infile1, infile2, outfile)
    elif infile1.endswith(".fam"):
        identify_mismatch_events_pe(infile1, infile2, outfile)
    else:
        raise RuntimeError("Unknown input type.")
    
        
if __name__ == '__main__':
    main()
    