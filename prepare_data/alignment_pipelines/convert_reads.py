#!/usr/bin/env python3

import sys
import argparse
import pysam
import random

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bed", required = True, help = "Bed file of methylation values")
    parser.add_argument("--bam", required = True, help = "Input bam or sam file")
    parser.add_argument("--out", required = True, help = "Name for output sam file")

    args = parser.parse_args()
    return args

def read_bedgraph_file(bed_file):
    my_bed_dict = {}

    file = open(bed_file, "r")
    next(file)
    for line in file:
        line = line.strip().split("\t")
        chr = line[0]
        pos1 = int(line[1])
        pos2 = int(line[2])
        meth_percent = float(line[3]) * 100

        if chr not in my_bed_dict:
            my_bed_dict[chr] = {}
        for x in range(pos1, pos2):
            my_bed_dict[chr][x] = meth_percent
    
    return my_bed_dict

def decide_to_convert(meth_percent):
    tmp = random.randrange(0,100)
    if tmp > meth_percent:
        return True
    else:
        return False

def check_alignment_positions(read, original, converted, bed_dict):

    my_seq = list(read.query_sequence)
    quals = read.query_qualities
    for x, alignment_pos in enumerate(read.get_reference_positions(full_length = True)):            
        if alignment_pos in bed_dict[read.reference_name] and decide_to_convert(bed_dict[read.reference_name][alignment_pos]) is True:
            if my_seq[x] == original:
                my_seq[x] = converted
    read.query_sequence = "".join(my_seq)
    read.query_qualities = quals

    return read


def read_bam_file(bam_file, bed_dict, out_file):
    my_bam = pysam.AlignmentFile(bam_file, "rb")
    out_bam = pysam.AlignmentFile(out_file, "wb", template=my_bam)

    for read in my_bam:
        if (read.is_reverse and read.is_read2) or (read.mate_is_reverse and read.is_read1): # 'Top' strand
            new_read = check_alignment_positions(read, "C", "T", bed_dict)
        elif (read.is_reverse and read.is_read1) or (read.mate_is_reverse and read.is_read2): # 'Bottom' strand
            new_read = check_alignment_positions(read, "G", "A", bed_dict)
        else:
            continue
        out_bam.write(new_read)



def run_script():
    args = argparser()

    bed_dict = read_bedgraph_file(args.bed)

    read_bam_file(args.bam, bed_dict, args.out)


if __name__ == "__main__":
    run_script()