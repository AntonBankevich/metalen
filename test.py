import os

import metalen.metalen_io
import random

import sys

from metalen import SeqIO


def RC(s):
    res = []
    rc = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    for a in s:
        res.append(rc[a])
    return "".join(res[::-1])


def RandomGenome(len):
    res = []
    for i in range(len):
        res.append(random.choice(['A', 'C', 'G', 'T']))
    return "".join(res)

def RandomPairedReads(fname1, fname2, genome, num, l):
    f1 = open(fname1, "w")
    f2 = open(fname2, "w")
    ins = random.randint(200, 300)
    n = len(genome)
    genome = genome + genome
    for i in range(num):
        pos = random.randint(0, n - 1)
        read1 = SeqIO.SeqRecord(genome[pos:pos + l], str(i) + "_" + str(pos), "B"*l)
        read2 = SeqIO.SeqRecord(RC(genome[pos + ins - l: pos + ins]), str(i) + "_" + str(pos), "B"*l)
        SeqIO.write(read1, f1, "fastq")
        SeqIO.write(read2, f2, "fastq")
    f1.close()
    f2.close()

def RandomLongReads(fname, genome, num, len_bounds):
    f = open(fname, "w")
    N = len(genome)
    genome = genome + genome
    for i in range(num):
        pos = random.randint(0, N - 1)
        l = random.randint(len_bounds[0], len_bounds[1])
        seq = genome[pos:pos + l]
        if random.choice([True, False]):
            seq = RC(seq)
        read = SeqIO.SeqRecord(seq, str(i) + "_" + str(pos))
        SeqIO.write(read, f, "fasta")
    f.close()

if __name__ == "__main__":
    genome_length = 1000000
    short_num = 100000
    long_num = 1000
    genome = RandomGenome(genome_length)
    dir = sys.argv[1]
    metalen.metalen_io.ensure_dir_existence(dir)
    RandomLongReads(os.path.join(dir, "long.fasta"), genome, long_num, [6000, 10000])
    RandomPairedReads(os.path.join(dir, "short_1.fastq"), os.path.join(dir, "short_2.fastq"), genome, short_num, 100)
