import math

from meta_length import SeqIO


def LimitSequence(min_limit, max_limit = 1000000000000):
    for lim in InnerLimitSequence(min_limit, max_limit):
        if lim >= max_limit:
            yield max_limit
            break
        else:
            yield lim

def InnerLimitSequence(min_limit, max_limit = 1000000000000):
    while min_limit < max_limit:
        yield min_limit
        yield min_limit * 3 / 2
        yield min_limit * 2
        yield 3 * min_limit
        yield 5 * min_limit
        yield 7 * min_limit
        min_limit *= 10
    yield max_limit

def Count(values, num, insert_size):
    sz = len(values)
    if sz == 0:
        return 0
    return sum([float(a[1] - insert_size) / a[0] for a in values]) * num / sz

def CountD(values, num, insert_size):
    sz = len(values)
    if sz <= 1:
        return 0
    E = Count(values, num, insert_size)
    E2 = sum([float(a[1] - insert_size) * (a[1] - insert_size) / a[0] / a[0] for a in values]) * num * num / sz
    return abs(E2 - E * E)

class LongReadRecord:
    def __init__(self, l, bins = 10):
        self.bin_size = (l + bins - 1) / bins
        self.cov = [0] * bins
        self.size = l

    def __len__(self):
        return self.size

    def update(self, pos):
        self.cov[pos / self.bin_size] += 1

    def get(self):
        if sum(self.cov) < 2:
            return sum(self.cov)
        tmp = sorted(self.cov)
        step = len(tmp) / 10
        return float(sum(tmp[step:-step])) / (len(tmp) - 2 * step) * len(tmp)

def ISCounter(max_value = 1000):
    return ShiftCounter(lambda rec: abs(rec.tlen) if abs(rec.tlen) < max_value else -1)

# def LenCounter(max_value = 1000):
#     return ShiftCounter(lambda rec: len(rec.seq), max)

class ShiftCounter:
    def __init__(self, rec_to_shift):
        self.shift_sum = 0
        self.shift_count = 0
        self.rec_to_shift = rec_to_shift

    def Process(self, rec):
        if not rec.proper_alignment:
            return
        shift = self.rec_to_shift(rec)
        if shift != -1:
            self.shift_count += 1
            self.shift_sum += shift

    def get(self):
        if self.shift_count == 0:
            return 0.
        else:
            return float(self.shift_sum) / self.shift_count

class EstimationResult:
    def __init__(self, est, disp, nonzero, insert_size):
        self.est = est
        self.disp = disp
        self.nonzero = nonzero
        self.insert_size = insert_size

class Calculator:
    def __init__(self, tslr_file, min_len, is_cnt, log):
        self.min_len = min_len
        log.info("Reading TSLRs")
        file_type = self.GetFileType(tslr_file)
        self.coverage_records = [LongReadRecord(len(rec)) for rec in SeqIO.parse(open(tslr_file, "r"), file_type)]
        self.is_cnt = is_cnt
        self.tslr_count = len([a for a in self.coverage_records if len(a) > self.min_len])

    def GetFileType(self, tslr_file):
        file_type = "fastq"
        if "fasta" in tslr_file.split(".") or "fa" in tslr_file.split("."):
            file_type = "fasta"
        return file_type

    def Process(self, rec):
        if not rec.is_unmapped:
            self.coverage_records[rec.tid].update(rec.pos)

    def Count(self, read_num, tslr_num):
        subcov = [(a.get(), a.size) for a in self.coverage_records if a.size > self.min_len][0:tslr_num]
        subsubcov = [a for a in subcov if a[0] > 0]
        insert_size = self.is_cnt.get()
        estimation = Count(subsubcov, read_num, insert_size)
        dispersion = math.sqrt(CountD(subsubcov, read_num, insert_size) / tslr_num)
        nonzero_fraction = float(len(subsubcov)) / tslr_num
        return EstimationResult(estimation, dispersion, nonzero_fraction, insert_size)
