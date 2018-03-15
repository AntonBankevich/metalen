import logging
import math
import os

from metalen import SeqIO
from metalen import sam_parser
from metalen import alignment
from metalen import io
from metalen import calculator

VERSION = "1.0"


def GetReadNum(read_id):
    tmp = read_id.split()[0].split(".")
    tmp1 = read_id.split("|")
    if len(tmp) > 1 and tmp[1].isdigit():
        return 2 * int(read_id.split()[0].split(".")[1])
    elif len(tmp1) > 2 and tmp1[2][1:].isdigit():
        return 2 * int(tmp1[2][1:])
    else:
        return 0


def estToString(num):
    if num > 1e9 / 2:
        return "{0:.2f}Gb".format(num * 1e-9)
    if num > 1e6 / 2:
        return "{0:.2f}Mb".format(num * 1e-6)
    if num > 1e3 / 2:
        return "{0:.2f}Kb".format(num * 1e-3)
    return "{0:.2f}".format(num)


class ResultPrinter:
    def __init__(self, calc, log, debug):
        self.calc = calc
        self.log = log
        self.debug = debug
        self.limits, self.tslr_limits = self.generate_output_params(calc.tslr_count)
        self.cur_limit = 0
        self.debug = debug

    def generate_output_params(self, tslr_count):
        from metalen.calculator import LimitSequence
        if self.debug:
            return list(LimitSequence(1000)), LimitSequence(10, self.calc.tslr_count)
        else:
            return [], LimitSequence(tslr_count, self.calc.tslr_count)

    def Process(self, rec):
        if rec.query_name != "" and GetReadNum(rec.query_name) > self.limits[self.cur_limit]:
            self.print_all_results(self.limits[self.cur_limit])
            self.cur_limit += 1

    def print_all_results(self, read_count):
        from metalen.calculator import LimitSequence
        if self.debug:
            tslr_limits = list(LimitSequence(self.tslr_limits, self.calc.tslr_count))
        else:
            tslr_limits = [self.calc.tslr_count]
        self.log.info("\nResults for " + str(read_count))
        for num_tslrs in tslr_limits:
            self.print_results(read_count, num_tslrs, self.log)

    def print_results(self, num_reads, num_tslrs, log):
        res = self.calc.Count(num_reads, num_tslrs)
        self.log.info("Total TSLRs: " + str(num_tslrs) + "(" + estToString(res.nonzero * 100) +
                      " % long reads were covered by short reads.")
        if res.nonzero < 0.7:
            self.log.info("WARNING: high fraction of uncovered long reads. Result is unreliable.")
        if res.nonzero > 0.01:
            log.info("Estimated metagenome size: " + estToString(res.est) + "+-" + estToString(res.disp))


class MetaLengthPipeline:
    def __init__(self, params):
        self.params = params

    def Run(self):
        # prepare to run
        io.ensure_dir_existence(self.params.output_dir)
        log_file = logging.StreamHandler(open(os.path.join(self.params.output_dir, "log.info"), "w"))
        self.params.log.addHandler(log_file)
        # find read count if not provided
        if self.params.read_count is None:
            self.params.log.info("Counting reads. To skip this step use --read-count option")
            self.params.read_count = self.CountReads()
            self.params.log.info("Number of reads: " + str(self.params.read_count))
        # preapare combined TSLRs file and construct index if not provided
        if self.params.tslr_index is not None or self.params.sam is not None:
            assert len(self.params.tslrs) == 1
            tslrs_file = self.params.tslrs[0]
        else:
            tslrs_file = self.ConcatTSLRs()
        # start alignment
        if self.params.sam is None:
            sam_handler = sam_parser.SamChain(map(sam_parser.Samfile, self.PerformAlignment(tslrs_file, self.params.log)))
        else:
            sam_handler = sam_parser.Samfile(open(self.params.sam, "r"))
        # calculate metagenome length
        self.ProcessSam(sam_handler, tslrs_file)
        # prepare to finish
        self.params.log.removeHandler(log_file)

    def PerformAlignment(self, tslrs_file, log):
        aligner = alignment.Bowtie2(self.params.bowtie_path, self.params.bowtie_params)
        if self.params.tslr_index is None:
            log.info("Creating index for TSLRs")
            alignment_calculator = alignment.AlignmentCalculator(self.params.output_dir, aligner, tslrs_file, log)
        else:
            alignment_calculator = alignment.AlignmentCalculator(self.params.output_dir, aligner, self.params.tslr_index, log)
        log.info("TSLR index ready")
        if self.params.save_sam:
            sam_files = alignment_calculator.align_bwa_pe_libs(zip(self.params.left_reads, self.params.right_reads), self.params.output_dir,
                                                               self.params.threads)
            sam_handler_list = [open(sam_file, "r") for sam_file in sam_files]
        else:
            sam_handler_list = [alignment_calculator.AlignPELibOnline(left, right, self.params.threads) for left, right in
                                zip(self.params.left_reads, self.params.right_reads)]
        return sam_handler_list

    def ProcessSam(self, sam_handler, tslrs_file):
        self.params.log.info("Preparing for estimation")
        is_counter = calculator.ISCounter()
        calc = calculator.Calculator(tslrs_file, self.params.min_len, is_counter, self.params.log)
        listeners = [is_counter, calc]
        printer = ResultPrinter(calc, self.params.log, self.params.debug)
        if self.params.debug:
            listeners.append(printer)
        self.params.log.info("Alignment analysis started")
        cnt = 0
        for rec in sam_handler:
            for listener in listeners:
                listener.Process(rec)
            cnt += 1
            if cnt % 10000000 == 0:
                self.params.log.info(str(cnt) + " alignments processed")
                self.params.log.info(rec.query_name)
        if calc.is_cnt.get() == 0:
            self.params.log.info("WARNING: Could not estimate insert size. Setting insert size value to 0.")
        self.params.log.info("Alignment analysis finished")
        self.params.log.info("\n")
        self.params.log.info("Insert size estimated as " + str(calc.is_cnt.get()))
        if self.params.output_coverages:
            self.print_coverages(calc.coverage_records, calc.seq_len, calc.is_cnt.get())
        printer.print_all_results(self.params.read_count)

    def print_coverages(self, coverage_records, seqs, ins):
        f = open(os.path.join(self.params.output_dir, "tslr_coverages.info"), "w")
        for rec, l in zip(coverage_records, seqs):
            if l > self.params.min_len:
                f.write(str(l) + " " + str(rec.get()) + " " + str(rec.get() / (l - ins)) + "\n")
        f.close()

    def CountReads(self):
        result = 0
        for f in self.params.left_reads:
            result += sum(1 for line in SeqIO.Open(f, "r"))
        return result / 2

    def ConcatTSLRs(self):
        file_name = os.path.join(self.params.output_dir, "tslrs.fasta")
        tslrs_file = open(file_name, "w")
        cnt = 0
        for f in self.params.tslrs:
            file_type = "fastq"
            if "fasta" in f.split(".") or "fa" in f.split("."):
                file_type = "fasta"
            for rec in SeqIO.parse(io.universal_read(f), file_type):
                if len(rec) > self.params.min_len:
                    rec.id = str(cnt) + "_" + rec.id
                    SeqIO.write(rec, tslrs_file, "fasta")
                    cnt += 1
        tslrs_file.close()
        return file_name
