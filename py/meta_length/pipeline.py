import logging
import math
import os
from logging import Logger
from typing import BinaryIO, List, Tuple

import histogram
from calculator import LongReadRecord, Calculator, ISCounter, ShiftCounter
from params import MetaLengthParameters
from sam_parser import Samfile
from . import SeqIO
from . import sam_parser
from . import alignment
from . import metalen_io
from . import calculator
import time

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
    def __init__(self, calc, is_calc, log, dir, debug):
        # type: (Calculator, ShiftCounter, Logger, str, bool) -> None
        self.calc = calc
        self.is_calc = is_calc
        self.log = log
        self.debug = debug
        self.limits, self.tslr_limits = self.generate_output_params(calc.tslr_count)
        print self.limits
        print self.tslr_limits
        self.cur_limit = 0
        self.dir = dir
        self.debug = debug

    def generate_output_params(self, tslr_count):
        from meta_length.calculator import LimitSequence
        if self.debug:
            return list(LimitSequence(1000)), LimitSequence(10, self.calc.tslr_count)
        else:
            return [], LimitSequence(tslr_count, self.calc.tslr_count)

    def Process(self, rec):
#        if rec.query_name != "" and GetReadNum(rec.query_name) > self.limits[self.cur_limit]:
        if rec.query_name != "" and self.calc.read_number == self.limits[self.cur_limit]:
            self.Point(str(self.calc.read_number))

    def Point(self, name):
        # type: (str) -> None
        self.log.info("Printing results for " + str(self.calc.read_number) + " reads")
        self.print_all_results()
        self.cur_limit += 1
        if histogram.Ready:
            hfname = os.path.join(self.dir, "histogram_" + str(self.calc.read_number) + ".pdf")
            self.calc.draw(hfname, self.is_calc.get())
            self.log.info("Histogram written to " + hfname)
        else:
            self.log.info("WARNING: can not draw histogram. " + histogram.Error)
        hfname = os.path.join(self.dir, "heights_" + str(self.calc.read_number) + ".txt")
        out = open(hfname, "w")
        self.calc.dump(out, self.is_calc.get())
        out.close()

    def print_all_results(self):
        from meta_length.calculator import LimitSequence
        tslr_limits = list(LimitSequence(self.tslr_limits, self.calc.tslr_count))
        self.log.info("")
        self.log.info("Results for " + str(self.calc.read_number) + " short reads")
        for num_tslrs in tslr_limits:
            self.print_results(num_tslrs)
        self.log.info("")

    def print_results(self, num_tslrs = None):
        if num_tslrs is None:
            num_tslrs = self.calc.tslr_count
        res = self.calc.Count(self.is_calc.get(), num_tslrs)
        self.log.info("Total TSLRs: " + str(num_tslrs) + ". " + estToString(res.nonzero * 100) +
                      " % long reads were covered by short reads.")
        if res.nonzero < 0.7:
            self.log.info("WARNING: high fraction of uncovered long reads. Result is unreliable.")
        if res.nonzero > 0.01:
            self.log.info("Estimated metagenome size: " + estToString(res.est) + "+-" + estToString(res.disp))


class MetaLengthPipeline:
    def __init__(self, params, log):
        # type: (MetaLengthParameters, Logger) -> None
        self.params = params
        self.log = log

    def Run(self):
        # prepare to run
        start = time.time()
        metalen_io.ensure_dir_existence(self.params.output_dir)
        alignment_dir = os.path.join(self.params.output_dir, "alignment")
        metalen_io.ensure_dir_existence(alignment_dir)
        log_file = logging.StreamHandler(open(os.path.join(self.params.output_dir, "meta_len.log"), "w"))
        self.log.addHandler(log_file)
        # write input params to file
        param_file = open(os.path.join(self.params.output_dir, "params.txt"), "w")
        param_file.write("\t".join(self.params.input))
        param_file.close()
        # prepare combined TSLRs file and construct index if not provided
        if self.params.tslr_index is not None or self.params.sam is not None:
            assert len(self.params.tslrs) == 1
            tslrs_file = self.params.tslrs[0]
        else:
            long_read_file = os.path.join(alignment_dir, "long.fasta")
            tslrs_file = self.ConcatTSLRs(long_read_file)
        # start alignment
        if self.params.sam is None:
            sam_handler = sam_parser.SamChain(map(sam_parser.Samfile, self.PerformAlignment(tslrs_file, alignment_dir)))
        else:
            sam_handler = sam_parser.Samfile(open(self.params.sam, "r"))
        # calculate metagenome length
        self.ProcessSam(sam_handler, tslrs_file)
        # prepare to finish
        self.log.info("Finished in " + str(time.time() - start) + " seconds\n")
        self.log.removeHandler(log_file)

    def PerformAlignment(self, tslrs_file, alignment_dir):
        # type: (str, str) -> List[BinaryIO]
        aligner = alignment.Bowtie2(self.params.bowtie_path, self.params.bowtie_params)
        metalen_io.ensure_dir_existence(alignment_dir)
        if self.params.tslr_index is None:
            self.log.info("Creating index for TSLRs")
            alignment_calculator = alignment.AlignmentCalculator(alignment_dir, aligner, tslrs_file, self.log)
        else:
            alignment_calculator = alignment.AlignmentCalculator(alignment_dir, aligner, self.params.tslr_index, self.log)
            self.log.info("TSLR index ready")
        if self.params.save_sam:
            self.log.info("Starting alignment")
            sam_files = alignment_calculator.align_bwa_pe_libs(zip(self.params.left_reads, self.params.right_reads), alignment_dir,
                                                               self.params.threads)
            sam_handler_list = [open(sam_file, "r") for sam_file in sam_files]
            self.log.info("Finished alignment")
        else:
            self.log.info("Starting alignment")
            sam_handler_list = [alignment_calculator.AlignPELibOnline(left, right, self.params.threads) for left, right in
                                zip(self.params.left_reads, self.params.right_reads)]
        return sam_handler_list

    def ProcessSam(self, sam_handler, tslrs_file):
        # type: (Samfile, str) -> None
        self.log.info("Preparing for estimation")
        is_counter = calculator.ISCounter()
        self.log.info("Reading TSLRs")
        calc = calculator.Calculator(tslrs_file, self.params.min_len)
        listeners = [is_counter, calc]
        printer = ResultPrinter(calc, is_counter, self.log, self.params.output_dir, self.params.debug)
        if self.params.debug:
            listeners.append(printer)
        print str(self.params.debug)
        self.log.info("Alignment analysis started")
        cnt = 0
        for rec in sam_handler:
            for listener in listeners:
                listener.Process(rec)
            cnt += 1
            if cnt % 10000000 == 0:
                self.log.info(str(cnt) + " alignments processed")
                self.log.info(rec.query_name)
        self.log.info("Alignment analysis finished")
        printer.Point("final")
        if is_counter.get() == 0:
            self.log.info("WARNING: Could not estimate insert size. Setting insert size value to 0.")
        else:
            self.log.info("Insert size estimated as " + str(is_counter.get()))
        if self.params.output_coverages:
            self.print_coverages(calc.coverage_records, is_counter.get(), calc.read_number)
        if self.params.debug:
            printer.print_all_results()
        else:
            printer.print_results()

    def print_coverages(self, coverage_records, ins, read_number):
        f = open(os.path.join(self.params.output_dir, "long_read_coverages.info"), "w")
        for rec in coverage_records:
            l = len(rec)
            if l > self.params.min_len:
                f.write(rec.id + " " + str(l) + " " + str(rec.get()) + " " + str(rec.get() / (l - ins) / read_number) + "\n")
        f.close()

    def CountReads(self):
        result = 0
        for f in self.params.left_reads:
            result += sum(1 for line in metalen_io.universal_open(f, "r")[0])
        return result / 2

    def ConcatTSLRs(self, file_name):
        tslrs_file = open(file_name, "w")
        cnt = 0
        for f in self.params.tslrs:
            handler, name = metalen_io.universal_open(f, "r")
            for rec in SeqIO.parse(handler, name):
                if len(rec) > self.params.min_len:
                    rec.id = str(cnt) + "_" + rec.id
                    SeqIO.write(rec, tslrs_file, "fasta")
                    cnt += 1
        tslrs_file.close()
        return file_name
