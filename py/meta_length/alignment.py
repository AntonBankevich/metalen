############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os

from . import metalen_io


class AlignmentCalculator:
    #if reference is a .fasta of .fa file it will be ised to construct index. Otherwise it is considered as index itself
    def __init__(self, work_dir, aligner, reference, log):
        metalen_io.ensure_dir_existence(work_dir)
        self.work_dir = work_dir
        self.aligner = aligner
        self.log = log
        log.info("Constructing " + aligner.name + " index")
        self.index = os.path.join(work_dir, "index")
        self.log_file = os.path.join(work_dir, "output.log")
        self.err_log_file = os.path.join(work_dir, "err_output")
        if reference.endswith(".fasta") or reference.endswith(".fa"):
            command_line = aligner.IndexCommand(reference, self.index)
            log.info(" ".join(command_line))
            metalen_io.universal_sys_call(command_line, log, self.log_file, self.err_log_file)
        else:
            self.index = reference
        log.info("Index constructed.")

    def AlignPELib(self, reads1, reads2, output, threads):
        self.log.info("Aligning paired-end library")
        self.log.info("Left reads: " + reads1)
        self.log.info("Right reads: " + reads2)
        self.log.info("Output directory: " + self.work_dir)
        self.log.info("Starting read alignment. See detailed log in " + self.log_file)
        metalen_io.universal_sys_call(self.aligner.PairedAlignCommand(self.index, reads1, reads2, threads), self.log, output, self.err_log_file)
        self.log.info("Done. See result in " + output)
        return output

    def AlignPELibOnline(self, reads1, reads2, threads):
        self.log.info("Aligning paired-end library")
        self.log.info("Left reads: " + reads1)
        self.log.info("Right reads: " + reads2)
        self.log.info("Output directory: " + self.work_dir)
        self.log.info("Starting read alignment. See detailed log in " + self.log_file)
        return metalen_io.online_sys_call(self.aligner.PairedAlignCommand(self.index, reads1, reads2, threads), self.log_file)

    def align_bwa_pe_libs(self, reads, output_dir, threads):
        self.log.info("===== Starting read alignment")
        result = []
        lib_num = 1
        metalen_io.ensure_dir_existence(output_dir)
        for left_reads, right_reads in reads:
            result.append(self.AlignPELib(left_reads, right_reads, os.path.join(output_dir, str(lib_num) + ".sam"), threads))
            lib_num += 1
        self.log.info("===== Read alignment finished. See result in " + output_dir)
        return result

class BWA:
    def __init__(self, path = None, algorithm = "is", additional_params = []):
        self.name = "bwa"
        self.algorithm = algorithm
        self.additional_params = additional_params
        self.align_command = "bwa mem"
        self.index_command = "bwa index"
        if path is not None:
            self.align_command = path + " mem"
            self.index_command = path + " index"

    def IndexCommand(self, reference, output):
        return " ".join([self.index_command, "-p", output, "-a", self.algorithm, reference]).split()

    def PairedAlignCommand(self, index, reads1, reads2, threads):
        return " ".join([self.align_command, "-t", str(threads), "-S", "-M", " ".join(self.additional_params), index, reads1, reads2]).split()

class Bowtie2:
    def __init__(self, path=None, additional_params = ["--fast", "--no-discordant", "--no-mixed"]):
        self.name = "bowtie2"
        self.align_command = "bowtie2"
        self.index_command = "bowtie2-build"
        self.additional_params = additional_params
        if path is not None:
            self.align_command = os.path.join(path, self.align_command)
            self.index_command = os.path.join(path, self.index_command)

    def IndexCommand(self, reference, output):
        return " ".join([self.index_command, reference, output]).split(" ")

    def PairedAlignCommand(self, index, reads1, reads2, threads):
        return " ".join([self.align_command, "-p", str(threads), "-x", index, "-1", reads1, "-2", reads2, " ".join(self.additional_params)])

