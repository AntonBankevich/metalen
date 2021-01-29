# bowtie_params = namedtuple("tslr", "pacbio")("--fast --no-discordant --no-mixed -a", "--very-sensitive --no-discordant --no-mixed -a")
import getopt
import os

import sys
from typing import List, BinaryIO

bowtie_params = {"tslr" : "--fast --no-discordant --no-mixed -a".split(), "pacbio" : "-D 40 -R 3 -N 0 -L 19 -i S,1,0.50 --rdg 1,3 --rfg 1,3 -k 100 --score-min L,-0.6,-1 --ignore-quals".split()}
version = 1.0

class MetaLengthParameters:
    left_reads = []
    right_reads = []
    single_reads = []
    tslrs = []
    combined_tslrs = None
    tslr_index = None
    hidden_read_count = None
    output_dir = None
    threads = 32
    bowtie_path = None
    bowtie_params = bowtie_params["tslr"]
    save_sam = False
    debug = True
    output_coverages = True
    min_len = 6000
    sam = None

    def print_usage_and_exit(self, err_log, code):
        err_log.write("MetaLen v" + str(version) +
                         ": metagenome size estimator that uses combination of long and short reads.\n" +
                         "MetaLen supports Illumina paired-end reads and Illumina Synthetic Long reads\n\n")
        err_log.write("Usage: " + str(self.program_name) + " [options] -o <output_dir>" + "\n")
        err_log.write("" + "\n")
        err_log.write("Basic options:" + "\n")
        err_log.write("-h/--help\t\t\tprints this usage message" + "\n")
        err_log.write("-v/--version\t\t\tprints version" + "\n")
        # err_log.write("--test\t\t\t\trun MetaLen on toy dataset" + "\n")
        err_log.write("-o/--output-dir\t\t\tdirectory to store all the output (required)" + "\n")
        err_log.write("-t/--threads\t<int>\t\tnumber of threads" + "\n")

        err_log.write("" + "\n")
        err_log.write("Input options:" + "\n")
        err_log.write("-1\tfile with forward paired-end reads\n")
        err_log.write("-2\tfile with reverse paired-end reads\n")
        err_log.write("--long\tfile with long reads\n")

        err_log.write("" + "\n")
        err_log.write("Additional options:" + "\n")
        err_log.write("--bowtie-path\tpath to location containing bowtie2 binary files. If this option is not specified we assume that bowtie2 binary directory is in system path\n")
        err_log.write("--min-len\tminimal length of long read to be used for estimation\n")

        sys.exit(code)

    def __init__(self, argv, err_log):
        # type: (List[str], BinaryIO) -> None
        self.program_name = argv[0]
        self.ParseValues(argv, err_log)
        self.CheckParamConsistency(err_log)

        # self.log = logging.getLogger("meta_size")
        # log = self.log
        #
        # log.setLevel(logging.DEBUG)
        # console = logging.StreamHandler(sys.stdout)
        # console.setFormatter(logging.Formatter('%(message)s'))
        # console.setLevel(logging.DEBUG)
        # log.addHandler(console)
        #
        # metalen_io.ensure_dir_existence(self.output_dir)
        # file_log = logging.StreamHandler(open(os.path.join(self.output_dir, "meta_len.log"), "w"))
        # file_log.setFormatter(logging.Formatter('%(message)s'))
        # file_log.setLevel(logging.DEBUG)
        # log.addHandler(file_log)

    def ParseValues(self, argv, err_log,):
        self.input = argv
        if len(argv) == 1:
            self.print_usage_and_exit(err_log, 1)
        long_params = "pacbio min-len= output-coverages bowtie-path= long= index= sam= save-sam read-count= threads= help output-dir=".split(" ")
        short_params = "s:1:2:ho:t:"
        try:
            options_list, tmp = getopt.gnu_getopt(argv[1:], short_params, long_params)
            if len(tmp) != 0:
                err_log.write("Parameter parsing error: unrecognized parameter " + str(tmp[0]) + "\n\n")
                self.print_usage_and_exit(err_log, 1)
        except getopt.GetoptError as e:
            err_log.write("Parameter parsing error: " + str(e) + "\n\n")
            self.print_usage_and_exit(err_log, 1)
        for (key, value) in options_list:
            if key == "--bowtie-path":
                self.bowtie_path = value
            elif key == "--pacbio":
                self.bowtie_params = bowtie_params["pacbio"]
            elif key == "--long":
                self.tslrs.append(value)
            elif key == "--index":
                self.tslr_index = value
            elif key == "--output-coverages":
                self.output_coverages = True
            elif key == "--min-len":
                self.min_len = int(value)
            elif key == "--sam":
                self.sam = value
            elif key == "--read-count":
                self.hidden_read_count = int(value)
            elif key == "--threads" or key == "-t":
                self.threads = int(value)
            elif key == "--help" or key == "-h":
                self.print_usage_and_exit(err_log, 0)
            elif key == "--output-dir" or key == "-o":
                self.output_dir = value
            elif key == "--save-sam":
                self.save_sam = True
            elif key == "-1":
                self.left_reads.append(value)
            elif key == "-2":
                self.right_reads.append(value)
            elif key == "-s":
                self.single_reads.append(value)
            else:
                err_log.write("No such key: " + key + "\n")
                self.print_usage_and_exit(err_log, 1)

    def CheckParamConsistency(self, err_log):
        if self.output_dir is None:
            err_log.write("Parameter error: no output directory provided\n\n")
            self.print_usage_and_exit(err_log, 1)
            return
        if len(self.left_reads) != len(self.right_reads):
            err_log.write("Parameter error: different number of files with left and right reads\n\n")
            self.print_usage_and_exit(err_log, 1)
            return
        if self.sam is not None and len(self.left_reads) == 0:
            err_log.write("Parameter error: at least one short read library should be provided\n\n")
            self.print_usage_and_exit(err_log, 1)
            return
        if len(self.tslrs) == 0:
            err_log.write("Parameter error: at least one long read library should be provided\n\n")
            self.print_usage_and_exit(err_log, 1)
            return
        if self.tslr_index and len(self.tslrs) > 1:
            err_log.write("Parameter error: long read index should be constructed from a single long read file filtered by read length\n\n")
            self.print_usage_and_exit(err_log, 1)

