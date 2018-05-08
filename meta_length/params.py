# bowtie_params = namedtuple("tslr", "pacbio")("--fast --no-discordant --no-mixed -a", "--very-sensitive --no-discordant --no-mixed -a")
import getopt
import os

import metalen_io
import sys
import logging

bowtie_params = {"tslr" : "--fast --no-discordant --no-mixed -a".split(), "pacbio" : "-D 40 -R 3 -N 0 -L 19 -i S,1,0.50 --rdg 1,3 --rfg 1,3 -k 100 --score-min L,-0.6,-1 --ignore-quals".split()}
version = 1.0

class MetaLengthParameters:
    left_reads = []
    right_reads = []
    single_reads = []
    tslrs = []
    combined_tslrs = None
    tslr_index = None
    read_count = None
    output_dir = None
    threads = 8
    bowtie_path = None
    bowtie_params = bowtie_params["tslr"]
    save_sam = False
    debug = True
    output_coverages = False
    min_len = 6000
    log = None
    sam = None

    def print_usage_and_exit(self, code):
        # sys.stderr.write("min-len= output-coverages bowtie-path= tslrs= tslr-index= sam= save-sam read-count= threads= help output-dir= 1:2:ho:t:\n")
        sys.stderr.write("MetaLen v" + str(version) +
                         ": metagenome size estimator that uses combination of long and short reads.\n" +
                         "MetaLen supports Illumina paired-end reads and Illumina Synthetic Long reads\n\n")
        sys.stderr.write("Usage: " + str(sys.argv[0]) + " [options] -o <output_dir>" + "\n")
        sys.stderr.write("" + "\n")
        sys.stderr.write("Basic options:" + "\n")
        sys.stderr.write("-h/--help\t\t\tprints this usage message" + "\n")
        sys.stderr.write("-v/--version\t\t\tprints version" + "\n")
        # sys.stderr.write("--test\t\t\t\trun MetaLen on toy dataset" + "\n")
        sys.stderr.write("-o\t\t<output_dir>\tdirectory to store all the output (required)" + "\n")
        sys.stderr.write("-t/--threads\t<int>\t\tnumber of threads" + "\n")

        sys.stderr.write("" + "\n")
        sys.stderr.write("Input options:" + "\n")
        sys.stderr.write("-1\tfile with forward paired-end reads\n")
        sys.stderr.write("-2\tfile with reverse paired-end reads\n")
        sys.stderr.write("--long\tfile with long reads\n")

        sys.stderr.write("" + "\n")
        sys.stderr.write("Additional options:" + "\n")
        sys.stderr.write("--bowtie-path\tpath to location containing bowtie2 binary files. If this option is not specified we assume that bowtie2 binary directory is in system path\n")
        sys.stderr.write("--min-len\tminimal length of long read to be used for estimation\n")

        sys.exit(code)

    def __init__(self, argv):
        self.ParseValues(argv)
        self.CheckParamConsistency()

        self.log = logging.getLogger("meta_size")
        log = self.log

        log.setLevel(logging.DEBUG)
        console = logging.StreamHandler(sys.stdout)
        console.setFormatter(logging.Formatter('%(message)s'))
        console.setLevel(logging.DEBUG)
        log.addHandler(console)

        metalen_io.ensure_dir_existence(self.output_dir)
        file_log = logging.StreamHandler(open(os.path.join(self.output_dir, "meta_length.log"), "w"))
        file_log.setFormatter(logging.Formatter('%(message)s'))
        file_log.setLevel(logging.DEBUG)
        log.addHandler(file_log)

    def ParseValues(self, argv):
        if len(argv) == 1:
            self.print_usage_and_exit(1)
        long_params = "pacbio min-len= output-coverages bowtie-path= long= index= sam= save-sam read-count= threads= help output-dir=".split(" ")
        short_params = "s:1:2:ho:t:"
        try:
            options_list, tmp = getopt.gnu_getopt(argv[1:], short_params, long_params)
            if len(tmp) != 0:
                self.print_usage_and_exit(1)
        except getopt.GetoptError:
            _, exc, _ = sys.exc_info()
            sys.stderr.write("Parameter parsing error: " + str(exc) + "\n\n")
            self.print_usage_and_exit(1)
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
            elif key == "min-len":
                self.min_len = int(value)
            elif key == "--sam":
                self.sam = value
            elif key == "--read-count":
                self.read_count = int(value)
            elif key == "--threads" or key == "-t":
                self.threads = int(value)
            elif key == "--help" or key == "-h":
                self.print_usage_and_exit(0)
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
                sys.stderr.write("No such key: " + key + "\n")
                self.print_usage_and_exit(1)

    def CheckParamConsistency(self):
        if (self.sam and self.read_count) or (
                    len(self.left_reads) == len(self.right_reads) and len(self.left_reads) != 0 and (
                len(self.tslrs) != 0 or self.tslr_index)):
            pass
        else:
            sys.stderr.write("Wrong set of arguments")
            self.print_usage_and_exit(1)
        if self.output_dir is None:
            self.print_usage_and_exit(1)
        if len(self.tslrs) == 0:
            self.print_usage_and_exit(1)
        if self.tslr_index and len(self.tslrs) > 1:
            self.print_usage_and_exit()

