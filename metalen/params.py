# bowtie_params = namedtuple("tslr", "pacbio")("--fast --no-discordant --no-mixed -a", "--very-sensitive --no-discordant --no-mixed -a")
import getopt
import os

import io
import sys
import logging

bowtie_params = {"tslr" : "--fast --no-discordant --no-mixed -a".split(), "pacbio" : "-D 40 -R 3 -N 0 -L 19 -i S,1,0.50 --rdg 1,3 --rfg 1,3 -k 100 --score-min L,-0.6,-1 --ignore-quals".split()}

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
        sys.stderr.write("min-len= output-coverages bowtie-path= tslrs= tslr-index= sam= save-sam read-count= threads= help output-dir= 1:2:ho:t:\n")
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
        io.ensure_dir_existence(self.output_dir)
        console = logging.StreamHandler(open(os.path.join(self.output_dir, "meta_length.log"), "w"))
        console.setFormatter(logging.Formatter('%(message)s'))
        console.setLevel(logging.DEBUG)
        log.addHandler(console)

    def ParseValues(self, argv):
        if len(argv) == 1:
            self.print_usage_and_exit(1)
        long_params = "pacbio min-len= output-coverages bowtie-path= tslrs= tslr-index= sam= save-sam read-count= threads= help output-dir=".split(" ")
        short_params = "s:1:2:ho:t:"
        try:
            options_list, tmp = getopt.gnu_getopt(argv[1:], short_params, long_params)
            if len(tmp) != 0:
                self.print_usage_and_exit(1)
        except getopt.GetoptError:
            _, exc, _ = sys.exc_info()
            sys.stderr.write(str(exc) + "\n")
            self.print_usage_and_exit(1)
        for (key, value) in options_list:
            if key == "--bowtie-path":
                self.bowtie_path = value
            elif key == "--pacbio":
                self.bowtie_params = bowtie_params["pacbio"]
            elif key == "--tslrs":
                self.tslrs.append(value)
            elif key == "--tslr-index":
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

