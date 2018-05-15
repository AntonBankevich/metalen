#!/usr/bin/env python
import sys
sys.path.append("py")

from meta_length.metalen_io import create_log
from meta_length.params import MetaLengthParameters
from meta_length.pipeline import MetaLengthPipeline

if __name__ == "__main__":
    parameters = MetaLengthParameters(sys.argv, sys.stderr)
    log = create_log("meta_size")
    pipeline = MetaLengthPipeline(parameters, log)
    pipeline.Run()
