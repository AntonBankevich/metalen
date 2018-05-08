import sys
import time
from meta_length.params import MetaLengthParameters
from meta_length.pipeline import MetaLengthPipeline

if __name__ == "__main__":
    start = time.time()
    parameters = MetaLengthParameters(sys.argv)
    pipeline = MetaLengthPipeline(parameters)
    pipeline.Run()
    parameters.log.info("Finished in " + str(time.time() - start) + " seconds")
