import sys
from metalen import sam_parser, SeqIO, metalen_io

class Record:
    def __init__(self, name, length, count = 0):
        self.name = name
        self.length = length
        self.count = count

def CollectRecords(fasta_file):
    res = dict()
    for rec in SeqIO.parse_fasta(metalen_io.universal_read(fasta_file)):
        res[rec.id] = Record(rec.id), len(rec)
    return res

if __name__ == "__main__":
    fasta_file = sys.argv[1]
    recs = CollectRecords(fasta_file)
    sam_file = sys.argv[2]
    for rec in sam_parser.Samfile(open(sam_file, "r")):
        if not rec.is_unmapped:
            recs[rec.tname].count += 1
    recs = recs.values()
    recs = sorted(recs, key = lambda rec: float(rec.count) / rec.length)[::-1]
    for rec in recs:
        print float(rec.count) / rec.length, rec.length, rec.count, rec.name