########################################################
# Process R1 for cellranger-atac pipeline
########################################################

from Bio.SeqIO.QualityIO import FastqGeneralIterator
from gzip import open as gzopen

import argparse

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input_R1F", required=True, help="input file")
ap.add_argument("-i2", "--input_R1R", required=True, help="input file")
ap.add_argument("-i3", "--input_R2F", required=True, help="input file")
ap.add_argument("-i4", "--input_R2R", required=True, help="input file")
ap.add_argument("-o", "--output", required=True, help="output barcode file")
args = vars(ap.parse_args())

input_file_R1F = args["input_R1F"]
input_file_R1R = args["input_R1R"]
input_file_R2F = args["input_R2F"]
input_file_R2R = args["input_R2R"]

output_file = args["output"]

with open(input_file_R1F, "r") as in_handle_R1F, open(input_file_R1R, "r") as in_handle_R1R, open(input_file_R2F, "r") as in_handle_R2F, open(input_file_R2R, "r") as in_handle_R2R, open(output_file, "w") as out_handle:
    for (title, seq, qual),(title2, seq2, qual2),(title3, seq3, qual3),(title4, seq4, qual4), in zip(FastqGeneralIterator(in_handle_R1F),FastqGeneralIterator(in_handle_R1R),FastqGeneralIterator(in_handle_R2F),FastqGeneralIterator(in_handle_R2R)):
        new_seq = seq + seq2 + seq3 +seq4
        new_qual = qual + qual2 + qual3 + qual4       
        out_handle.write("@%s\n%s\n+\n%s\n" % (title, new_seq, new_qual))
