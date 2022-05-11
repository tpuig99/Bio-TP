from Bio import SeqIO
import sys

with open(sys.argv[1], "rU") as input_handle:
    with open(sys.argv[2], "w") as output_handle:
        sequences = SeqIO.parse(input_handle, "genbank")
        count = SeqIO.write(sequences, output_handle, "fasta")

print("Converted %i records" % count)