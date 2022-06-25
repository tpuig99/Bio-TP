import chunk
from cmath import inf
from itertools import count
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import sys

if len(sys.argv) != 3:
    print("Wrong arguments amount")
    exit(1)

parser = argparse.ArgumentParser()
parser.add_argument("input", metavar="INPUT", help="input file, must be a GENBANK file (.gb)")
parser.add_argument("output", metavar="OUTPUT", help="output file is a FASTA file with longest ORF(.fas)")
args = parser.parse_args()

try:
    f = open(args.input, "r")
except OSError:
    print("Could not open/read file:", args.input)
    exit()

with f as input_handle:

    try:
        record = SeqIO.read(input_handle, "genbank")
    except ValueError:
        print("Error reading as a genbank file")
        exit()
    

    table = 1
    min_pro_len = 100
    END_CODON =['TAG','TAA','TGA']
    START_CODON = "ATG"
    sequences = []
    sequences_f = []

    for strand, nuc in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
        for frame in range(3):
            length = 3 * ((len(record)-frame) // 3) #Multiple of three
            my_seq = Seq(nuc[frame:frame+length])
            l_start = 0
            l_end = 0
            codon_len = 3 # chunk length
            chunks = [my_seq[i:i+codon_len] for i in range(0, len(my_seq),codon_len)]
            
            start_codon = -1
            while True:
                try:
                    start_codon = chunks.index(START_CODON,start_codon+1)
                except ValueError:
                    break
                end_codon = inf
                for codon in END_CODON:
                    try:
                        n = chunks.index(codon,start_codon)
                    except ValueError:
                        continue
                    if n < end_codon:
                        end_codon = n
                if end_codon != inf and (end_codon - start_codon) > (l_end - l_start):
                    l_start = start_codon
                    l_end = end_codon

            l_start *= 3
            l_end *= 3
            sequences.append(my_seq[l_start:l_end+3])
    
    if len(sequences) == 0:
        print("None ORF where found on file")
        exit()

    maxSeq = max(sequences, key=lambda x:len(x))
    max_seq = SeqRecord(
        Seq(maxSeq),
        id=record.id,
        description=record.description
    )

    with open(args.output, "w") as output_handle:
        SeqIO.write(max_seq, output_handle, "fasta")
