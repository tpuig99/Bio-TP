from cmath import inf
from itertools import count
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys


with open(sys.argv[1], "rU") as input_handle:

    record = SeqIO.read(input_handle, "genbank")

    table = 1
    min_pro_len = 100
    END_CODON =['TAG','TAA','TGA']
    sequences = []
    sequences_f = []

    for strand, nuc in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
        for frame in range(3):
            length = 3 * ((len(record)-frame) // 3) #Multiple of three
            my_seq = Seq(nuc[frame:frame+length])
            l_start = 0
            l_end = 0
            f_start = 0
            f_end = 0
            start_codon = my_seq.find("ATG")
            while start_codon  != -1:
                end_codon = inf
                for codon in END_CODON:
                    n = my_seq.find(codon,start_codon+2)
                    if n < end_codon:
                        end_codon = n
                if (end_codon - start_codon) > (l_end - l_start):
                    l_start = start_codon
                    l_end = end_codon
                if f_end == 0:
                    f_start = start_codon
                    f_end = end_codon
                start_codon = my_seq.find("ATG",start_codon+2)
            sequences.append(my_seq[l_start:l_end+3])
            sequences_f.append(my_seq[f_start:f_end+3])

    maxSeq = max(sequences, key=lambda x:len(x))
    max_seq = SeqRecord(
        Seq(maxSeq),
        id=record.id,
        description=record.description
    )

    maxSeqF = max(sequences_f, key=lambda x:len(x))
    max_seq_f = SeqRecord(
            Seq(maxSeqF),
            id=record.id,
            description=record.description
    )
    with open(sys.argv[2], "w") as output_handle:
        SeqIO.write(max_seq, output_handle, "fasta")
    with open(f'first_{sys.argv[2]}', "w") as output_handle:
        SeqIO.write(max_seq_f, output_handle, "fasta")