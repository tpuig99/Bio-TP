from Bio.Blast import NCBIWWW
from Bio import SeqIO
import sys

sequence_data = SeqIO.read(sys.argv[1], format="fasta")

result_handle = NCBIWWW.qblast(
    program="blastn",
    database="nt",
    sequence=sequence_data.seq,
)

with open(sys.argv[2], 'w') as save_file:
    blast_results = result_handle.read()
    save_file.write(blast_results)