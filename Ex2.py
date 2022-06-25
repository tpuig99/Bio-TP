from Bio.Blast import NCBIWWW
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
import sys
import argparse

LOCAL_BLAST = 'local'
REMOTE_BLAST = 'remote'

parser = argparse.ArgumentParser()
parser.add_argument("input", metavar="INPUT", help="input file for BLAST, must be a FASTA file (.fas)", type=argparse.FileType('r'))
parser.add_argument("output", metavar="OUTPUT", help="output file for BLAST in XML format (.xml)")
parser.add_argument('-t', '--type', dest='blast_type', required=True, help="BLAST execution type, local or remote", choices=[LOCAL_BLAST, REMOTE_BLAST])
args = parser.parse_args()

with args.input as input:
    try:
        sequence_data = SeqIO.read(input, format="fasta")
        if not any(sequence_data):
            raise Exception()
    except:
        sys.exit(f"error: file {input.name} is empty or is not in FASTA format")

if (args.blast_type == LOCAL_BLAST):
    try:
        cline = NcbiblastnCommandline(query=sequence_data, db="nt", out=args.output, outfmt=5)
        cline()
    except:
        sys.exit("error: local BLAST failed, is blastn installed? If so, check env variables PATH and BLASTDB are correctly setted and that the local database is configured")

else:
    result_handle = NCBIWWW.qblast(
        program="blastn",
        database="nt",
        sequence=sequence_data.seq,
    )

    with open(args.output, 'w') as output_file:
        blast_results = result_handle.read()
        output_file.write(blast_results)