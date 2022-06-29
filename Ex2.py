"""Executes blastp against a FASTA input and returns a result in XML format containing all matches

Help:
    python3 ./Ex2 -h

Example:
    pyhton3 ./Ex2.py input.fas output.xml --> runs remote blastp against input.fas with default parameters and saves result in output.xml
"""

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
parser.add_argument('-t', '--type', dest='blast_type', help="BLAST execution type, local or remote, default is remote", choices=[LOCAL_BLAST, REMOTE_BLAST], default=REMOTE_BLAST)
parser.add_argument('--database', '-d', dest='database', help="database aginst which BLAST will be executed, default is 'nt'", default="nt")
parser.add_argument('--expect', '-e', dest='expect', help="e-value cutoff, default is 10.0", default=10.0, type=float)
parser.add_argument('--hitlist_size', '-H', dest='hitlist_size', help="number of hits to return, ignored if BLAST is local, default is 50", default=50, type=int)
parser.add_argument('--word_size', '-w', dest='word_size', help="size of word for query, default is None", default=None, type=int)

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
        cline = NcbiblastnCommandline(
            query=sequence_data,
            db=args.database,
            evalue=args.expect,
            out=args.output,
            word_size=args.word_size,
            outfmt=5
        )
        cline()
    except:
        sys.exit("error: local BLAST failed, is blastn installed? If so, check env variables PATH and BLASTDB are correctly setted and that the local database is configured")

else:
    try:
        result_handle = NCBIWWW.qblast(
            program="blastn",
            database=args.database,
            sequence=sequence_data.seq,
            expect=args.expect,
            hitlist_size=args.hitlist_size,
            word_size=args.word_size, 
        )
    except:
        sys.exit("error: remote BLAST failed, check parameters and re-run program")

    with open(args.output, 'w') as output_file:
        blast_results = result_handle.read()
        output_file.write(blast_results)