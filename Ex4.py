from cmath import inf
from itertools import count
from Bio.Blast import NCBIXML
from Bio import Entrez
import sys
import argparse

if len(sys.argv) != 4:
    print("Wrong arguments amount")
    exit(1)

parser = argparse.ArgumentParser()
parser.add_argument("input", metavar="INPUT", help="input file, must be a BLAST file in XML format (.xml)")
parser.add_argument("pattern", metavar="PATTERN", help="pattern to look for")
parser.add_argument("output", metavar="OUTPUT", help="output file")
args = parser.parse_args()

try:
    f = open(args.input, "r")
except OSError:
    print("Could not open/read file:", args.input)
    exit()

pattern = args.pattern
out_file = open(args.output, "w")
Entrez.email = "@"

with f as input_handle:
     records= NCBIXML.parse(input_handle)
     i_hit = 1
     for rec in records:
          for hit in rec.alignments:
               if(pattern.lower() in hit.title.lower()):
                    i_hsp = 1
                    for hsp in hit.hsps:
                         out_file.write(f'~~~~~~Alignment {i_hit} - hsp {i_hsp}~~~~~~\n')
                         out_file.write('*sequence: ' + hit.title + '\n')
                         out_file.write('*accession: ' + hit.accession + '\n')
                         out_file.write('*length: ' + str(hit.length) + '\n')
                         out_file.write('*score: ' + str(hsp.score) + '\n')
                         out_file.write('*gaps: ' + str(hsp.gaps) + '\n')
                         out_file.write(Entrez.efetch(db="protein", rettype="fasta", id=hit.accession).read())
                         i_hsp += 1
                    i_hit += 1

out_file.close()