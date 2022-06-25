"""Executes Emboss to get all ORFs from a sequence and then finds motifs using PROSITE database

[MUST] Download Emboss and configure PROSITE:
    1. (for Ubuntu) sudo apt-get install emboss
    2. Download `prosite.dat` and `prosite.doc` from https://ftp.expasy.org/databases/prosite/ and save them int a directory (e.g. $HOME/prosite)
    3. sudo prosextract -prositedir $PATH_TO_PROSITE_DIR (this last param should be chosen dir from the previous step)

Help:
    ./Ex5 -h

Example:
    ./Ex5.py input.fas output.patmatmotifs --> runs emboss getorf and then finds motifs with patmatmotifs to generate a report 
"""

import argparse
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument("input", metavar="INPUT", help="input file for 'getorf', must be a GENBANK file (.gb)", type=argparse.FileType('r'))
parser.add_argument("output_getorf", metavar="OUTPUT_GETORF", help="intermediate output file for 'getorf' (.fas)")
parser.add_argument("output_patmatmotifs", metavar="OUTPUT_PATMATMOTIFS", help="final output file for 'patmatmotifs' (.patmatmotifs)")
parser.add_argument('-f', '--full', dest='full', help="executes a full (documented) output for matpatmotifs, default true", default=True)
parser.add_argument('-np', '--noprune', dest='noprune', help="allows simple motifs to be included in output, default true", default=True)
parser.add_argument('-min', '--minsize', dest='minsize', help="min size of nt ORF to be reported, default is 30", type=int, default=30)
parser.add_argument('-max', '--maxsize', dest='maxsize', help="max size of nt ORF to be reported, default is 1000000", type=int, default=1000000)

args = parser.parse_args()

getorf = f"getorf -minsize {args.minsize} -maxsize {args.maxsize} {args.input.name} {args.output_getorf}"
patmatmotifs = f"patmatmotifs {'-full' if args.full else ''} {'-noprune' if args.noprune else ''} {args.output_getorf} {args.output_patmatmotifs}"

print('Executing getorf with following command:', getorf)
if (os.system(getorf) != 0):
    sys.exit('error: failed execution of getorf, please check parameters and input file')

print('Executing patmatmotifs with following command:', getorf)
if (os.system(patmatmotifs) != 0):
    sys.exit('error: failed execution of patmatmotifs, please check parameters')
