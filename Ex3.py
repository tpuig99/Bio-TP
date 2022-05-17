from copyreg import constructor
from Bio import SeqIO

from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline

from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator , DistanceTreeConstructor

import matplotlib
import matplotlib.pyplot as plt
import glob
import sys
import re
import os

#Read all fasta files
def readFiles():
    fasta_filenames = glob.glob("./fasta_inputs/*")
    fasta_files = [ SeqIO.read(f, format="fasta") for f in fasta_filenames]
    for idx , s in enumerate( fasta_files):
        print(os.path.basename(fasta_filenames[idx])[:-6])
        s.id = os.path.basename(fasta_filenames[idx])[:-6]
    return fasta_files

def main():
    if not os.path.exists("./output_ex3"):
        os.makedirs("./output_ex3")
    #first arg: aln filename
    combinedFile= "./output_ex3/" +sys.argv[1] + ".fasta"
    fasta_files = readFiles()

    #Combine all files and send it to clustal to align
    SeqIO.write( fasta_files, combinedFile,  "fasta")
    cline=ClustalwCommandline("clustalw2" , infile=combinedFile)
    cline()


    #read alignment
    align_file = combinedFile[:-6] + '.aln' 
    alignment=AlignIO.read( align_file ,  "clustal")
    print("Alignment: ")
    print(alignment)

    #create distance matrix
    calculator = DistanceCalculator('identity')
    matrix = calculator.get_distance(alignment)
    print("Distance matrix: ")
    print(matrix)
    matrix_filepath = f"{combinedFile[:-6]}_matrix.txt"
    with open(matrix_filepath, "w") as matrix_file:
        matrix_file.write(str(matrix))

    
    #create tree
    constructor= DistanceTreeConstructor(calculator)
    tree = constructor.build_tree(alignment)
    tree.rooted = True
    print("Tree: ")
    print(tree)
    #Save tree into a file
    Phylo.write(tree , "output_ex3/ex3_tree.xml" , "phyloxml")
    print("Tree figure: ")
    Phylo.draw_ascii(tree)
    fig = plt.figure(figsize=(13,5) , dpi=100)
    axes = fig.add_subplot(1,1,1)
    Phylo.draw(tree ,axes=axes)
    fig.savefig("./output_ex3/tree")
    return

if __name__ == "__main__":
    main()