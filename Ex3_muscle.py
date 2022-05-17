from copyreg import constructor
from msilib.schema import Directory
from Bio import SeqIO

from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline

from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator , DistanceTreeConstructor
from Bio.Phylo.Consensus import bootstrap_trees , majority_consensus 

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

    output_folder = "output_ex3_muscle"
    if not os.path.exists(f"./{output_folder}"):
        os.makedirs(f"./{output_folder}")
    #first arg: aln filename
    combinedFile= f"./{output_folder}/{sys.argv[1]}.fasta"
    
    fasta_files = readFiles()

    #Combine all files and send it to muscle to align
    SeqIO.write( fasta_files, combinedFile,  "fasta")
    align_file = combinedFile[:-6] + '_aligned.fata' 
    cline=MuscleCommandline ( input=combinedFile , out=align_file)
    cline()


    #read alignment
  
    alignment=AlignIO.read( align_file ,  "fasta")
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
    Phylo.write(tree , f"./{output_folder}/tree.xml" , "phyloxml")
    print("Tree figure: ")
    Phylo.draw_ascii(tree)
    fig = plt.figure(figsize=(13,5) , dpi=100)
    axes = fig.add_subplot(1,1,1)
    Phylo.draw(tree ,axes=axes)
    fig.savefig(f"./{output_folder}/tree")

    #Boostrap
    trees = bootstrap_trees(alignment , 1000 ,constructor)
    tree = majority_consensus(trees)
    tree.rooted = True
    print("Bootstrapped Tree: ")
    print(tree)
    #Save tree into a file
    Phylo.write(tree , f"./{output_folder}/{sys.argv[1]}_bootstrap_tree.xml" , "phyloxml")
    print("Tree figure: ")
    Phylo.draw_ascii(tree)
    fig = plt.figure(figsize=(13,5) , dpi=100)
    axes = fig.add_subplot(1,1,1)
    Phylo.draw(tree ,axes=axes)
    fig.savefig(f"./{output_folder}/{sys.argv[1]}_bootstrap_tree")
    return
    return

if __name__ == "__main__":
    main()