from Bio import SeqIO
import sys
import os

filename = sys.argv[1]

count = 0

if filename.endswith('.gbk'):
    filetype = "genbank"
elif filename.endswith('.fasta'):
    filetype = "fasta"

for record in SeqIO.parse(filename, filetype):
    count = count + 1


print("There were " + str(count) + " records in file " + filename)
