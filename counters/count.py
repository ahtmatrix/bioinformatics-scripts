import sys
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import pairwise2
import warnings
from Bio import BiopythonWarning
import re


def count_records(fullpath, filename, filetype):

    num_records = 0
    for record in SeqIO.parse(fullpath, filetype):
        num_records +=1
    return num_records

# creates a list of the files in this directory
raw_datadir_listing = os.listdir(os.getcwd())

# loops over the list of files
for files in raw_datadir_listing:
    full_path = os.path.join(os.getcwd(), files)
    filename = os.path.splitext(files)[0]
    
    if files.endswith('.gbk'):
        print filename +" record count= "+ str(count_records(full_path, filename, "genbank"))
        
    elif files.endswith('.fasta'):
        print filename +" record count= "+ str(count_records(full_path, filename, "fasta"))