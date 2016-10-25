import sys
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import pairwise2
import warnings
from Bio import BiopythonWarning
import re
import csv

# make a dictionary of the csv file
reader = csv.DictReader(open('TIS_annotation.csv'))
csv_dictionary = {}
for row in reader:
    csv_dictionary[row['Accession_#']] = row['primary_TIS_type']


print csv_dictionary.get("NC_005179.1")


def set_TIS(fullpath, filename):

    annotated_extracted_TIS_list = []

    for record in SeqIO.parse(fullpath, "fasta"):
        
        try:
            
            print record.id + " : " + csv_dictionary.get(record.id)
        except TypeError:
        
            print record.id
    # SeqIO.write(annotated_extracted_TIS_list, "annotated_extracted_TIS_" +
    #             filename + ".TIS.fasta", "fasta")
    return

# creates a list of the files in this directory
raw_datadir_listing = os.listdir(os.getcwd())

# loops over the list of files
for files in raw_datadir_listing:
    if files.startswith('extracted_TIS'):
        full_path = os.path.join(os.getcwd(), files)
        filename = os.path.splitext(files)[0]
        set_TIS(full_path, filename)