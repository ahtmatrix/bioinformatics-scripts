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


# dictionary[new_key] = dictionary[old_key]
# del dictionary[old_key]


def fix_accession(fullpath, filename):

    for record in SeqIO.parse(fullpath, "genbank"):
    
        for k, v  in csv_dictionary.items():
            if k == record.id:
                continue
            elif k == record.name:
                csv_dictionary[record.id] = csv_dictionary[k]
                del csv_dictionary[k]
            
    return

# creates a list of the files in this directory
raw_datadir_listing = os.listdir(os.getcwd())

# loops over the list of files
for files in raw_datadir_listing:
    if files.endswith(".gbk"):
        full_path = os.path.join(os.getcwd(), files)
        filename = os.path.splitext(files)[0]
        fix_accession(full_path, filename)


with open('fixed_TIS_annotation.csv', 'wb') as csv_file:
    writer = csv.writer(csv_file)
    for k,v in csv_dictionary.items():
        writer.writerow([k,v])
