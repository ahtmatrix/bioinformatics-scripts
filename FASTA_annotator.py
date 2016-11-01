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

def set_viral_TIS(fullpath, filename):

    annotated_extracted_TIS_list = []

    for record in SeqIO.parse(fullpath, "fasta"):
        record_description = str(record.description)
        annotated_record = SeqRecord(record.seq, record.name,description = record_description + str(csv_dictionary.get(record.id)))

        annotated_extracted_TIS_list.append(annotated_record)
        
    SeqIO.write(annotated_extracted_TIS_list, "annotated_" + filename +".fasta", "fasta")
    return

def set_rna_TIS(fullpath, filename):

    annotated_extracted_TIS_list = []

    for record in SeqIO.parse(fullpath, "fasta"):
        record_description = str(record.description)
                                                                                                #you can change this
        annotated_record = SeqRecord(record.seq, record.name,description = record_description + "Scanning")

        annotated_extracted_TIS_list.append(annotated_record)
        
    SeqIO.write(annotated_extracted_TIS_list, "annotated_" + filename +".fasta", "fasta")
    return

# creates a list of the files in this directory
raw_datadir_listing = os.listdir(os.getcwd())

# loops over the list of files
for files in raw_datadir_listing:
    if files.startswith('extracted_TIS_viral'):
        full_path = os.path.join(os.getcwd(), files)
        filename = os.path.splitext(files)[0]
        set_viral_TIS(full_path, filename)
    elif files.startswith('extracted_TIS_rna'):
        full_path = os.path.join(os.getcwd(), files)
        filename = os.path.splitext(files)[0]
        set_rna_TIS(full_path, filename)