import sys
import os
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import pairwise2

# Usage
# python datavalidation.py

def get_upstream_cds(fullpath, filename):
    
    extracted_cds_list = []
    # reads in a gbk and creates a Seqrecordord object
    for record in SeqIO.parse(fullpath, "genbank"):
        if record.features:
            for feature in record.features:
                if feature.type == "CDS":
                    
                    
                    translated_protein = str(feature.qualifiers.get('translation', 'no_translation')).strip('\'[]')
                    cds_to_protein = str(feature.extract(record).seq.translate(to_stop = True))
                    
                    
                    if translated_protein != cds_to_protein:
                        
                        alignment = pairwise2.align.globalxx(translated_protein, cds_to_protein, one_alignment_only=True)
                        print(pairwise2.format_alignment(*alignment[0]))
                        
    return

# creates a list of the files in this directory
raw_datadir_listing = os.listdir(os.getcwd())
# loops over the list of files
for files in raw_datadir_listing:
    if files.endswith('.gbk'):
        full_path = os.path.join(os.getcwd(), files)
        filename = os.path.splitext(files)[0]
        get_upstream_cds(full_path, filename)
