import sys
import os
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import pairwise2
import warnings
from Bio import BiopythonWarning

# Usage
# python ExtractUpstreamCDS.py [number of bases upstream to cut]

# to combine multiple .gbk files into 1 gbk
# navigate to directory containing gbk files
# cat *.gbk > filename.gbk

# python -m pdb -Werror myprogram.py
warnings.filterwarnings('error')
num_bp_upstreamcds = int(sys.argv[1])

def validate_cds(record, feature):
    
    try:
        
        protein_in_file = str(feature.qualifiers.get('translation', 'no_translation')).strip('\'[]')
        
        cds_to_protein = str(feature.extract(record).seq.translate(to_stop = True))
    except BiopythonWarning:
        print record.id
        
    
        
    
    # if protein_in_file != cds_to_protein:
    #     print "oh snap"
    #     #alignment = pairwise2.align.globalxx(translated_protein, cds_to_protein, one_alignment_only=True)
        #print(pairwise2.format_alignment(*alignment[0]))
    return


def get_upstream_cds(fullpath, filename):

    extracted_cds_list = []
    
    # reads in a gbk and creates a Seqrecordord object
    for record in SeqIO.parse(fullpath, "genbank"):
        
        if record.features:
            for feature in record.features:
                if feature.type == "CDS":
                    
                    validate_cds(record, feature)
                    
                    #get the CDS nucleotide locations
                    cds_location = feature.location
                    start = cds_location.start.position
                    end = cds_location.end.position
            
                    
                    #create a SeqFeature object containing the location of where to extract
                    upstream_cds = SeqFeature(FeatureLocation(start-num_bp_upstreamcds, start))
                    
                    extracted_upstream_cds = upstream_cds.extract(record)
                    
                    if len(extracted_upstream_cds.seq) != num_bp_upstreamcds:
                        asldkfjalksdjflasjdf = 0
                        #print "upstream cds length is too short = " + str(len(extracted_upstream_cds.seq))
                    else:
                        extracted_cds_list.append(extracted_upstream_cds)

        #SeqIO.write(extracted_cds_list, filename +".CDS.fasta", "fasta")
    return

# creates a list of the files in this directory
raw_datadir_listing = os.listdir(os.getcwd())

# loops over the list of files
for files in raw_datadir_listing:
    if files.endswith('.gbk'):
        full_path = os.path.join(os.getcwd(), files)
        filename = os.path.splitext(files)[0]
        get_upstream_cds(full_path, filename)