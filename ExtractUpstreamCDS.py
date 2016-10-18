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

# python -m pdb -Werror myprogram.py to run and stop pdb at the warning
# python -Werror -m pdb ExtractUpstreamCDS.py 30
#warnings.filterwarnings('error')
num_bp_upstreamcds = int(sys.argv[1])

def validate_cds(record, feature):
    
    try:
        
        if str(feature.qualifiers.get('protein_id')).strip('\'[]') == "YP_680427.1":
        
        # ACG codon for AAV YP_680427.1 will be counted
        # all sequences that differ by just 1 start codon will be counted
            print 
        
        
        
        #for splicing a seq
        #print feature.extract(record).seq[3::3]
            
        # if str(feature.qualifiers.get('protein_id')).strip('\'[]') == "YP_680427.1":
        #     for x in range(21, 27):
        #         print "table = " + str(x)
        #         print "\n"
        #         print str(feature.extract(record).seq.translate(table = x,to_stop = True))
        #         print "\n"
        
        protein_in_file = str(feature.qualifiers.get('translation', 'no_translation')).strip('\'[]')
        #diff extracted CDS compare with FASTA nucleotdie on NCBI
        
        #is the problem with 
        cds_to_protein = str(feature.extract(record).seq.translate(to_stop = True))
        
    #    if protein_in_file != cds_to_protein:
            # print feature.location_operator
    #        print "protein check fail: " + record.id +" "+ str(feature.qualifiers.get('protein_id')).strip('\'[]')
        
    except BiopythonWarning:
        print record.id +" -->  "+ str(feature.qualifiers.get('protein_id')).strip('\'[]')
        
    
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
                    
                    # need to test if its taking + or - 1 off the location
                    # genbank starts with 1
                    upstream_cds = SeqFeature(FeatureLocation(start-num_bp_upstreamcds, start))
                    
                    extracted_upstream_cds = upstream_cds.extract(record)
                    
                    if len(extracted_upstream_cds.seq) != num_bp_upstreamcds:
                        continue
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
