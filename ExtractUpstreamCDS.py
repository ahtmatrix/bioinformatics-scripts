import sys
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
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
#turns warnings into errors so it can be caught


#grep '>' [filename] | wc -l

warnings.filterwarnings('error')


num_bp_upstreamcds = int(sys.argv[1])


def validate_cds(record, feature):
    feature_validity = None
    
    try:
        protein_in_file = str(feature.qualifiers.get('translation', 'no_translation')).strip('\'[]')
        #diff extracted CDS compare with FASTA nucleotdie on NCBI
        
        #is the problem with
        cds_to_protein = str(feature.extract(record).seq.translate(to_stop = True))
    
       
        #print "protein_check_fail: |" + record.id +"|"+ str(feature.qualifiers.get('transl_except')).strip('\'[]')  +"| " +str(feature.qualifiers.get('note')).strip('\'[]')  + " |"+str(feature.qualifiers.get('protein_id')).strip('\'[]') + "| " + protein_in_file + " |" + cds_to_protein
       
        temp_fix_protein = list(cds_to_protein)
        temp_fix_protein[0] = protein_in_file[0]
        
        fixed_cds_to_protein = "".join(temp_fix_protein)
        
        if fixed_cds_to_protein != protein_in_file:
            
            feature_validity = False
        else:
            feature_validity = True
            
    except BiopythonWarning:
        print "Biopythonwarning:" + record.id +" -->  "+ str(feature.qualifiers.get('protein_id')).strip('\'[]')# + " " + protein_in_file
        
    return feature_validity


def get_upstream_cds(fullpath, filename):

    extracted_cds_list = []
    
    # reads in a gbk and creates a Seqrecordord object
    for record in SeqIO.parse(fullpath, "genbank"):
        
        if record.features:
            for feature in record.features:
                if feature.type == "CDS":
                    if validate_cds(record, feature) == True:
                        #get the CDS nucleotide locations
                        cds_start = feature.location.start.position
                        cds_end = feature.location.end.position
                        
                    
                        #create a SeqFeature object containing the location of where to extract
                        # need to test if its taking + or - 1 off the location
                        # genbank starts with 1
                        upstream_and_cds            = SeqFeature(FeatureLocation(cds_start-num_bp_upstreamcds, cds_end))
                        extracted_upstream_and_cds  = upstream_and_cds.extract(record)
                        
                        
                        #only used for length culling
                        upstream_only               = SeqFeature(FeatureLocation(cds_start-num_bp_upstreamcds, cds_start))
                        extracted_upstream_only     = upstream_only.extract(record)
                        
                        
                        cds_protein_id = str(feature.qualifiers.get('protein_id')).strip('\'[]')
                        
                        
                        if len(extracted_upstream_only.seq) == num_bp_upstreamcds:
                            
                            annotated_record = SeqRecord(extracted_upstream_and_cds.seq, extracted_upstream_and_cds.id, description = "|" + cds_protein_id +"|")
                            
                            extracted_cds_list.append(annotated_record)

    #extraction is using the GENBANK protein for all
    SeqIO.write(extracted_cds_list, filename +"_"+str(num_bp_upstreamcds)+"upstream_"+".CDS.fasta", "fasta")
    return

# creates a list of the files in this directory
raw_datadir_listing = os.listdir(os.getcwd())

# loops over the list of files
for files in raw_datadir_listing:
    if files.endswith('.gbk'):
        full_path = os.path.join(os.getcwd(), files)
        filename = os.path.splitext(files)[0]
        get_upstream_cds(full_path, filename)