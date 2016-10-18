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
#turns warnings into errors so it can be caught
warnings.filterwarnings('error')


num_bp_upstreamcds = int(sys.argv[1])


def validate_cds(record, feature):
    
    try:
        
        # ACG codon for AAV YP_680427.1 will be counted
        # all sequences that differ by just 1 start codon will be counted
        
        #check the start codon
        #   if start codon only differs by 1 let it through
        #
        
        #for splicing a seq
        #print feature.extract(record).seq[3::3]
        
        protein_in_file = str(feature.qualifiers.get('translation', 'no_translation')).strip('\'[]')
        #diff extracted CDS compare with FASTA nucleotdie on NCBI
        
        #is the problem with
        cds_to_protein = str(feature.extract(record).seq.translate(to_stop = True))
        
        
        #may have to reconstruct a Seq object for each extraction
       
        temp_fix_protein = list(cds_to_protein)
        temp_fix_protein[0] = protein_in_file[0]
        
        fixed_cds_to_protein = "".join(temp_fix_protein)
        
        if fixed_cds_to_protein != protein_in_file:
            print "protein_check_fail: |" + record.id +"|"+ str(feature.qualifiers.get('transl_except')).strip('\'[]')  +"| " +str(feature.qualifiers.get('note')).strip('\'[]')  + " |"+str(feature.qualifiers.get('protein_id')).strip('\'[]') + "| " + protein_in_file + " |" + cds_to_protein
            
            # if protein_in_file in cds_to_protein:
            #     print "protein_in_file is within cds_to_protein"
    
    
    
    except BiopythonWarning:
        print "Biopythonwarning:" + record.id +" -->  "+ str(feature.qualifiers.get('protein_id')).strip('\'[]')# + " " + protein_in_file
        
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