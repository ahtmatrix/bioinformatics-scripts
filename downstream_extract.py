import sys
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import pairwise2
import warnings
from Bio import BiopythonWarning

# Usage
# python SeqExtract.py [number of bases upstream to cut]

# to combine multiple .gbk files into 1 gbk
# navigate to directory containing gbk files
# cat *.gbk > filename.gbk

# python -m pdb -Werror myprogram.py to run and stop pdb at the warning
# python -Werror -m pdb SeqExtract.py 30
# turns warnings into errors so it can be caught

# grep '>' [filename] | wc -l

warnings.filterwarnings('error')

#number of base pairs upstream to extract
num_bp_upstreamSTOP = int(sys.argv[1])

#number of base pairs downstream of the cds to extract
num_bp_downstreamSTOP = int(sys.argv[2])


def validate_cds(record, feature):
    feature_validity = None

    try:
        protein_in_file = str(feature.qualifiers.get('translation', 'no_translation')).strip('\'[]')
        # diff extracted CDS compare with FASTA nucleotide on NCBI

        # is the problem with?
        cds_to_protein = str(feature.extract(record).seq.translate(to_stop=True))

        # print "protein_check_fail: |" + record.id +"|"+
        # str(feature.qualifiers.get('transl_except')).strip('\'[]')  +"| "
        # +str(feature.qualifiers.get('note')).strip('\'[]')  + "
        # |"+str(feature.qualifiers.get('protein_id')).strip('\'[]') + "| " +
        # protein_in_file + " |" + cds_to_protein

        temp_fix_protein = list(cds_to_protein)
        temp_fix_protein[0] = protein_in_file[0]

        fixed_cds_to_protein = "".join(temp_fix_protein)

        if fixed_cds_to_protein != protein_in_file:
            feature_validity = False
        else:
            feature_validity = True

    except BiopythonWarning:
        print "Biopythonwarning:" + record.id + " -->  " + str(feature.qualifiers.get('protein_id')).strip('\'[]')

    return feature_validity


def extract_upstream_of_STOP_and_downstream_of_stop(fullpath, filename):

    extracted_cds_list = []

    # reads in a gbk and creates a SeqRecord object
    for record in SeqIO.parse(fullpath, "genbank"):
        if record.features:
            for feature in record.features:
                if feature.type == "CDS":
                    if validate_cds(record, feature) == True:
                        # get the CDS nucleotide locations
                        cds_start_location = feature.location.start.position
                        cds_end_location = feature.location.end.position
                        
                        #get the 3'UTR sequence coordinates and extract
                        ThreeUTR_location = SeqFeature(FeatureLocation(cds_end_location, cds_end_location + num_bp_downstreamSTOP))
                        extracted_3UTR = ThreeUTR_location.extract(record)
                        
                        if len(extracted_3UTR.seq) == num_bp_downstreamSTOP:
                            
                            #extract -num_bp_upstreamSTOP + STOP + num_bp_downstreamSTOP    
                            extract_location = SeqFeature(FeatureLocation(cds_end_location- num_bp_upstreamSTOP, cds_end_location + num_bp_downstreamSTOP))
                        
                            #need to check if complement
                            #if it is complemement, then reverse complement it
                            
                            if "+" in str(feature.location):
                                extracted_seq = extract_location.extract(record)
                                print "reverse complement disengaged" + str(feature.location)
                                
                            elif "-" in str(feature.location):
                                extracted_seq = extract_location.extract(record).reverse_complement()
                                print "reverse complement engaged   " + str(feature.location)
                                
                            cds_protein_id = str(feature.qualifiers.get('protein_id')).strip('\'[]')
                            annotated_record = SeqRecord(extracted_seq.seq, extracted_seq.name, description="|" + cds_protein_id + "|")
                            extracted_cds_list.append(annotated_record)

                            # create a SeqFeature object containing the location of where to extract
                            # need to test if its taking + or - 1 off the location
                            # genbank starts with 1
                #upstream_cds_downstream_location = SeqFeature(FeatureLocation(cds_start_location - num_bp_upstreamcds, cds_end_location + num_bp_downstreamcds))

    # extraction is using the GENBANK protein for all
    SeqIO.write(extracted_cds_list, filename + str(num_bp_upstreamSTOP) + "upstream_STOP_"+ str(num_bp_downstreamSTOP)+"downstream.fasta", "fasta")
    return

# creates a list of the files in this directory
raw_datadir_listing = os.listdir(os.getcwd())

# loops over the list of files
for files in raw_datadir_listing:
    if files.endswith('.gbk'):
        full_path = os.path.join(os.getcwd(), files)
        filename = os.path.splitext(files)[0]
        extract_upstream_of_STOP_and_downstream_of_stop(full_path, filename)