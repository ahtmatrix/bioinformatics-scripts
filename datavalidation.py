import sys
import os
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from difflib import SequenceMatcher

# Usage
# python [FULL PATH OF ExtractedUstreamCDS.py] ["fasta" or "genbank"]
# [FULL PATH OF GBK FILES TO EXTRACT]

# to combine multiple .gbk files into 1 gbk
#   navigate to directory containing gbk files
#   cat *.gbk > filename.gbk

num_bp_upstreamcds = int(sys.argv[1])



# fasta_string = open("sequence.CDS.fasta").read()
# result_handle = NCBIWWW.qblast("blastn", "nt", fasta_string)
# save_file = open("my_blast.xml", "w")
# save_file.write(result_handle.read())
# save_file.close()
# result_handle.close()


def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()

    
def get_upstream_cds(fullpath, filename):
    #rint "Working on " + files + "..."
    text_file = open("Output.csv", "a")
    text_file.write("record ID," + "protein ID,"+"gbk protein," + "translated protein," + "similarity %")
    text_file.write("\n")
    
    extracted_cds_list = []
    # reads in a gbk and creates a Seqrecordord object
    for record in SeqIO.parse(fullpath, "genbank"):
        if record.features:
            for feature in record.features:
                if feature.type == "CDS":
                    
                    
                    translated_protein = str(feature.qualifiers.get('translation', 'no_translation')).strip('\'[]')
                    cds_to_protein = str(feature.extract(record).seq.translate(to_stop = True))
                    
                    
                    if translated_protein != cds_to_protein:
                        
                        
                        similarity_ratio = str(similar(translated_protein, cds_to_protein))
                        
                        

                        # cds = feature.extract(record)
                    
                        # extracted_cds_list.append(cds)
                        print similarity_ratio
                        
                        

                        text_file.write(record.id +"," +str(feature.qualifiers.get('protein_id')).strip('\'[]')+","+ translated_protein + "," + cds_to_protein + "," + similarity_ratio)
                        text_file.write("\n")
                
    text_file.close()

    return

# creates a list of the files in this directory
raw_datadir_listing = os.listdir(os.getcwd())
# loops over the list of files
for files in raw_datadir_listing:
    if files.endswith('.gbk'):
        full_path = os.path.join(os.getcwd(), files)
        filename = os.path.splitext(files)[0]
        get_upstream_cds(full_path, filename)


# finished_listing = os.listdir(directory)
# for files in finished_listing:
#     if files.endswith('.CDS.gbk'):
#         os.rename(files, os.getcwd()+"/extracted_seq/"+files)
