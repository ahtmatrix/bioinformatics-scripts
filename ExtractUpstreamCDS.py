import sys
import os
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

# Usage
# python [FULL PATH OF ExtractedUstreamCDS.py] ["fasta" or "genbank"]
# [FULL PATH OF GBK FILES TO EXTRACT]

# to combine multiple .gbk files into 1 gbk
#   navigate to directory containing gbk files
#   cat *.gbk > filename.gbk

# either fasta or genbank
output_filetype = sys.argv[1]
num_bp_upstreamcds = int(sys.argv[2])



    
def get_upstream_cds(fullpath, filename):
    #rint "Working on " + files + "..."

    count = 0
    extracted_cds_list = []
    # reads in a gbk and creates a Seqrecordord object
    for record in SeqIO.parse(fullpath, "genbank"):
        
        if record.features:
            for feature in record.features:
                if feature.type == "CDS":
                    
                    
                    translated_protein = str(feature.qualifiers.get('translation', 'no_translation')).strip('\'[]')
                    cds_to_protein = str(feature.extract(record).seq.translate(to_stop = True))
                    
                    if translated_protein != cds_to_protein:
                        
                        count = count + 1
                        print record.id
                        #print feature.location
                        print count
                        print translated_protein
                        print cds_to_protein
                    else:
                        continue
                    
                    cds_location = feature.location
                    start = cds_location.start.position
                    end = cds_location.end.position
                    
                    
                    
                    cds = SeqFeature(FeatureLocation(start, end))
                    
                    
                    
                    #print "\n"
                    
                    
                    
                    #translated sequence in ncbi is different from CDS extracted then translated
                    
                    
                    upstream_cds = SeqFeature(FeatureLocation(start-num_bp_upstreamcds, start))
                    
                    extracted_upstream_cds = upstream_cds.extract(record)
                    
                    if len(extracted_upstream_cds.seq) != num_bp_upstreamcds:
                        print "upstream cds length is too short = " + str(len(extracted_upstream_cds.seq))
                    else:
                        extracted_cds_list.append(extracted_upstream_cds)
                    # translate to double check  ?? WHy?
                    
    if output_filetype == "fasta":
        SeqIO.write(extracted_cds_list, filename +
                    ".CDS.fasta", output_filetype)
    elif output_filetype == "genbank":
        SeqIO.write(extracted_cds_list, filename + ".CDS.gbk", output_filetype)
    else:
        print "use either 'fasta' or 'genbank' in first argument"
    #print "Done with " + files + "..."
    
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
