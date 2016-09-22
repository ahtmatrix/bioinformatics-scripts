from Bio import SeqIO
from Bio import SeqFeature
from Bio.SeqFeature import SeqFeature, FeatureLocation
import sys
import os

#Usage
#python [FULL PATH OF ExtractedUstreamCDS.py] ["fasta" or "genbank"] [FULL PATH OF GBK FILES TO EXTRACT]

#either fasta or genbank
filetype = sys.argv[1]
rawdata_location = sys.argv[2]

def getUpstreamCDS(filename):

    CDSRecordList = []
# reads in a gbk and creates a Seqrecordord object
    for record in SeqIO.parse(filename, "genbank"):
        
        
        
        if record.features:
            for feature in record.features:
                if feature.type == "CDS":
                    cds_location=feature.location
                    start=cds_location.start.position
                    upstream_cds = SeqFeature(FeatureLocation(0, start))
                    
                    CDSRecordList.append(upstream_cds.extract(record))
                    
                    #translate to double check  ?? WHy?
                    #
    if filetype == "fasta":
        SeqIO.write(CDSRecordList, filename+"CDS.fasta", filetype)
    elif filetype == "genbank":
        SeqIO.write(CDSRecordList, filename+"CDS.gbk", filetype)
    else:
        print("use either 'fasta' or 'genbank' in first argument")


#creates a list of the files in this directory
directory = rawdata_location
listing = os.listdir(directory)

#loops over the list of files
for files in listing:
    
    if files.endswith('.gbk'):  
        full_name = os.path.join(directory,files)
        
        getUpstreamCDS(full_name)
    else:
        print("this program requires .gbk files to run")