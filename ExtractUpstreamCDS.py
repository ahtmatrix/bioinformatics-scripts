from Bio import SeqIO
from Bio import SeqFeature
from Bio.SeqFeature import SeqFeature, FeatureLocation
import sys
import os

upstream = int(sys.argv[1])


def getUpstreamCDS(filename):

    CDSRecordList = []
# reads in a gbk and creates a Seqrecordord object
    for record in SeqIO.parse(filename, "genbank"):
        
        
        
        if record.features:
            for feature in record.features:
                if feature.type == "CDS":
                    cds_location=feature.location
                    start=cds_location.start.position
                    upstream_cds = SeqFeature(FeatureLocation(start-upstream, start))
                    
                    CDSRecordList.append(upstream_cds.extract(record))
                    
                    
                    
    SeqIO.write(CDSRecordList, filename+"CDS.fasta", "fasta")


#creates a list of the files in this directory
directory = "/home/ubuntu/workspace/biopython-scripts/"
listing = os.listdir(directory)

#loops over the list of files
for files in listing:
    
    if files.endswith('.gbk'):  
        full_name = os.path.join(directory,files)
        getUpstreamCDS(full_name)