import sys
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import pairwise2
import warnings
from Bio import BiopythonWarning
import re

# Usage
# python ExtractUpstreamCDS.py [number of bases upstream to cut]

# to combine multiple .gbk files into 1 gbk
# navigate to directory containing gbk files
# cat *.gbk > filename.gbk

# python -m pdb -Werror myprogram.py to run and stop pdb at the warning
# python -Werror -m pdb ExtractUpstreamCDS.py 30
# turns warnings into errors so it can be caught
warnings.filterwarnings('error')

# positive number upstream of start codon
num_bp_upstream_start = int(sys.argv[1])

# positive number downstream of start codon
num_bp_downstream_start = int(sys.argv[2])


def get_TIS(fullpath, filename):

    list_of_numbers_in_filename = re.findall('\d+', filename)
    num_bp_upstreamcds = int(list_of_numbers_in_filename[0])

    extracted_TIS_list = []

    TIS_coordinates = SeqFeature(FeatureLocation(
        num_bp_upstreamcds - num_bp_upstream_start, num_bp_upstreamcds + num_bp_downstream_start))
    # reads in a gbk and creates a Seqrecordord object
    for record in SeqIO.parse(fullpath, "fasta"):
        TIS_only_record = TIS_coordinates.extract(record)

        #annotated_TIS_only_record = SeqRecord(TIS_only_record.seq, TIS_only_record.id, description = "|" + cds_protein_id +"|")

        extracted_TIS_list.append(TIS_only_record)

    SeqIO.write(extracted_TIS_list, "extracted_TIS_" +
                filename + ".TIS.fasta", "fasta")
    return

# creates a list of the files in this directory
raw_datadir_listing = os.listdir(os.getcwd())

# loops over the list of files
for files in raw_datadir_listing:
    if files.endswith('upstream_.CDS.fasta'):
        full_path = os.path.join(os.getcwd(), files)
        filename = os.path.splitext(files)[0]
        get_TIS(full_path, filename)
