import sys
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import pairwise2
import warnings
from Bio import BiopythonWarning
import string
import random
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from random import randint
from random import choice

from string import ascii_uppercase
import random
from itertools import islice
from Bio import Seq




def random_chars(size, chars="ATCG"):
    selection = iter(lambda: random.choice(chars), object())
    while True:
        yield ''.join(islice(selection, size))

random_gen_up = random_chars(9)
random_gen_down = random_chars(1000)


104360

2528

# random_SeqRecords=[]
# for i in range(0, 2528):
#     random_seq = str(next(random_gen_up) + "ATG" + next(random_gen_down))
#     record = SeqRecord(id = 'random_noise_%i' % (i+1), seq = Seq.Seq(random_seq))
#     random_SeqRecords.append(record)
    
    
pure_random_SeqRecords=[]
for i in range(0, 2528):
    random_seq = str(next(random_gen_down))
    record = SeqRecord(id = 'random_noise_%i' % (i+1), seq = Seq.Seq(random_seq))
    pure_random_SeqRecords.append(record)



SeqIO.write(pure_random_SeqRecords,"random_noise_viral.fasta","fasta")