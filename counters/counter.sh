#!/bin/bash


echo "viral_30upstream_.CDS.TIS.fasta record count = "
grep ">" viral_30upstream_.CDS.TIS.fasta| wc -l

echo "viral_30upstream_.CDS.fasta record count = "
grep ">" viral_30upstream_.CDS.fasta| wc -l

echo "rna_30upstream_.CDS.TIS.fasta record count = "
grep ">" rna_30upstream_.CDS.TIS.fasta| wc -l 

echo "rna_30upstream_.CDS.fasta record count = "
grep ">" rna_30upstream_.CDS.fasta| wc -l 


echo "rna.gbk record count = "
grep "//" rna.gbk | wc -l

echo "viral.gbk record count = "
grep "//" viral.gbk | wc -l