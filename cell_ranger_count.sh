#!/usr/bin/bash
date
#The directory to which the output must be saved 
cd /opt/localdata/msomada/gdm/multi_outs

#Add the directory path for the fastq files to the --fastqs argument, unique identifier for the samples and the path to the reference
/opt/localdata/softwares/cellranger-7.0.0/cellranger count --id=1000 \
   --fastqs=///  \
   --sample=-- \
   --transcriptome=//refdata-gex-GRCh38-2020-A/ \
   --localmem=80 --localcores=16

date