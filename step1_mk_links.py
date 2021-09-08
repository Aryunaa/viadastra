#/media/ElissarDisk/ADASTRA/envs/belok_py36/bin/python

import pandas as pd
import os
from os.path import dirname as up

updir=up(os.getcwd())
print(updir)

dest = updir+'/data/'
print('dest '+dest)
os.mkdir(dest)

source = updir+'/BAM/'
print('source '+source)
met = updir + '/parameters/metadata.tsv'
print('met '+met)
metadata = pd.read_csv(met,sep='\t')
df = metadata[['BAM','ID']]
id_bam = dict(zip(df.BAM, df.ID))
for bam in id_bam:
    print('idbam '+id_bam[bam])
    #os.mkdir(dest + id_bam[bam],mode=0o777, dir_fd=None)
    #os.symlink(source +bam, dest + id_bam[bam]+'/'+id_bam[bam])
    os.symlink(source + bam, dest + '/' + id_bam[bam]+'.bam')
print('symlinks created')