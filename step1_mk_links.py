#/media/ElissarDisk/ADASTRA/envs/belok_py36/bin/python

import pandas as pd
import os
import sys
import configparser

def readConfig(path, dir):
    config = configparser.ConfigParser()
    config.read(path)

    maindir = config.get("Directories", "maindir")
    bam = config.get("Directories", "bam")
    data_in = config.get("Directories", "data_in")
    metadata = config.get("Files", "metadata")

    if (dir == 'maindir'):
        return (maindir)
    elif (dir == 'metadata'):
        return (maindir+metadata)
    elif (dir == 'bam'):
        return (maindir+bam)
    elif (dir == 'data_in'):
        return (maindir+data_in)

path = sys.argv[1]
met = readConfig(path,'metadata')
maindir = readConfig(path,'maindir')
source = readConfig(path, 'bam')
dest = readConfig(path,'data_in')

print('met '+met)
metadata = pd.read_csv(met,sep='\t')
df = metadata[['BAM','ID']]
id_bam = dict(zip(df.BAM, df.ID))
for bam in id_bam:
    bai = bam.replace('.bam','.bai')
    print('idbam '+id_bam[bam])
    #os.mkdir(dest + id_bam[bam],mode=0o777, dir_fd=None)
    #os.symlink(source +bam, dest + id_bam[bam]+'/'+id_bam[bam])
    os.symlink(source + bam, dest + '/' + id_bam[bam]+'.bam')
    #os.symlink(source + bai, dest + '/' + id_bam[bam] + '.bai')
print('symlinks created')