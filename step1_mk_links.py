#/media/ElissarDisk/ADASTRA/envs/belok_py36/bin/python

import pandas as pd
import os
import sys
import configparser
import pysam


def read_head(BAM):
    samfile = pysam.view("-H", BAM)

    # for row in samfile:
    # print(row)

    # print(samfile)
    return (samfile)

def inters(samfile):
    proper_head = ['SN:chr1',
                   'SN:chr2',
                   'SN:chr3',
                   'SN:chr4',
                   'SN:chr5',
                   'SN:chr6',
                   'SN:chr7',
                   'SN:chr8',
                   'SN:chr9',
                   'SN:chr10',
                   'SN:chr11',
                   'SN:chr12',
                   'SN:chr13',
                   'SN:chr14',
                   'SN:chr15',
                   'SN:chr16',
                   'SN:chr17',
                   'SN:chr18',
                   'SN:chr19',
                   'SN:chr20',
                   'SN:chr21',
                   'SN:chr22',
                   'SN:chrX',
                   'SN:chrY'
                   ]
    x = pd.Series(samfile.split('\n'))
    x2 = x.str.split('\t')
    df = pd.DataFrame(x2.tolist())
    heads_chr = df.iloc[:, 1].tolist()
    inters = 0
    for i in heads_chr:
        for j in proper_head:
            if i == j:
                inters += 1
    print(inters)
    return (inters)


# reading config -----------------------------
path = sys.argv[1]
#path = 'CONFIG.cfg'
config = configparser.ConfigParser()
config.read(path)
maindir = config["Directories"]["maindir"]
print(maindir)
source = os.path.join(maindir,config["Directories"]["bam"])
print(source)
dest = os.path.join(maindir,config["Directories"]["data_in"])
print(dest)
met = os.path.join(maindir,config["Files"]["metadata"])
path_processing_list = os.path.join(maindir,config["Files"]["processing_list"])
path_exception_list = os.path.join(maindir,config["Files"]["exception_list"])


# creating exception list --------------------------
# creating processing list -------------------------
to_process_id = []
to_process_bam = []
exceptions_id = []
exceptions_bam = []
exceptions_cause = []
# reading metadata, filtrating by no RNA------------------------
print('met '+met)
metadata = pd.read_csv(met,sep='\t')
norna = metadata[metadata["Extra1"] != 'RNA_seq']
rna = metadata[metadata["Extra1"]=='RNA_seq']
for i in range(rna.shape[0]):
    exceptions_id.append(rna.iloc[i,1])
    exceptions_bam.append(rna.iloc[i,0])
    exceptions_cause.append('RNA_seq data')
# creating filtrating by headers + append to lists -------------
for i in range(norna.shape[0]):
    print(source)
    print(norna.iloc[i,0])
    pathfile = os.path.join(source,norna.iloc[i,0])
    print(pathfile)
    #pysam read -h
    bam = read_head(pathfile)
    #intersection by chr1,ch2 ...
    bam_inters = inters(bam)
    if(bam_inters>0):
        print('true')
        to_process_bam.append(norna.iloc[i,0])
        to_process_id.append(norna.iloc[i,1])
    else:
        exceptions_bam.append(norna.iloc[i,0])
        exceptions_id.append(norna.iloc[i,1])
        exceptions_cause.append('bams are not appropriate (not UCSC assembly)')

# to dicts, to dataframes
to_process = dict(zip(to_process_bam,to_process_id))
#exceptions = dict(zip(exceptions_bam,exceptions_id))
exceptions = pd.DataFrame(
    {'ID': exceptions_id,
     'cause': exceptions_cause,
     'bam': exceptions_bam
    })

# dataframes /  lists to files -------------------------------------
exceptions.to_csv(path_exception_list,sep='\t')

txt_processing_list = open(path_processing_list,"w")
for key in to_process:
    txt_processing_list.write(to_process[key]+"\n")
txt_processing_list.close()

# creating dict for symlinks ----------------------------------------------
id_bam = to_process

print(id_bam)
print('start_cycle')
for bam in id_bam:
    bai = bam.replace('.bam','.bai')
    print('idbam '+id_bam[bam])
    #os.mkdir(dest + id_bam[bam],mode=0o777, dir_fd=None)
    if(os.path.exists(dest + '/' + id_bam[bam]+'.bam')):
        os.remove(dest + '/' + id_bam[bam]+'.bam')
        os.symlink(source + bam, dest + '/' + id_bam[bam]+'.bam')
    else:
        os.symlink(source + bam, dest + '/' + id_bam[bam]+'.bam')
print('symlinks created')

