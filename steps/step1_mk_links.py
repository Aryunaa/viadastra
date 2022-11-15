#/media/ElissarDisk/ADASTRA/envs/belok_py36/bin/python

import pandas as pd
import os
import sys
import configparser
import pysam

#читает заголовки у бамки
def read_head(BAM):
    samfile = pysam.view("-H", BAM)

    # for row in samfile:
    # print(row)

    # print(samfile)
    return (samfile)

#сравнивает нужный хэдер и тот что выходит из пайсэма
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

def mklinks(parameter):# reading config -----------------------------
    if(os.path.exists(parameter)):
        try:
            path = parameter
            #path = 'CONFIG.cfg'
            '''
            create directories from config file
            '''


            config = configparser.ConfigParser()

            config.read(path)
            dirs = config["Directories"]
            maindir = config["Directories"]["maindir"]
            print(maindir)
            for dir in dirs:
                print(dir)
                dirp = os.path.join(maindir, config["Directories"][dir])

                if (os.path.exists(dirp) == False):
                    try:
                        os.makedirs(dirp, exist_ok=False)
                        print("Directory '%s' created successfully" % dirp)
                    except OSError as error:
                        print("Directory '%s' can not be created")


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
            # reading metadata------------------------
            print('met '+met)
            metadata = pd.read_csv(met,sep='\t')

            pathloc = metadata.columns.get_loc("path")
            idloc = metadata.columns.get_loc("ID")
            # creating filtrating by headers + append to lists -------------
            for i in range(metadata.shape[0]):
                print(source)
                print(metadata.iloc[i,pathloc])
                pathfile = os.path.join(source,metadata.iloc[i,pathloc])
                print(pathfile)
                #pysam read -h
                if(os.path.exists(pathfile)):
                    bam = read_head(pathfile)
                    #intersection by chr1,ch2 ...

                    bam_inters = inters(bam)
                    if(bam_inters>0):
                        print('true')
                        to_process_bam.append(metadata.iloc[i,idloc])
                        to_process_id.append(metadata.iloc[i,idloc])
                    else:
                        exceptions_bam.append(metadata.iloc[i,idloc])
                        exceptions_id.append(metadata.iloc[i,idloc])
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
            return (0)
        except Exception:
            print('failed to create ' + str(Exception))
            return (2)
    else:
        print('config '+ parameter+ ' does not exist!')
        return (2)




