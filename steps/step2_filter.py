import subprocess32 as subprocess
import os
import configparser
import sys as sys
import pathlib
import pandas as pd
from pybedtools import BedTool
import shlex
def read_cfg(configpath):
    config = configparser.ConfigParser()
    config.read(configpath)
    maindir = config["Directories"]["maindir"]
    print(maindir)
    filt_bed=os.path.join(maindir,config["Directories"]["vcf_filtered"])
    source_vcf=os.path.join(config["Directories"]["final_data_out"])
    filt_bed_rs=os.path.join(maindir,config["Directories"]["rssnps"])
    mainlogs = os.path.join(maindir, config["Directories"]["mainlogs"])
    all_log = os.path.join( mainlogs, 'babachilogs')
    met = os.path.join(maindir, config["Files"]["metadata"])
    statfile = os.path.join(mainlogs, 'stats.tsv')
    metadata = pd.read_csv(met, sep='\t')
    with open(all_log, "w") as log:
        log.write('STARTING! all'+ '\n')

    metadata['starting_snps'] = 0
    metadata['filtrated_snps'] = 0
    metadata['rssnps'] = 0
    st_index = metadata.get_loc("starting_snps")
    filt_index = metadata.get_loc("filtrated_snps")
    rs_index = metadata.get_loc("rssnps")
    #metadata['bamsize MB'] = 0

def filter(configpath,trs,jobs): #filter and get stats
    config = configparser.ConfigParser()
    config.read(configpath)
    maindir = config["Directories"]["maindir"]
    print(maindir)
    filt_bed=os.path.join(maindir,config["Directories"]["vcf_filtered"])
    source_vcf=os.path.join(config["Directories"]["final_data_out"])
    filt_bed_rs=os.path.join(maindir,config["Directories"]["rssnps"])
    mainlogs = os.path.join(maindir, config["Directories"]["mainlogs"])
    all_log = os.path.join( mainlogs, 'babachilogs')
    met = os.path.join(maindir, config["Files"]["metadata"])
    statfile = os.path.join(mainlogs, 'stats.tsv')
    metadata = pd.read_csv(met, sep='\t')
    with open(all_log, "w") as log:
        log.write('STARTING! all'+ '\n')

    metadata['starting_snps'] = 0
    metadata['filtrated_snps'] = 0
    metadata['rssnps'] = 0
    st_index = metadata.columns.get_loc("starting_snps")
    filt_index = metadata.columns.get_loc("filtrated_snps")
    rs_index = metadata.columns.get_loc("rssnps")
    #metadata['bamsize MB'] = 0

    # parallel -j 8 babachi filter -O /media/ElissarDisk/ADASTRA/neuro/processed_neuro/vcf_filtered
    os.chdir(source_vcf)
    process = subprocess.Popen(
        [shlex.split('find -type f -name "*.vcf" | parallel -j '+ jobs+' babachi filter -a '+ trs, ' -O ', filt_bed)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True)
    stderr, stdout = process.communicate()
    with open(all_log, "a") as log:
        log.write(stdout)
        log.write(stderr)
    myids = list(metadata['ID'])
    for i in myids:
        # get starting vcf
        startvcf = os.path.join(source_vcf, i + '.vcf')
        file_rs = open(startvcf, "r")
        line_rs = file_rs.readline()
        n = 0
        while line_rs.startswith("##"):
            n += 1
            line_rs = file_rs.readline()
        file_rs.close()

        vcf = pd.read_csv(startvcf,
                          sep='\t', skiprows=n)

        # get bed file
        header_list = ['#chr', 'start', 'end', 'ID', 'ref', 'alt', 'ref_counts', 'alt_counts', 'sample_id']
        strfb = os.path.join(filt_bed, i + '.snps.bed')

        bedf = pd.read_csv(strfb, sep='\t', names=header_list)
        bedf['sample_id'] = i
        print(i)
        bedf.to_csv(strfb, sep='\t', index=False, header=None)

        # make bed file with rs
        bedrs = bedf[bedf['ID'].str.contains("rs")]
        bedrs.to_csv(os.path.join(filt_bed_rs, i + '.snps.bed'))

        # get stats
        a = list(metadata.index[metadata['ID'] == i])
        loc = a[0]

        metadata.iloc[loc, st_index] = vcf.shape[0]
        metadata.iloc[loc, filt_index] = bedf.shape[0]
        metadata.iloc[loc, rs_index] = bedrs.shape[0]

    metadata.to_csv(statfile, index=False, sep='\t')

    '''
    try:
        #parallel -j 8 babachi filter -O /media/ElissarDisk/ADASTRA/neuro/processed_neuro/vcf_filtered
        os.chdir(source_vcf)
        process = subprocess.Popen(['find -type f -name "*.vcf"','|','parallel','-j', jobs,'babachi filter','-a',trs ,'-O' ,filt_bed],
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE,
                                        universal_newlines=True)
        stderr, stdout = process.communicate()
        with open(all_log, "a") as log:
            log.write(stdout)
            log.write(stderr)
        myids = list(metadata['ID'])
        for i in myids:
            #get starting vcf
            startvcf = os.path.join(source_vcf, i + '.vcf')
            file_rs = open(startvcf, "r")
            line_rs = file_rs.readline()
            n = 0
            while line_rs.startswith("##"):
                n += 1
                line_rs = file_rs.readline()
            file_rs.close()

            vcf = pd.read_csv(startvcf,
                              sep='\t', skiprows=n)

            #get bed file
            header_list = ['#chr', 'start', 'end', 'ID', 'ref', 'alt', 'ref_counts', 'alt_counts', 'sample_id']
            strfb = os.path.join(filt_bed, i + '.snps.bed')

            bedf = pd.read_csv(strfb, sep='\t', names=header_list)
            bedf['sample_id'] = i
            print(i)
            bedf.to_csv(strfb, sep='\t', index=False, header=None)

            #make bed file with rs
            bedrs = bedf[bedf['ID'].str.contains("rs")]
            bedrs.to_csv(os.path.join(filt_bed_rs,i+'.snps.bed'))

            #get stats
            a = list(metadata.index[metadata['ID'] == i])
            loc = a[0]

            metadata.iloc[loc, st_index] = vcf.shape[0]
            metadata.iloc[loc, filt_index] = bedf.shape[0]
            metadata.iloc[loc, rs_index] = bedrs.shape[0]

        metadata.to_csv(statfile, index=False, sep='\t')

        return(0)
    except Exception:
        with open(all_log, "a") as log:
            log.write(str(Exception) + '\n')
            log.write('Some exception emerged'+ '\n')
        sys.exit(10)
        return(2)
    '''

def pullsort(configpath):
    config = configparser.ConfigParser()
    config.read(configpath)
    maindir = config["Directories"]["maindir"]
    filt_bed_rs = os.path.join(maindir, config["Directories"]["rssnps"])
    tmp_path = os.path.join(maindir,config["Directories"],["babachi"])
    mainlogs = os.path.join(maindir, config["Directories"]["mainlogs"])
    all_log = os.path.join(mainlogs, 'babachilogs')
    met = os.path.join(maindir, config["Files"]["metadata"])
    metadata = pd.read_csv(met, sep='\t')
    lst = metadata.BADgroup.unique()
    lst[0]
    header_list = ['#CHROM', 'POS1', 'POS2', 'ID', 'REF', 'ALT', 'REF_COUNT', 'ALT_COUNT', 'SAMPLEID']
    #list of bads
    for i in lst:
        pulltmp = metadata[metadata['BADgroup'] == i]

        intersect = list(pulltmp['ID'])
        paths_rs = []
        for my_id in intersect:
            if (os.path.exists(os.path.join(filt_bed_rs, my_id + '.snps.bed'))):
                paths_rs.append(os.path.join(filt_bed_rs, my_id + '.snps.bed'))
        print(paths_rs)

        # читаем чипсеки, смотрим распределение
        # chip_list = []

        with open(os.path.join(tmp_path, 'bedshapes'), "w") as log:
            for my_id in list(metadata['ID']):
                if (os.path.exists(os.path.join(filt_bed_rs, my_id + '.snps.bed'))):
                    tempdf = pd.read_csv(os.path.join(filt_bed_rs, my_id + '.snps.bed'), sep='\t', names=header_list)
                    log.write(my_id + '\t' + str(int(tempdf.shape[0])) + '\n')

        pulledtmps_rs = pd.concat([pd.read_csv(f, sep='\t', names=header_list) for f in paths_rs])
        pulledtmps_rs.to_csv(os.path.join(tmp_path, 'ppulled' + i + '.tsv'), mode='w', header=False, index=False,
                             sep='\t')
        print(i + ' pulled')
        process = subprocess.run(['bedtools', 'sort', '-i',
                                  os.path.join(tmp_path, 'ppulled' + i + '.tsv')],
                                 stdout=open(os.path.join(tmp_path, 'pulled' + i + '.tsv'), "w"),
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True
                                 )
        stderr, stdout = process.communicate()
        with open(all_log, "a") as log:
            log.write(stdout)
            log.write(stderr)
        print(i + ' sorted')
        print(pulledtmps_rs.shape[0])

def pullsortv2(configpath):
    #vcf_list = vcf_data.values.tolist()
    #test = BedTool(vcf_list)
    config = configparser.ConfigParser()
    config.read(configpath)
    maindir = config["Directories"]["maindir"]
    filt_bed_rs = os.path.join(maindir, config["Directories"]["rssnps"])
    tmp_path = os.path.join(maindir, config["Directories"], ["babachi"])
    mainlogs = os.path.join(maindir, config["Directories"]["mainlogs"])
    all_log = os.path.join(mainlogs, 'babachilogs')
    met = os.path.join(maindir, config["Files"]["metadata"])
    metadata = pd.read_csv(met, sep='\t')
    lst = metadata.BADgroup.unique()
    lst[0]
    header_list = ['#CHROM', 'POS1', 'POS2', 'ID', 'REF', 'ALT', 'REF_COUNT', 'ALT_COUNT', 'SAMPLEID']
    # list of bads
    for i in lst:
        pulltmp = metadata[metadata['BADgroup'] == i]

        intersect = list(pulltmp['ID'])
        paths_rs = []
        for my_id in intersect:
            if (os.path.exists(os.path.join(filt_bed_rs, my_id + '.snps.bed'))):
                paths_rs.append(os.path.join(filt_bed_rs, my_id + '.snps.bed'))
        print(paths_rs)

        # читаем чипсеки, смотрим распределение
        # chip_list = []

        with open(os.path.join(tmp_path, 'bedshapes'), "w") as log:
            for my_id in list(metadata['ID']):
                if (os.path.exists(os.path.join(filt_bed_rs, my_id + '.snps.bed'))):
                    tempdf = pd.read_csv(os.path.join(filt_bed_rs, my_id + '.snps.bed'), sep='\t', names=header_list)
                    log.write(my_id + '\t' + str(int(tempdf.shape[0])) + '\n')

        pulledtmps_rs = pd.concat([pd.read_csv(f, sep='\t', names=header_list) for f in paths_rs])
        pulledtmps_rs.to_csv(os.path.join(tmp_path, 'ppulled' + i + '.tsv'), mode='w', header=False, index=False,
                             sep='\t')
        print(i + ' pulled')

        bed_list = pulledtmps_rs.values.tolist()
        test = BedTool(bed_list)
        sorteddf = test.sort().to_dataframe()
        sorteddf.to_csv(os.path.join(tmp_path,'pulled'+i+'.tsv'), header=False, index=False,
                             sep='\t')
        print(i + ' sorted')
        print(pulledtmps_rs.shape[0])



#def babachi(configpath):

'''
        with open(all_log, "a") as log:
            log.write("script has been performed successfully")
'''
