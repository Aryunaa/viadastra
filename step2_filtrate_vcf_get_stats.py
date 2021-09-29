import pandas as pd
import os
import pysam
import sys
import configparser

def stats_cov(BAM):
    samfile = pysam.coverage(BAM)

    # for row in samfile:
    # print(row)

    # print(samfile)
    return (samfile)


def stats_read_num(BAM):
    samfile = pysam.view("-c", BAM)

    # for row in samfile:
    # print(row)

    # print(samfile)
    num = int(samfile.replace('\n', ''))
    return (num)

def vcf_filter(my_id,treshold):
    print(my_id)
    file = open(os.path.join(processed_data, my_id + '/' + my_id + '.vcf'), "r")
    line = file.readline()
    n = 0
    while line.startswith("##"):
        n += 1
        line = file.readline()
    file.close()

    vcf_data = pd.read_csv(os.path.join(processed_data, my_id + '/' + my_id + '.vcf'), sep='\t', skiprows=n)
    vcf_filtrated = pd.DataFrame(columns=vcf_data.columns)
    #treshold = 5.5
    for i in range(vcf_data.shape[0]):
        temp = vcf_data.iloc[i, 9]
        print(temp)
        temp_list = temp.split(':')
        ad = temp_list[1].split(',')
        fad = int(ad[0])
        sad = int(ad[1])
        print(temp_list)
        if ((temp_list[0] == '1/0' or temp_list[0] == '0/1') and (fad > treshold and sad > treshold)
                and (vcf_data.iloc[i,3] in nucliotides) and (vcf_data.iloc[i,4] in nucliotides)) :
            vcf_filtrated = vcf_filtrated.append(vcf_data.iloc[i, :], ignore_index=False)
    vcf_filtrated.to_csv(os.path.join(processed_data, my_id + '/' + my_id + '_bad_filtrated.vcf'), sep='\t')
    return(vcf_filtrated.shape[0])

path = sys.argv[2]
#path = 'CONFIG.cfg'
config = configparser.ConfigParser()
config.read(path)
maindir = config["Directories"]["maindir"]
processed_data = os.path.join(maindir,config["Directories"]["data_out"])
met = os.path.join(maindir,config["Files"]["metadata"])
processing_list_path = os.path.join(maindir,config["Files"]["processing_list"])
indir = os.path.join(maindir,config["Directories"]["data_in"])


list_file = open(processing_list_path, "r")
listi = list_file.readlines()
processing_list = []
for i in listi:
    i = i.replace('\n','')
    processing_list.append(i)

treshold = int(sys.argv[1])
nucliotides = {'A', 'T', 'G', 'C'}
table = pd.DataFrame(columns=['ID','Dataset reads number', 'Duplicate reads number', 'SNP number','SNP passed BAD treshold'])
for my_id in processing_list:
    num = stats_read_num(os.path.join(indir,my_id + '.bam'))
    num_dup = stats_read_num(os.path.join(processed_data, my_id + '/' + my_id + '_ready.bam'))
    vcf_data = pd.read_csv(os.path.join(processed_data, my_id + '/' + my_id + '.vcf'), sep='\t', skiprows=48)
    snp_num = vcf_data.shape[0]
    snp_badfiltr_num = vcf_filter(my_id,treshold)
    to_append = [my_id, num, num_dup,snp_num,snp_badfiltr_num]
    a_series = pd.Series(to_append, index=table.columns)
    table = table.append(a_series, ignore_index=True)
    print(a_series)
    #print(table)

print(table)
table.to_csv(maindir+'logs/table.tsv',sep='\t')

