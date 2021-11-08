import vcf
import sys
import os
import configparser
import pandas as pd
nucliotides = {'A','T','G','C'}


def vcf_filter_bad(my_id,treshold):
    read_shape = 0
    write_shape = 0
    print(my_id + ' vcf_filter_bad started')

    vcf_reader = vcf.Reader(open(os.path.join(processed_data, my_id + '/' + my_id + '.vcf'), 'r'))
    vcf_writer = vcf.Writer(open(os.path.join(processed_data, my_id + '/' + my_id + '_bad_rs_nucli_getero_filtrated.vcf'), 'w'), vcf_reader)
    for record in vcf_reader:
        read_shape += 1
        if ((record.genotype('20')['GT'] == '0/1') and
            record.genotype('20')['AD'][0]>=treshold and record.genotype('20')['AD'][1]>=treshold and
            (record.REF in nucliotides) and (str(record.ALT[0]) in nucliotides) and (record.ID != None)):
            vcf_writer.write_record(record)
            write_shape += 1
    print(my_id + ' vcf_filter_bad finished')


def ref_alt(my_id):
    vcf_filtrated = vcf.Reader(open(os.path.join(processed_data, my_id + '/' + my_id + '_bad_rs_nucli_getero_filtrated.vcf'), 'r'))
    # vcf_filtrated_wr = vcf.Writer(open('data/BAM00030_alt_ref', 'w'), vcf_reader)
    ref_alt = open(os.path.join(processed_data, my_id + '/' + my_id + '_ad_alt_ref.tsv'), "w")
    ref_alt.write('alt_counts' + '\t' + 'ref_counts' + '\n')
    for record in vcf_filtrated:
        ref_alt.write(str(record.genotype('20')['AD'][1]))
        ref_alt.write('\t')
        ref_alt.write(str(record.genotype('20')['AD'][0]))
        ref_alt.write('\n')
    ref_alt.close()


#path = sys.argv[4]
threshold = int(sys.argv[1])
path = sys.argv[2]
config = configparser.ConfigParser()
config.read(path)
maindir = config["Directories"]["maindir"]
processed_data = os.path.join(maindir,config["Directories"]["data_out"])
met = os.path.join(maindir,config["Files"]["metadata"])
indir = os.path.join(maindir,config["Directories"]["data_in"])

processing_list_path = os.path.join(maindir,config["Files"]["processing_list"])
list_file = open(processing_list_path, "r")
listi = list_file.readlines()
processing_list = []
for i in listi:
    i = i.replace('\n','')
    processing_list.append(i)

ref_alt_paths_atac = []
ref_alt_paths_chip = []

metadata = pd.read_csv(met,sep='\t')

for my_id in processing_list:

    print(my_id)

    vcf_filter_bad(my_id,threshold)
    ref_alt(my_id)
    t = metadata[metadata['ID'] == my_id]
    badgr = t['Extra1']
    #badgr = badser.iloc[0]
    print(badgr)
    if (badgr.iloc[0]=='ChIP_seq' or badgr.iloc[0]=='DRIP_seq'): # ChIP_seq

        ref_alt_paths_chip.append(os.path.join(processed_data, my_id + '/' + my_id + '_ad_alt_ref.tsv'))
    elif (badgr.iloc[0]=='ATAC_seq'):
        ref_alt_paths_atac.append(os.path.join(processed_data, my_id + '/' + my_id + '_ad_alt_ref.tsv'))

output_file_chip = os.path.join(processed_data,'_chip_ref_alt_count.tsv')
output_file_atac = os.path.join(processed_data,'_atac_ref_alt_count.tsv')
print(ref_alt_paths_chip)
print(ref_alt_paths_atac)
#combine all files in the list
combined_csv_chip = pd.concat([pd.read_csv(f,sep='\t') for f in ref_alt_paths_chip])
combined_csv_atac = pd.concat([pd.read_csv(f,sep='\t') for f in ref_alt_paths_atac])
#export to csv
#combined_csv.to_csv(output_file, index=False)

grp_chip = (combined_csv_chip.groupby(['alt_counts', 'ref_counts']).size()
       .reset_index(name='counts'))
grp_atac = (combined_csv_atac.groupby(['alt_counts', 'ref_counts']).size()
       .reset_index(name='counts'))

grp_chip.to_csv(output_file_chip, index=False,sep='\t')
grp_atac.to_csv(output_file_atac, index=False,sep='\t')
