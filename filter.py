import vcf
import sys
import os
import configparser
import pandas as pd
nucliotides = {'A','T','G','C'}

def vcf_filter_nucli_getero(my_id,rs):
    if (rs!='rs'):
        vcf_reader = vcf.Reader(open(os.path.join(processed_data, my_id + '/' + my_id + '.vcf'), 'r'))
        vcf_writer = vcf.Writer(open(os.path.join(processed_data, my_id + '/' + my_id + '_nucli_getero_filtrated.vcf'), 'w'), vcf_reader)
        for record in vcf_reader:
            if((record.genotype('20')['GT']=='0/1')
               and (str(record.ALT[0]) in nucliotides) and
               (record.REF in nucliotides)):
                vcf_writer.write_record(record)
    else:
        vcf_reader = vcf.Reader(open(os.path.join(processed_data, my_id + '/' + my_id + '.vcf'), 'r'))
        vcf_writer = vcf.Writer(open(os.path.join(processed_data, my_id + '/' + my_id + '_rs_nucli_getero_filtrated.vcf'), 'w'), vcf_reader)
        for record in vcf_reader:
            if ((record.genotype('20')['GT'] == '0/1')
                    and (str(record.ALT[0]) in nucliotides) and
                    (record.REF in nucliotides) and (record.ID != None)):
                vcf_writer.write_record(record)

def ref_alt(my_id,rs):
    if (rs != 'rs'):
        vcf_filtrated = vcf.Reader(open(os.path.join(processed_data, my_id + '/' + my_id + '_nucli_getero_filtrated.vcf'), 'r'))
        # vcf_filtrated_wr = vcf.Writer(open('data/BAM00030_alt_ref', 'w'), vcf_reader)
        ref_alt = open(os.path.join(processed_data, my_id + '/' + my_id + '_alt_ref.tsv', "w"))
        ref_alt.write('ref' + '\t' + 'alt' + '\n')
        for record in vcf_filtrated:
            ref_alt.write(str(record.genotype('20')['AD'][0]))
            ref_alt.write('\t')
            ref_alt.write(str(record.genotype('20')['AD'][1]))
            ref_alt.write('\n')
        ref_alt.close()
    else:
        vcf_filtrated = vcf.Reader(open(os.path.join(processed_data, my_id + '/' + my_id + '_rs_nucli_getero_filtrated.vcf'), 'r'))
        # vcf_filtrated_wr = vcf.Writer(open('data/BAM00030_alt_ref', 'w'), vcf_reader)
        ref_alt = open(os.path.join(processed_data, my_id + '/' + my_id + '_rs_alt_ref.tsv'), "w")
        ref_alt.write('ref' + '\t' + 'alt' + '\n')
        for record in vcf_filtrated:
            ref_alt.write(str(record.genotype('20')['AD'][0]))
            ref_alt.write('\t')
            ref_alt.write(str(record.genotype('20')['AD'][1]))
            ref_alt.write('\n')
        ref_alt.close()


#path = sys.argv[4]
rs = sys.argv[1]
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

ref_alt_paths = []

for my_id in processing_list:

    print(my_id)

    vcf_filter_nucli_getero(my_id,rs)
    ref_alt(my_id,rs)

    if(rs!='rs'):
        ref_alt_paths.append(os.path.join(processed_data, my_id + '/' + my_id + '_alt_ref.tsv'))
    else:
        ref_alt_paths.append(os.path.join(processed_data, my_id + '/' + my_id + '_rs_alt_ref.tsv'))


output_file = os.path.join(processed_data,'ref_alt_count.tsv')

print(ref_alt_paths)
#combine all files in the list
combined_csv = pd.concat([pd.read_csv(f,sep='\t') for f in ref_alt_paths])
#export to csv
#combined_csv.to_csv(output_file, index=False)

grp = (combined_csv.groupby(['ref', 'alt']).size()
       .reset_index(name='count'))

grp.to_csv(output_file, index=False,sep='\t')
