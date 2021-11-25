import pandas as pd
import os
import pysam
import sys
import configparser
import vcf

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

def vcf_filter_nucli_getero(my_id):
    print(my_id+' vcf_filter_nucli_getero started')
    read_shape = 0
    write_shape_nors = 0
    write_shape_rs = 0
    vcf_reader = vcf.Reader(open(os.path.join(processed_data, my_id + '/' + my_id + '.vcf'), 'r'))
    vcf_writer_nors = vcf.Writer(open(os.path.join(processed_data, my_id + '/' + my_id + '_nucli_getero_filtrated.vcf'), 'w'), vcf_reader)
    vcf_writer_rs = vcf.Writer(open(os.path.join(processed_data, my_id + '/' + my_id + '_rs_nucli_getero_filtrated.vcf'), 'w'), vcf_reader)
    for record in vcf_reader:
        read_shape += 1
        if((record.genotype('20')['GT']=='0/1')
        and (str(record.ALT[0]) in nucliotides) and
        (record.REF in nucliotides)):
            vcf_writer_nors.write_record(record)
            write_shape_nors += 1
        if ((record.genotype('20')['GT'] == '0/1')
            and (str(record.ALT[0]) in nucliotides) and
            (record.REF in nucliotides) and (record.ID != None)):
            vcf_writer_rs.write_record(record)
            write_shape_rs += 1
    print(my_id+' vcf_filter_nucli_getero finished')
    vcf_writer_rs.close()
    vcf_writer_nors.close()
    dict = {'read_shape': read_shape,'write_shape_nors':write_shape_nors,'write_shape_rs':write_shape_rs}
    return (dict)

def vcf_to_tsv(my_id):
    print(my_id + ' vcf to tsv started')
    vcf_reader = vcf.Reader(open(os.path.join(processed_data, my_id + '/' + my_id + '_rs_nucli_getero_filtrated.vcf'), 'r'))
    vcf_writer = open(os.path.join(processed_data, my_id + '/' + my_id + '_rs_nucli_getero_filtrated.tsv'),"w")
    for record in vcf_reader:
        vcf_writer.write(record.CHROM)
        vcf_writer.write('\t')
        vcf_writer.write(str(record.POS))
        vcf_writer.write('\t')
        vcf_writer.write(record.ID)
        vcf_writer.write('\t')
        vcf_writer.write(record.REF)
        vcf_writer.write('\t')
        vcf_writer.write(str(record.ALT[0]))
        vcf_writer.write('\t')
        vcf_writer.write(str(record.genotype('20')['AD'][0]))
        vcf_writer.write('\t')
        vcf_writer.write(str(record.genotype('20')['AD'][1]))
        vcf_writer.write('\n')
    vcf_writer.close()
    print(my_id + ' vcf to csv finished')

def filter_bad(my_id,threshold_bad):
    print(my_id + ' filterbad started')
    rs_tsv_table = pd.read_csv(os.path.join(processed_data, my_id + '/' + my_id + '_rs_nucli_getero_filtrated.tsv'),sep='\t',header=None)
    rs_tsv_table.columns = header_list
    #chip_c = cool_df[cool_df['BADgroup']=='chipseq']
    rs_filtered_table = rs_tsv_table[(rs_tsv_table['REF_COUNTS']>=threshold_bad) & (rs_tsv_table['ALT_COUNTS']>=threshold_bad)]
    rs_filtered_table.to_csv(os.path.join(processed_data, my_id + '/' + my_id + '_bad_rs_nucli_getero_filtrated.tsv'),index=False,sep='\t')

    '''
    nors_tsv_table = pd.read_csv(os.path.join(processed_data, my_id + '/' + my_id + '_nucli_getero_filtrated.tsv'),
                               sep='\t', header=None)
    nors_tsv_table.columns = header_list
    # chip_c = cool_df[cool_df['BADgroup']=='chipseq']
    nors_filtered_table = nors_tsv_table[(nors_tsv_table['REF_COUNTS'] >= threshold_bad) & (nors_tsv_table['ALT_COUNTS'] >= threshold_bad)]
    nors_filtered_table.to_csv(os.path.join(processed_data, my_id + '/' + my_id + '_bad_nucli_getero_filtrated.tsv'),
                             index=False, sep='\t')
    
    dict = {'write_shape_nors':nors_filtered_table.shape[0],'write_shape_rs':rs_filtered_table.shape[0]}
    '''
    print(my_id + ' filterbad finished')
    return rs_filtered_table.shape[0]




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

threshold_bad = int(sys.argv[1])
header_list = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'REF_COUNTS', 'ALT_COUNTS']
nucliotides = {'A', 'T', 'G', 'C'}

table = pd.DataFrame(columns=['ID', 'BAM reads count', 'Total num. of SNVs', 'Total num. of rsSNPs',
                              'rsSNPs passing ' + str(threshold_bad) + ' coverage'])
ref_alt_paths = []
for my_id in processing_list:
    # filtrate and get stats
    num = stats_read_num(os.path.join(indir, my_id + '.bam'))
    dict_vcf_filter = vcf_filter_nucli_getero(my_id)
    snv_num = dict_vcf_filter['read_shape']
    rssnp_num = dict_vcf_filter['write_shape_rs']

    vcf_to_tsv(my_id)
    snp_bad_rs_num = filter_bad(my_id,threshold_bad)

    to_append = [my_id, num, snv_num, rssnp_num, snp_bad_rs_num]
    a_series = pd.Series(to_append, index=table.columns)
    table = table.append(a_series, ignore_index=True)
    print(a_series)

    # get alt ref table
    rs_filtered_table = pd.read_csv(os.path.join(processed_data, my_id + '/' + my_id + '_bad_rs_nucli_getero_filtrated.tsv'), sep='\t')
    #rs_filtered_table.columns = header_list
    ref_alt = rs_filtered_table[['REF_COUNTS','ALT_COUNTS']]
    ref_alt.columns = ['ref','alt']
    ref_alt.to_csv(os.path.join(processed_data, my_id + '/' + my_id + 'ref_alt.tsv'))
    ref_alt_paths.append(os.path.join(processed_data, my_id + '/' + my_id + 'ref_alt.tsv'))




print(table)
metadata = pd.read_csv(met, sep='\t')
df = pd.merge(table, metadata, on="ID")
df.to_csv(maindir + 'logs/vcf_fast_filtration_stats.tsv', sep='\t',index=False)


# make ref alt count for all
output_file = os.path.join(processed_data, 'ref_alt_count.tsv')
print(ref_alt_paths)
combined_csv = pd.concat([pd.read_csv(f, sep='\t') for f in ref_alt_paths])
grp = (combined_csv.groupby(['ref', 'alt']).size()
       .reset_index(name='count'))
grp.to_csv(output_file, index=False, sep='\t')