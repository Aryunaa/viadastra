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

def filter_bad(my_id,threshold_bad):




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
                              'SNVs passing ' + str(treshold_bad) + ' coverage',
                              'rsSNPs passing' + str(treshold_bad) + ' coverage' ])
bad_paths_rs = []
ref_alt_paths = []
for my_id in processing_list:
    # filtrate and get stats
    num = stats_read_num(os.path.join(indir, my_id + '.bam'))
    dict_vcf_filter = vcf_filter_nucli_getero(my_id)
    snp_num = dict_vcf_filter['read_shape']
    rssnp_num = dict_vcf_filter['write_shape_rs']


    snp_bad_rs_num = dict_bad['write_shape']

    snp_badfiltr_num = dict_bad['write_shape']

    to_append = [my_id, num, snp_num, rssnp_num, snp_badfiltr_num, snp_bad_rs_num]
    a_series = pd.Series(to_append, index=table.columns)
    table = table.append(a_series, ignore_index=True)
    print(a_series)

    bad_paths_rs.append(os.path.join(processed_data, my_id + '/' + my_id + 'bad_qual_rs_nucli_getero_filtrated.vcf'))

    # get alt ref table
    rs = 'rs'
    ref_alt(my_id, rs)
    if (rs != 'rs'):
        ref_alt_paths.append(os.path.join(processed_data, my_id + '/' + my_id + '_alt_ref.tsv'))
    else:
        ref_alt_paths.append(os.path.join(processed_data, my_id + '/' + my_id + '_rs_alt_ref.tsv'))

print(table)
metadata = pd.read_csv(met, sep='\t')
df = pd.merge(table, metadata, on="ID")
df.to_csv(maindir + 'logs/vcf_filtration_stats.tsv', sep='\t')

# make rsID_count_bad_qual
for i in bad_paths_rs:
    print(i)
header_list = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
               '20']
combined_bad_rs = pd.concat([pd.read_csv(f, sep='\t', comment='#', names=header_list) for f in bad_paths_rs])
grp = (combined_bad_rs.groupby(['ID']).size()
       .reset_index(name='count'))
output_file = os.path.join(maindir, 'logs/rsID_count_bad_qual.tsv')
grp.to_csv(output_file, index=False, sep='\t')

# make ref alt count for all
output_file = os.path.join(processed_data, 'ref_alt_count.tsv')
print(ref_alt_paths)
combined_csv = pd.concat([pd.read_csv(f, sep='\t') for f in ref_alt_paths])
grp = (combined_csv.groupby(['ref', 'alt']).size()
       .reset_index(name='count'))
grp.to_csv(output_file, index=False, sep='\t')