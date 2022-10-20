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

def vcf_filter_bad(my_id,treshold,qthreshold,rs):
    read_shape = 0
    write_shape = 0
    print(my_id + ' vcf_filter_bad started')
    if (rs!='rs'):
        vcf_reader = vcf.Reader(open(os.path.join(processed_data, my_id + '/' + my_id + '.vcf'), 'r'))
        vcf_writer = vcf.Writer(open(os.path.join(processed_data, my_id + '/' + my_id + 'bad_qual_nucli_getero_filtrated.vcf'), 'w'), vcf_reader)
        for record in vcf_reader:
            read_shape+=1
            if((record.genotype('20')['GT']=='0/1') and record.QUAL>= qthreshold and
                record.genotype('20')['AD'][0]>=treshold and record.genotype('20')['AD'][1]>=treshold and
               (record.REF in nucliotides) and (str(record.ALT[0]) in nucliotides)):
                vcf_writer.write_record(record)
                write_shape+=1
    else:
        vcf_reader = vcf.Reader(open(os.path.join(processed_data, my_id + '/' + my_id + '.vcf'), 'r'))
        vcf_writer = vcf.Writer(open(os.path.join(processed_data, my_id + '/' + my_id + 'bad_qual_rs_nucli_getero_filtrated.vcf'), 'w'), vcf_reader)
        for record in vcf_reader:
            read_shape += 1
            if ((record.genotype('20')['GT'] == '0/1') and record.QUAL>= qthreshold and
                 record.genotype('20')['AD'][0]>=treshold and record.genotype('20')['AD'][1]>=treshold and
                (record.REF in nucliotides) and (str(record.ALT[0]) in nucliotides) and (record.ID != None)):
                vcf_writer.write_record(record)
                write_shape += 1
    print(my_id + ' vcf_filter_bad finished')
    dict = {'write_shape':write_shape, 'read_shape':read_shape}
    return(dict)


def vcf_filter_asb(my_id,treshold,rs):
    print(my_id + ' vcf_filter_asb started')
    read_shape = 0
    write_shape = 0
    if (rs != 'rs'):
        vcf_reader = vcf.Reader(open(os.path.join(processed_data, my_id + '/' + my_id + '.vcf'), 'r'))
        vcf_writer = vcf.Writer(
            open(os.path.join(processed_data, my_id + '/' + my_id + 'asb_nucli_getero_filtrated.vcf'), 'w'), vcf_reader)
        for record in vcf_reader:
            read_shape += 1
            if ((record.genotype('20')['GT'] == '0/1')
                    and (str(record.ALT[0]) in nucliotides) and (record.REF in nucliotides)
                and record.genotype('20')['DP']>=treshold):
                vcf_writer.write_record(record)
                write_shape += 1
    else:
        vcf_reader = vcf.Reader(open(os.path.join(processed_data, my_id + '/' + my_id + '.vcf'), 'r'))
        vcf_writer = vcf.Writer(
            open(os.path.join(processed_data, my_id + '/' + my_id + 'asb_rs_nucli_getero_filtrated.vcf'), 'w'),
            vcf_reader)
        for record in vcf_reader:
            read_shape += 1
            if ((record.genotype('20')['GT'] == '0/1')
                    and (str(record.ALT[0]) in nucliotides) and (record.REF in nucliotides)
                    and record.genotype('20')['DP']>=treshold and (record.ID != None)):
                vcf_writer.write_record(record)
                write_shape += 1
    print(my_id + ' vcf_filter_asb finished')
    dict = {'write_shape': write_shape, 'read_shape': read_shape}
    return (dict)

#main filter
def vcf_filter_nucli_getero(my_id,rs):
    print(my_id+' vcf_filter_nucli_getero started')
    read_shape = 0
    write_shape = 0
    if (rs!='rs'):
        vcf_reader = vcf.Reader(open(os.path.join(processed_data, my_id + '/' + my_id + '.vcf'), 'r'))
        vcf_writer = vcf.Writer(open(os.path.join(processed_data, my_id + '/' + my_id + '_nucli_getero_filtrated.vcf'), 'w'), vcf_reader)
        for record in vcf_reader:
            read_shape += 1
            if((record.genotype('20')['GT']=='0/1')
               and (str(record.ALT[0]) in nucliotides) and
               (record.REF in nucliotides)):
                vcf_writer.write_record(record)
                write_shape += 1
    else:
        vcf_reader = vcf.Reader(open(os.path.join(processed_data, my_id + '/' + my_id + '.vcf'), 'r'))
        vcf_writer = vcf.Writer(open(os.path.join(processed_data, my_id + '/' + my_id + '_rs_nucli_getero_filtrated.vcf'), 'w'), vcf_reader)
        for record in vcf_reader:
            read_shape += 1
            if ((record.genotype('20')['GT'] == '0/1')
                    and (str(record.ALT[0]) in nucliotides) and
                    (record.REF in nucliotides) and (record.ID != None)):
                vcf_writer.write_record(record)
                write_shape += 1
    print(my_id+' vcf_filter_nucli_getero finished')
    dict = {'write_shape': write_shape, 'read_shape': read_shape}
    return (dict)


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


path = sys.argv[4]
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

treshold_bad = int(sys.argv[1])
treshold_asb = int(sys.argv[2])
treshold_qual = int(sys.argv[3])
nucliotides = {'A', 'T', 'G', 'C'}
table = pd.DataFrame(columns=['ID','BAM reads count', 'Total num. of SNVs','Total num. of rsSNPs',
                              'SNVs passing ' + str(treshold_bad) +' coverage and QUAL '+ str(treshold_qual)+ ' filtering',
                              'rsSNPs passing' + str(treshold_bad) +' coverage and QUAL '+ str(treshold_qual)+ ' filtering',
                              'Deep-covered SNVs' + ' (>='+ str(treshold_asb) +' reads)',
                              'Deep-covered rsSNPs (>=' + str(treshold_asb)+ ' reads)'])
bad_paths_rs = []
ref_alt_paths = []
for my_id in processing_list:
    #filtrate and get stats
    num = stats_read_num(os.path.join(indir,my_id + '.bam'))
    #num_dup = stats_read_num(os.path.join(processed_data, my_id + '/' + my_id + '_ready.bam'))
    dict_rs=vcf_filter_nucli_getero(my_id,'rs')
    rssnp_num=dict_rs['write_shape']

    
    dict_bad_rs = vcf_filter_bad(my_id,treshold_bad,treshold_qual,'rs')
    snp_num = dict_bad_rs['read_shape']
    snp_bad_rs_num = dict_bad_rs['write_shape']
    dict_bad_nors = vcf_filter_bad(my_id,treshold_bad,treshold_qual,'nors')
    snp_badfiltr_num = dict_bad_nors['write_shape']

    dict_asb_rs = vcf_filter_asb(my_id,treshold_asb,'rs')
    dict_asb_nors = vcf_filter_asb(my_id, treshold_asb, 'nors')
    snp_asb_nors = dict_asb_nors['write_shape']
    snp_asb_rs = dict_asb_rs['write_shape']
    to_append = [my_id, num,snp_num,rssnp_num,snp_badfiltr_num, snp_bad_rs_num,snp_asb_nors,snp_asb_rs]
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
metadata = pd.read_csv(met,sep='\t')
df = pd.merge(table, metadata, on="ID")
df.to_csv(maindir+'logs/vcf_filtration_stats.tsv',sep='\t')


#make rsID_count_bad_qual
for i in bad_paths_rs:
    print(i)
header_list = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
       '20']
combined_bad_rs = pd.concat([pd.read_csv(f,sep='\t',comment='#',names=header_list) for f in bad_paths_rs])
grp = (combined_bad_rs.groupby(['ID']).size()
       .reset_index(name='count'))
output_file = os.path.join(maindir,'logs/rsID_count_bad_qual.tsv')
grp.to_csv(output_file, index=False,sep='\t')

#make ref alt count for all
output_file = os.path.join(processed_data,'ref_alt_count.tsv')
print(ref_alt_paths)
combined_csv = pd.concat([pd.read_csv(f,sep='\t') for f in ref_alt_paths])
grp = (combined_csv.groupby(['ref', 'alt']).size()
       .reset_index(name='count'))
grp.to_csv(output_file, index=False,sep='\t')
