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
    dict = {'read_shape': read_shape,'write_shape_nors':write_shape_nors,'write_shape_rs':write_shape_rs}
    return (dict)

def vcf_to_tsv(my_id):
    print(my_id + ' vcf to tsv started')
