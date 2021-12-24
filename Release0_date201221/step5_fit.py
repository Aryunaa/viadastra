import pandas as pd
import os
import pysam
import sys
import configparser
import pybedtools
from pybedtools import BedTool
import subprocess32 as subprocess
import shlex

#path = '/media/ElissarDisk/ADASTRA/parameters/CONFIG.cfg'
threshold = int(sys.argv[1])
path = sys.argv[2]
force = sys.argv[3]


config = configparser.ConfigParser()
config.read(path)
maindir = config["Directories"]["maindir"]
processed_data = os.path.join(maindir,config["Directories"]["data_out"])
met = os.path.join(maindir,config["Files"]["metadata"])
processing_list_path = os.path.join(maindir,config["Files"]["processing_list"])
indir = os.path.join(maindir,config["Directories"]["data_in"])
fit = os.path.join(maindir,config["Directories"]["fit"])
babachi = os.path.join(maindir,config["Directories"]["babachi"])
######################################################
def annotate_by_bad(i,threshold):
    if(os.path.exists(os.path.join(fit, i + '_BAD_annotated.tsv'))):
        pass
    else:
    #'pulled_'+i+'_tobabachi.tsv', 'pulled_'+i+'_tobabachi.bed', i+'_', threshold
        path_vcf = 'pulled_'+i+'_tobabachi.tsv'
        path_bed = 'pulled_'+i+'_tobabachi.bed'

        header_list = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'REF_COUNTS', 'ALT_COUNTS']
        vcf_data = pd.read_csv(os.path.join(processed_data, path_vcf), sep='\t',
                                    names=header_list)
        vcf_data['POS2'] = vcf_data['POS']
        vcf_data['POS2']+=1
        vcf_data = vcf_data[['#CHROM', 'POS', 'POS2','ID', 'REF', 'ALT', 'REF_COUNTS', 'ALT_COUNTS']]
        #vcf_list = vcf_data.values.tolist()
        test = BedTool.from_dataframe(vcf_data)

        bed_data = pd.read_csv(os.path.join(babachi, path_bed), sep='\t')
        bed_data = bed_data[['#chr','start','end','BAD']]
        #bed_list = bed_data.values.tolist()
        annotations = BedTool.from_dataframe(bed_data)
        inters = test.intersect(annotations, wb=True)
        df = inters.to_dataframe()
        df.columns = ['#CHROM', 'POS','POS2','ID', 'REF', 'ALT', 'REF_COUNTS', 'ALT_COUNTS','chr','start','end','BAD']
        annotated_vcf = df[['#CHROM', 'POS','ID', 'REF', 'ALT', 'REF_COUNTS', 'ALT_COUNTS','BAD']]
        annotated_vcf = annotated_vcf[(annotated_vcf.REF_COUNTS >= threshold) & (annotated_vcf.ALT_COUNTS >=threshold)]
        annotated_vcf.to_csv(os.path.join(fit, i + '_BAD_annotated.tsv'),header=True, index=False, sep='\t')
    #annotated_vcf = annotated_vcf[((annotated_vcf.ref + annotated_vcf.alt) >= threshold)]
    '''
    vcf_bad1 = annotated_vcf[annotated_vcf.BAD == 1]
    vcf_bad1 = vcf_bad1.drop(['bad'], axis=1)
    vcf_bad2 = annotated_vcf[annotated_vcf.BAD == 2]
    vcf_bad2 = vcf_bad2.drop(['bad'], axis=1)
    vcf_bad3 = annotated_vcf[annotated_vcf.BAD == 3]
    vcf_bad3 = vcf_bad3.drop(['bad'], axis=1)
    vcf_bad4 = annotated_vcf[annotated_vcf.BAD == 4]
    vcf_bad4 = vcf_bad4.drop(['bad'], axis=1)
    vcf_bad5 = annotated_vcf[annotated_vcf.BAD == 5]
    vcf_bad5 = vcf_bad5.drop(['bad'], axis=1)
    vcf_bad6 = annotated_vcf[annotated_vcf.BAD == 6]
    vcf_bad6 = vcf_bad6.drop(['bad'], axis=1)

    vcf_bad1.to_csv(os.path.join(processed_data, out_path + '_table_BAD1.tsv'), header=True, index=False, sep='\t')
    vcf_bad2.to_csv(os.path.join(processed_data, out_path + '_table_BAD2.tsv'), header=True, index=False, sep='\t')
    vcf_bad3.to_csv(os.path.join(processed_data, out_path + '_table_BAD3.tsv'), header=True, index=False, sep='\t')
    vcf_bad4.to_csv(os.path.join(processed_data, out_path + '_table_BAD4.tsv'), header=True, index=False, sep='\t')
    vcf_bad5.to_csv(os.path.join(processed_data, out_path + '_table_BAD5.tsv'), header=True, index=False, sep='\t')
    vcf_bad6.to_csv(os.path.join(processed_data, out_path + '_table_BAD6.tsv'), header=True, index=False, sep='\t')
    '''
    '''
    grp_bad1 = (vcf_bad1.groupby(['ref','alt']).size()
           .reset_index(name='counts'))
    grp_bad2 = (vcf_bad2.groupby(['ref','alt']).size()
           .reset_index(name='counts'))
    grp_bad3 = (vcf_bad3.groupby(['ref','alt']).size()
           .reset_index(name='counts'))
    grp_bad4 = (vcf_bad4.groupby(['ref','alt']).size()
           .reset_index(name='counts'))
    grp_bad5 = (vcf_bad5.groupby(['ref','alt']).size()
           .reset_index(name='counts'))
    grp_bad6 = (vcf_bad6.groupby(['ref','alt']).size()
           .reset_index(name='counts'))

    grp_bad1.to_csv(os.path.join(processed_data,out_path + '_BAD1.tsv'),header=True,index=False,sep='\t')
    grp_bad2.to_csv(os.path.join(processed_data, out_path + '_BAD2.tsv'),header=True,index=False,sep='\t')
    grp_bad3.to_csv(os.path.join(processed_data,out_path + '_BAD3.tsv'),header=True,index=False,sep='\t')
    grp_bad4.to_csv(os.path.join(processed_data, out_path + '_BAD4.tsv'),header=True,index=False,sep='\t')
    grp_bad5.to_csv(os.path.join(processed_data, out_path + '_BAD5.tsv'),header=True,index=False,sep='\t')
    grp_bad6.to_csv(os.path.join(processed_data,out_path + '_BAD6.tsv'),header=True,index=False,sep='\t')
    '''



def annotate_by_bad_myid(i,my_id,threshold):
    if (os.path.exists(os.path.join(fit, i+'_annotated/'+my_id + '.tsv'))):
        pass
    else:
        if(not os.path.exists(os.path.join(fit, i+'_annotated/'))):
            os.mkdir((os.path.join(fit, i+'_annotated/')))

        # 'pulled_'+i+'_tobabachi.tsv', 'pulled_'+i+'_tobabachi.bed', i+'_', threshold
        path_vcf = my_id + '/' + my_id + '_rs_nucli_getero_filtrated.tsv'
        path_bed = 'pulled_' + i + '_tobabachi.bed'

        header_list = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'REF_COUNTS', 'ALT_COUNTS']
        vcf_data = pd.read_csv(os.path.join(processed_data, path_vcf), sep='\t',
                               names=header_list)
        vcf_data['POS2'] = vcf_data['POS']
        vcf_data['POS2'] += 1
        vcf_data = vcf_data[['#CHROM', 'POS', 'POS2', 'ID', 'REF', 'ALT', 'REF_COUNTS', 'ALT_COUNTS']]
        # vcf_list = vcf_data.values.tolist()
        test = BedTool.from_dataframe(vcf_data)

        bed_data = pd.read_csv(os.path.join(babachi, path_bed), sep='\t')
        bed_data = bed_data[['#chr', 'start', 'end', 'BAD']]
        # bed_list = bed_data.values.tolist()
        annotations = BedTool.from_dataframe(bed_data)
        inters = test.intersect(annotations, wb=True)
        df = inters.to_dataframe()
        df.columns = ['#CHROM', 'POS', 'POS2', 'ID', 'REF', 'ALT', 'REF_COUNTS', 'ALT_COUNTS', 'chr', 'start', 'end',
                      'BAD']
        annotated_vcf = df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'REF_COUNTS', 'ALT_COUNTS', 'BAD']]
        annotated_vcf = annotated_vcf[(annotated_vcf.REF_COUNTS >= threshold) & (annotated_vcf.ALT_COUNTS >= threshold)]
        annotated_vcf.to_csv(os.path.join(fit, i+'_annotated/'+my_id + '.tsv'), header=True, index=False, sep='\t')
    #return(os.path.join(i+'_annotated/',my_id + '.tsv'))

def loggi(tmp_log,tmp_err,stdout,stderr,k):
    stdout.split('\n')
    stderr.split('\n')

    if (k=='a'):
        with open(tmp_log, "a") as log:
            log.write(stdout)
        with open(tmp_err, "a") as err:
            err.write(stderr)
    else:
        with open(tmp_log, "w") as log:
            log.write(stdout)
        with open(tmp_err, "w") as err:
            err.write(stderr)
    print('logged')

def negbinfit_ids(i,force_pam):

    bad = i
    ser = metadata[metadata['BADgroup'] == bad]
    ids = list(ser['ID'])
    ids_list = list(set(ids) & set(processing_list))
    txt_step1 = open(os.path.join(fit, 'processing_list_' + bad), "w")

    for my_id in ids_list:
        txt_step1.write(my_id + '.tsv' + "\n")
    txt_step1.close()




    tmp_log = os.path.join(fit,i+'_ids_log')
    tmp_err = os.path.join(fit,i+'_ids_err')
    with open(tmp_log, "w") as log:
        log.write('STARTING!\n')
    with open(tmp_err, "w") as err:
        err.write('STARTING!\n')

    if (os.path.exists(os.path.join(fit,i+"_fit_ids"))):
        with open(tmp_log, "a") as log:
            log.write(os.path.join(fit,i+"_fit_ids") + ' exists\n')
            print(os.path.join(fit,i+"_fit_ids") + ' exists')
    else:
        os.mkdir(os.path.join(fit,i+"_fit_ids"))
        os.chdir(os.path.join(fit, i+'_annotated/'))
        collect = f'negbin_fit collect -f {os.path.join(fit,"processing_list_"+i)} -O {os.path.join(fit,i+"_fit_ids")}'
        process = subprocess.Popen(shlex.split(collect),
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   universal_newlines=True
                                   )
        stderr, stdout = process.communicate()
        loggi(tmp_log, tmp_err, stdout, stderr, 'a')
        print('negbin_fit collect done')
        with open(tmp_log, "a") as log:
            log.write('negbin_fit collect done\n')
        fit_nb = f'negbin_fit -O {os.path.join(fit,i+"_fit_ids")} --visualize'
        process = subprocess.Popen(shlex.split(fit_nb),
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   universal_newlines=True
                                   )
        stderr, stdout = process.communicate()
        loggi(tmp_log, tmp_err, stdout, stderr, 'a')
        print('negbin_fit done')
        with open(tmp_log, "a") as log:
            log.write('negbin_fit done\n')

    if(os.path.exists(os.path.join(fit,i+ "_ids_pvals"))):
        with open(tmp_log, "a") as log:
            log.write(os.path.join(fit,i+ "_ids_pvals") + ' exists')
    else:
        os.mkdir(os.path.join(fit,i+ "_ids_pvals"))
        os.chdir(os.path.join(fit, i + '_annotated/'))
        calc_pval = f'calc_pval -f {os.path.join(fit,"processing_list_"+i)} -O {os.path.join(fit,i+ "_ids_pvals")} -w {os.path.join(fit,i+"_fit_ids")}'
        process = subprocess.Popen(shlex.split(calc_pval),
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   universal_newlines=True
                                   )
        stderr, stdout = process.communicate()
        loggi(tmp_log, tmp_err, stdout, stderr, 'a')
        print('calc_pval done')

    txt_step2 = open(os.path.join(fit, 'aggr_list_' + bad), "w")
    for my_id in ids_list:
        #txt_step1.write(my_id + '.tsv' + "\n")
        if(os.path.isfile(os.path.join(fit,i+ "_ids_pvals/") + my_id + '.pvalue_table')):
            txt_step2.write(my_id + '.pvalue_table' + "\n")
    txt_step2.close()

    if(os.path.exists(os.path.join(fit,i+ "ids_aggregated.tsv")))  :
        with open(tmp_log, "a") as log:
            log.write(os.path.join(fit,i+ "ids_aggregated.tsv") + ' exists')
    else:
        os.chdir(os.path.join(fit,i+ "_ids_pvals"))
        aggr = f'calc_pval aggregate -f {os.path.join(fit, "aggr_list_" + i)} -O {os.path.join(fit,i+ "ids_aggregated.tsv")}'
        process = subprocess.Popen(shlex.split(aggr),
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   universal_newlines=True
                                   )
        stderr, stdout = process.communicate()
        loggi(tmp_log, tmp_err, stdout, stderr, 'a')
        print('calc_pval aggregate done')



metadata = pd.read_csv(met,sep='\t')
grp = (metadata.groupby(['BADgroup']).size()
       .reset_index(name='count'))
grp = grp[grp.BADgroup!='.']
bad_list = []
for i in grp.iloc[:,0]:
    bad_list.append(i)
print(bad_list)
grp = (metadata.groupby(['ASBgroup']).size()
       .reset_index(name='count'))
grp = grp[grp.BADgroup!='.']
asb_list = []
for i in grp.iloc[:,0]:
    asb_list.append(i)
print(asb_list)

if (not os.path.exists(fit)):
    os.mkdir(fit)

processing_list = []
with open(processing_list_path, 'r') as fp:
    for line in fp:
        processing_list.append(line.strip())

for my_id in processing_list:
    for my_id in processing_list:
        ser = metadata[metadata['ID'] == my_id]
        k = ser.BADgroup.astype(str)
        i = k.iloc[0]
        if(i!='.'):
            annotate_by_bad_myid(i,my_id,threshold)

for i in bad_list:
    print(i+' start')
    negbinfit_ids(i,force)

    print(i+' end')


######################################################



