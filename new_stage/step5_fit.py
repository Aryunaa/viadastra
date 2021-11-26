import pandas as pd
import os
import pysam
import sys
import configparser
import pybedtools
from pybedtools import BedTool
import subprocess32 as subprocess
#path = '/media/ElissarDisk/ADASTRA/parameters/CONFIG.cfg'
threshold = int(sys.argv[1])
path = sys.argv[2]
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

def negbinfit(out_path,badn,badt, threshold):
    #badt = '_BAD1.tsv'
    print(os.path.join(processed_data,out_path + badt))
    subprocess.run(['negbin_fit', os.path.join(processed_data,out_path + badt),
                                '--bad',str(badn), '--allele-reads-tr', str(threshold),
                                '-O',os.path.join(processed_data,'grouped_bads/fitted_negbin'), '--visualize', '-r'
                              ],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               universal_newlines=True
                               )
    subprocess.run(['negbin_fit', os.path.join(processed_data,out_path + badt),
                                '--bad',str(badn), '--allele-reads-tr', str(threshold),
                                '-O',os.path.join(processed_data,'grouped_bads/fitted_negbin_lin'), '--visualize', '-r',
                              ],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               universal_newlines=True
                               )
    subprocess.run(['negbin_fit', os.path.join(processed_data,out_path + badt),
                                '--bad',str(badn), '--allele-reads-tr', str(threshold),
                                '-O',os.path.join(processed_data,'grouped_bads/fitted_negbin_lin'), '--visualize', '-r','-l'
                              ],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               universal_newlines=True
                               )



if(not os.path.isdir(os.path.join(processed_data,'grouped_bads'))):
    os.mkdir(os.path.join(processed_data,'grouped_bads'))

metadata = pd.read_csv(met,sep='\t')
grp = (metadata.groupby(['BADgroup']).size()
       .reset_index(name='count'))
grp = grp[grp.BADgroup!='.']
bad_list = []
for i in grp.iloc[:,0]:
    bad_list.append(i)
print(bad_list)

#group_by_bad('pulled_atacseq_tobabachi.tsv','pulled_atacseq_tobabachi.bed', 'grouped_bads/atacseq_',threshold)
if (not os.path.exists(fit)):
    os.mkdir(fit)

for i in bad_list:
    print(i+' start')
    #processed data, babachi, fit
    annotate_by_bad(i, threshold)
    '''
    if (not os.path.exists(os.path.join(fit,i+'/'))):
        os.mkdir(fit)
    subprocess.run(['negbin_fit', 'collect','-I',
                    os.path.join(fit, i + '_BAD_annotated.tsv'),
                    '-O', os.path.join(fit,i+'/')
                    ],
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   universal_newlines=True
                   )
    subprocess.run(['negbin_fit',
                    '-O', os.path.join(fit,i+'/'),
                    '--visualize'
                    ],
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   universal_newlines=True
                   )
    '''
    print(i+' end')



######################################################




'''
for i in range(6):
    print(str(i+1)+'_BAD.tsv')
    negbinfit('grouped_bads/atacseq_altref_counts',i+1,'_BAD'+str(i+1)+'.tsv',threshold)
    negbinfit('grouped_bads/chipseq_altref_counts', i + 1, '_BAD' + str(i + 1) + '.tsv', threshold)
'''


