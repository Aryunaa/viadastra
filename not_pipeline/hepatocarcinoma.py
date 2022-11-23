import pandas as pd
import os
import sys
import configparser
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
from pybedtools import BedTool

#babachi
#get cnv id
# annotate with bads and ads

def writechrs_merg(df, i):
    fig, ax = plt.subplots(2, figsize=(13, 8))
    if (os.path.exists("./37/BAM00037.sns.badmap.visualization/BAM00037.sns.badmap_" + i + ".png")):
        im = plt.imread("./37/BAM00037.sns.badmap.visualization/BAM00037.sns.badmap_" + i + ".png")

        im = ax[0].imshow(im)
        ax[0].axis("off")
        ax[0].autoscale()

    plt.style.use('ggplot')

    for Start, End, ad in zip(df['Start'].values, df['End'].values, df['AD'].values):
        ax[1].add_patch(plt.Rectangle((Start, ad), End - Start, 0.3))

    ax[1].autoscale()
    ax[1].set_ylim(0.8, 6)
    ax[1].set_xlabel('Chromosome position, bp')
    ax[1].set_ylabel('AD', rotation='vertical')
    ax[1].set_title(i)

    # ax.set_xlim(df['Start'].min(), df['End'].max())
    # default_x_ticks = range(len(df['End']))
    # labels []
    # plt.xticks(x, labels, rotation ='vertical')
    plt.yticks([1, 2, 3, 4, 5, 6], rotation='horizontal')
    # plt.show()
    plt.tight_layout()
    fig.savefig('./37/adcn_viz/' + i + '.png')
    # plt.show()
def runbabachi(rssnpspath):

    return 0

def annotatebybadandad(rssnpspath, path_badmap, path_cnvs, out_path):
    rssnps = pd.read_csv(rssnpspath, sep='\t')
    rssnpslist = rssnps.values.tolist()
    test = BedTool(rssnpslist)
    print(rssnps.shape[0])
    # bed_data = pd.read_csv(os.path.join(processed_data,'pulled_atacseq_tobabachi.bed'),sep='\t')
    # badmaps
    bed_data_old = pd.read_csv(path_badmap, sep='\t')
    bed_data = bed_data_old[['#chr', 'start', 'end']]
    bed_data = bed_data[bed_data.end >= bed_data.start]
    bed_list = bed_data.values.tolist()
    print(bed_data.shape[0])
    annotations = BedTool(bed_list)
    i = test.intersect(annotations, wb=True)

    # cnvs
    bed2 = pd.read_csv(path_cnvs, sep='\t')
    bed2 = bed2[bed2['Minor_Copy_Number'] != 0]
    bed2['knownAD'] = bed2['Major_Copy_Number'] * 1.0 / bed2['Minor_Copy_Number']
    bed2data = bed2[['Chromosome', 'Start', 'End']]
    bed2data = bed2data[bed2data.End >= bed2data.Start]
    bed2_list = bed2data.values.tolist()
    annotations_2 = BedTool(bed2_list)
    i2 = i.intersect(annotations_2, wb=True)

    df = i2.to_dataframe()
    # print(df.iloc[0,:])
    df.columns = ['#CHROM', 'POS1', 'POS2', 'ID', 'REF', 'ALT', 'REF_COUNTS', 'ALT_COUNTS', 'SAMPLEID',
                  '#chr', 'start', 'end', 'Chromosome', 'Start', 'End']
    df = pd.merge(df, bed_data_old[['#chr', 'start', 'end', 'BAD']], how='left', on=['#chr', 'start', 'end'])
    df = pd.merge(df, bed2[['Chromosome', 'Start', 'End','knownAD']])

    print(df.shape[0])
    # annotated_vcf = annotated_vcf[(annotated_vcf.ref >= threshold) & (annotated_vcf.alt >=threshold)]
    # annotated_vcf = annotated_vcf[((annotated_vcf.ref + annotated_vcf.alt) >= threshold)]
    # df['POS2'] = df['POS'] + 1
    # df['POS']=df['POS2']
    annotated_snps = df[['#CHROM', 'POS1', 'POS2', 'ID', 'REF', 'ALT', 'REF_COUNTS', 'ALT_COUNTS', 'BAD', 'knownAD','SAMPLEID']]
    print('shape= ' + str(annotated_snps.shape[0]))
    annotated_snps.to_csv(out_path, header=True, index=False, sep='\t')
    print('annotated')
    return (annotated_snps)

config = configparser.ConfigParser()
pwd = os.getcwd()
config.read('/home/ariuna/viadastra_short/configs/petson_config.cfg')
maindir = config["Directories"]["maindir"]
babachi = '/home/ariuna/tcga_atacseq/babachi_all'
rssnpsdir = os.path.join(config["Directories"]["maindir"], config["Directories"]["rssnps"])
metapath = config["Files"]["metadata"]
metadata = pd.read_csv(metapath, sep = '\t')
cnvdata = pd.read_csv('/home/ariuna/viadastra_short/metas/stats_cnvstoo.tsv', sep = '\t')
cnvs_dir = '/home/ariuna/tcga_atacseq/cnv_data/'
idloc = cnvdata.columns.get_loc("ID")
ascid_loc = cnvdata.columns.get_loc("file_id_ascnv")
aspath_loc = cnvdata.columns.get_loc("file_name_ascnv")
cnvdata['tau'] = np.NaN
cnvdata['p_value'] = np.NaN
tauloc = cnvdata.columns.get_loc("tau")
ploc = cnvdata.columns.get_loc("p_value")
for i in range(cnvdata.shape[0]):
    print(cnvdata.iloc[i, idloc])
    id = cnvdata.iloc[i,idloc]
    print(cnvdata.shape[1] + ' '+ str(ascid_loc)+ " " + str(aspath_loc))
    cnvsfile = cnvdata.iloc[i,ascid_loc]+'/'+cnvdata.iloc[i,aspath_loc]
    path_cnvs = os.path.join(cnvs_dir,cnvsfile)
    print(path_cnvs)
    rssnpspath = os.path.join(rssnpsdir,id + '.snps.bed')
    pathbadmap = os.path.join(babachi,id + '.snps.badmap.bed')
    outputpath = '/home/ariuna/tcga_atacseq/processed_data/annotated_snps/' + id + '.snps.annotated.bed'
    if( os.path.exists(path_cnvs) and os.path.exists(rssnpspath) and os.path.exists(pathbadmap)):
        annotated_snps = annotatebybadandad(rssnpspath, pathbadmap, path_cnvs,outputpath)
        tau, p_value = stats.kendalltau(annotated_snps['knownAD'], annotated_snps['BAD'])
        cnvdata.iloc[i,tauloc] = tau
        cnvdata.iloc[i,ploc] = p_value


cnvdata.to_csv('/home/ariuna/tcga_atacseq/logs/stats.tsv',sep = '\t', index=False)





