import subprocess32 as subprocess
import os
import pandas as pd
import configparser
import sys as sys
'''
maindir = '/media/ElissarDisk/ADASTRA/'
indir = maindir + 'data/'
processed_ref = maindir + 'processed_ref/genome-norm.fasta'
ref_vcf = maindir + 'reference/00-common_all.vcf.gz'

outdir = maindir + 'processed_data/'
logdir = maindir + 'logs/data_processing/'
javapars = '-Xmx12G -XX:ParallelGCThreads=4'
'''

def readConfig_SNP(path, dir):
    config = configparser.ConfigParser()
    config.read(path)

    """
    maindir = '/media/ElissarDisk/ADASTRA/'
    indir = maindir + 'data/'
    processed_ref = maindir + 'processed_ref/genome-norm.fasta'
    ref_vcf = maindir + 'reference/00-common_all.vcf.gz'

    outdir = maindir + 'processed_data/'
    logdir = maindir + 'logs/data_processing/'
    javapars = '-Xmx12G -XX:ParallelGCThreads=4'
    """
    # Читаем некоторые
    #    значения из конфиг. файла.
    maindir = config.get("Directories", "maindir")
    indir = config.get("Directories", "data_in")
    processed_ref = config.get("Files", "ref_out1")
    ref_vcf = config.get("Files", "ref_vcf")
    metadata = config.get("Files","metadata")
    data_out = config.get("Directories", "data_out")
    data_log = config.get("Directories", "data_log")

    javapars = config.get("Parameters", "JavaParameters")

    if (dir == 'maindir'):
        return (maindir)
    elif (dir == 'indir'):
        return (maindir + indir)
    elif (dir == 'processed_ref'):
        return (maindir + processed_ref)
    elif (dir == 'ref_vcf'):
        return (maindir + ref_vcf)
    elif (dir == 'outdir'):
        return (maindir + data_out)
    elif (dir == 'logdir'):
        return (maindir + data_log)
    elif (dir == 'javapars'):
        return (javapars)
    elif (dir == 'metadata'):
        return (metadata)

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

#___________snp for one file
def SNP_func(my_id):
    print('start SNP')
    tmp_log = logdir+my_id+'_log'
    tmp_err = logdir+my_id+'_err'
    # ______

    print('samtools index')
    process = subprocess.Popen(['samtools', 'index', indir + my_id + '.bam',
                                outdir + my_id + '/' + my_id + '.bai'],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               universal_newlines=True
                               )
    stderr, stdout = process.communicate()
    loggi(tmp_log, tmp_err, stdout, stderr, 'w')

    print('mk bai symlink')
    os.symlink(outdir + my_id + '/' + my_id + '.bai', indir + my_id + '.bai')

    print('samtools view -b')
    process = subprocess.Popen(['samtools', 'view', '-b', indir + my_id + '.bam',
                                'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
                                'chr11', 'chr12',
                                'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21',
                                'chr22', 'chrX', 'chrY',
                                '-o' + outdir + my_id + '/' + my_id + '_chop.bam'],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               universal_newlines=True)

    stderr, stdout = process.communicate()
    loggi(tmp_log, tmp_err, stdout, stderr, 'a')



    print('picard addorreplace')
    process = subprocess.Popen(['picard', 'AddOrReplaceReadGroups', 'I=' + outdir + my_id+'/'+my_id+'_chop.bam',
                                'O=' + outdir + my_id+'/'+my_id+'_formatted.bam',
                                'RGLB=lib1', 'RGPL=seq1', 'RGPU=unit1','RGSM=20', 'RGID=1'],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               universal_newlines=True)
    stderr, stdout = process.communicate()
    loggi(tmp_log, tmp_err, stdout, stderr, 'a')



    # ______
    print('picard markduplicates')
    process = subprocess.Popen(['picard', 'MarkDuplicates',
                                'I=' + outdir + my_id + '/' + my_id + '_formatted.bam',
                                'O=' + outdir + my_id + '/' + my_id + '_ready.bam',
                                'REMOVE_DUPLICATES=true',
                                'M=' + outdir + my_id + '/' + my_id + '_metrics.txt'],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               universal_newlines=True)
    stderr, stdout = process.communicate()
    loggi(tmp_log, tmp_err, stdout, stderr, 'a')

    # ______

    print('gatk baserecalibrator')
    process = subprocess.Popen(['gatk', 'BaseRecalibrator',
                                '--java-options',javapars,
                                '-R', processed_ref,
                                '-I', outdir + my_id + '/' + my_id + '_ready.bam',
                                '--known-sites', ref_vcf,
                                '-O', outdir + my_id + '/' + my_id +'.table'],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               universal_newlines=True)
    stderr, stdout = process.communicate()
    loggi(tmp_log, tmp_err, stdout, stderr, 'a')
    # ______

    print('gatk applyBQSR')
    process = subprocess.Popen(['gatk', 'ApplyBQSR',
                                '-R',  processed_ref,
                                '--java-options', javapars,
                                '-I', outdir + my_id + '/' + my_id + '_ready.bam',
                                '--bqsr-recal-file', outdir + my_id + '/' + my_id +'.table',
                                '-O', outdir + my_id + '/' + my_id +'_final.bam' ],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               universal_newlines=True)
    stderr, stdout = process.communicate()
    loggi(tmp_log, tmp_err, stdout, stderr, 'a')


    # ______
    print('gatk haplotypecaller')
    process = subprocess.Popen(['gatk', 'HaplotypeCaller',
                                '--java-options', javapars,
                                '-R', processed_ref,
                                '-I', outdir + my_id + '/' + my_id +'_final.bam',
                                '--dbsnp', ref_vcf,
                                '-O', outdir + my_id + '/' + my_id +'.vcf'],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               universal_newlines=True)
    stderr, stdout = process.communicate()
    loggi(tmp_log, tmp_err, stdout, stderr, 'a')


def stats(my_id,fini ):
    if (fini == 2 ):
        strf = outdir + my_id + '/' + my_id +'_ready.bam'
        statfilecov = outdir + my_id + '/'  + 'stats_nodup_cov.txt'
        statfile = outdir + my_id + '/'  +'stats_nidup_num.txt'

    elif (fini == 1 ):
        strf = outdir + my_id + '/' + my_id +'_final.bam'
        statfilecov = outdir + my_id + '/'  + 'stats_final_cov.txt'
        statfile = outdir + my_id + '/'  +'stats_final_num.txt'

    elif (fini == 0):
        strf = indir + my_id + '.bam'
        statfilecov = outdir + my_id + '/' + 'stats_start_cov.txt'
        statfile = outdir + my_id + '/' + 'stats_start_num.txt'
    elif (fini == 3): #vcf
        strf = outdir + my_id + '/' + my_id +'.vcf'
        statfile = outdir + my_id + '/' +'vcf_stats.txt'

    ##########start#########################################

    if (fini == 3): #vcf_number of peaks
        with open(strf) as f:
            text = f.readlines()
            size = len(text)
        with open(statfile, "w") as g:
            g.write('number of peaks '+ str(size))


    else:
        # __coverage###################################
        process = subprocess.run(['samtools', 'coverage',
                                    strf,
                                    '-o',statfilecov])

        # __number of reads##############################
        process = subprocess.run(['samtools', 'view', '-c',
                                     strf,
                                    '-o',statfile])

def rm(my_id):
    subprocess.run(['rm', outdir + my_id + '/' + my_id + '_final.bai',
                    outdir + my_id + '/' + my_id + '_final.bam'])



########my_id = 'BAM00030'#######################
def pipe_my_id(my_id):
    print('start SNP '+ my_id)
    #my_id = 'BAM00030'
    os.mkdir(outdir+my_id)

    print("Begginng with " + my_id)
    SNP_func(my_id)
    stats(my_id,0)
    stats(my_id,1)
    stats(my_id,2)
    stats(my_id,3)
    rm(my_id)



#path = "/media/ElissarDisk/ADASTRA/parameters/CONFIG.cfg"
path = sys.argv[1]
maindir = readConfig_SNP(path,'maindir')
indir = readConfig_SNP(path,'indir')
processed_ref = readConfig_SNP(path,'processed_ref')
ref_vcf = readConfig_SNP(path,'ref_vcf')
outdir = readConfig_SNP(path,'outdir')
logdir = readConfig_SNP(path,'logdir')
javapars = readConfig_SNP(path,'javapars')
met = readConfig_SNP(path,'metadata')
'''
met = maindir + 'parameters/metadata.tsv'
metadata = pd.read_csv(met,sep='\t')
norna = metadata[metadata["Extra1"] != 'RNA-seq']
idid = norna['ID']
idid = list(idid)
with open(maindir + 'parameters/idid', "w") as outfile:
    outfile.write("\n".join(idid))
'''
#for i in idid:
#    pipe_my_id(i)
my_id = sys.argv[2]
pipe_my_id(my_id)

#subprocess.run(['parallel', '-j', '4', 'step2_snp_calling.py','::::',maindir + 'parameters/idid'],capture_output=True)
