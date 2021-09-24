import subprocess32 as subprocess
import os
from config import readConfig_SNP
import configparser
import sys as sys

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
def process_bam(my_id):
    print('start SNP')
    tmp_log = logdir+my_id+'_log'
    tmp_err = logdir+my_id+'_err'
    with open(tmp_log, "w") as log:
        log.write('STARTING!')
    with open(tmp_err, "w") as err:
        err.write('STARTING!')

    # ______
    print('samtools sort')
    if(os.path.exists(os.path.join(outdir,my_id) + '/' + my_id + '_sortsam')):
        print('already exists, go further')
    else:
        process = subprocess.Popen(['samtools', 'sort', indir + my_id + '.bam',
                                    '>', os.path.join(outdir,my_id) + '/' + my_id + '_sortsam'],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   universal_newlines=True
                                   )
        stderr, stdout = process.communicate()
        loggi(tmp_log, tmp_err, stdout, stderr, 'a')
    print('done')

    # ______
    print('samtools index')
    if(os.path.exists(os.path.join(outdir,my_id) + '/' + my_id + '_sortsam.bai')):
        print('already exists, go further')
    else:
        process = subprocess.Popen(['samtools', 'index',os.path.join(outdir,my_id) + '/' + my_id + '_sortsam',
                                    os.path.join(outdir,my_id) + '/' + my_id + '_sortsam.bai'],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   universal_newlines=True
                                   )
        stderr, stdout = process.communicate()
        loggi(tmp_log, tmp_err, stdout, stderr, 'a')
    print('done')

    # ______
    print('samtools view -b')
    if(os.path.exists(outdir + my_id + '/' + my_id + '_chop.bam')):
        print('already exists, go further')
    else:
        process = subprocess.Popen(['samtools', 'view', '-b', os.path.join(outdir,my_id) + '/' + my_id + '_sortsam',
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
    print('done')

    # ______
    print('picard SortSam')
    if(os.path.exists(outdir + my_id + '/' + my_id + '_sorted.bam')):
        print('already exists, go further')
    else:
        process = subprocess.Popen(['picard', 'SortSam', 'I=' + outdir + my_id + '/' + my_id + '_chop.bam',
                                    'O=' + outdir + my_id + '/' + my_id + '_sorted.bam',
                                    'SORT_ORDER=coordinate','VALIDATION_STRINGENCY=LENIENT'],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   universal_newlines=True)
        stderr, stdout = process.communicate()
        loggi(tmp_log, tmp_err, stdout, stderr, 'a')
    print('done')

    # ______
    print('picard addorreplace')
    if(os.path.exists(outdir + my_id+'/'+my_id+'_formatted.bam')):
        print('already exists, go further')
    else:
        process = subprocess.Popen(['picard', 'AddOrReplaceReadGroups', 'I=' + outdir + my_id+'/'+my_id+'_sorted.bam',
                                    'O=' + outdir + my_id+'/'+my_id+'_formatted.bam', 'VALIDATION_STRINGENCY=LENIENT',
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
                                'REMOVE_DUPLICATES=true','VALIDATION_STRINGENCY=LENIENT',
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


def stats(my_id,clause ):
    if (clause == '_ready.bam' ):
        strf = outdir + my_id + '/' + my_id +'_ready.bam'
        statfilecov = outdir + my_id + '/'  + 'stats_nodup_cov.txt'
        statfile = outdir + my_id + '/'  +'stats_nidup_num.txt'

    elif (clause == '_final.bam' ):
        strf = outdir + my_id + '/' + my_id +'_final.bam'
        statfilecov = outdir + my_id + '/'  + 'stats_final_cov.txt'
        statfile = outdir + my_id + '/'  +'stats_final_num.txt'

    elif (clause == '.bam'):
        strf = indir + my_id + '.bam'
        statfilecov = outdir + my_id + '/' + 'stats_start_cov.txt'
        statfile = outdir + my_id + '/' + 'stats_start_num.txt'
    elif (clause == '.vcf'): #vcf
        strf = outdir + my_id + '/' + my_id +'.vcf'
        statfile = outdir + my_id + '/' +'vcf_stats.txt'

    ##########start#########################################

    if (clause == '.vcf'): #vcf_number of peaks
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
    tmp_path = os.path.join(outdir,my_id)
    if (os.path.exists(tmp_path)):
        process_bam(my_id)
    else:
        os.mkdir(tmp_path)
        process_bam(my_id)

    print("Beginnig with " + my_id)
    process_bam(my_id)
    stats(my_id,'.bam')
    stats(my_id,'_final.bam')
    stats(my_id,'_ready.bam')
    stats(my_id,'.vcf')
    rm(my_id)



#path = "/media/ElissarDisk/ADASTRA/parameters/CONFIG.cfg"
# reading config ------------------------------------------
path = sys.argv[1]

dicti = readConfig_SNP(path)

maindir = dicti['maindir']
print(maindir)
indir = dicti['indir']
print(indir)
processed_ref = dicti['processed_ref']
ref_vcf = dicti['ref_vcf']
outdir = dicti['outdir']
logdir = dicti['logdir']
javapars = dicti['javapars']
met = dicti['metadata']


# pipe ----------------------------------------------------
my_id = sys.argv[2]
pipe_my_id(my_id)