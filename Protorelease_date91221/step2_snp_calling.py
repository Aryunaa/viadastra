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
    all_log = os.path.join(maindir, 'logs/whole_log')
    with open(all_log, "a") as log:
        log.write('STARTING! '+ my_id + '\n')

    tmp_log = logdir+my_id+'_log'
    tmp_err = logdir+my_id+'_err'
    with open(tmp_log, "w") as log:
        log.write('STARTING!')
    with open(tmp_err, "w") as err:
        err.write('STARTING!')

    # ______
    print('samtools sort')
    with open(all_log, "a") as log:
        log.write('samtools sort'+ '\n')
    if(os.path.exists(os.path.join(outdir,my_id) + '/' + my_id + '_sortsam')):
        print('already exists, go further')
        with open(all_log, "a") as log:
            log.write('already exists, go further'+ '\n')
    else:
        process = subprocess.Popen(['samtools', 'sort', indir + my_id + '.bam',
                                  '-o',os.path.join(outdir,my_id) + '/' + my_id + '_sortsam'],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   universal_newlines=True
                                   )
        stderr, stdout = process.communicate()
        loggi(tmp_log, tmp_err, stdout, stderr, 'a')
        print('done')

    # ______
    print('samtools index')
    with open(all_log, "a") as log:
        log.write('samtools index'+ '\n')
    if(os.path.exists(os.path.join(outdir,my_id) + '/' + my_id + '_sortsam.bai')):
        print('already exists, go further')
        with open(all_log, "a") as log:
            log.write('already exists, go further'+ '\n')
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
        with open(all_log, "a") as log:
            log.write('done'+ '\n')

    # ______
    print('samtools view -b')
    with open(all_log, "a") as log:
        log.write('samtools view -b'+ '\n')
    if(os.path.exists(outdir + my_id + '/' + my_id + '_chop.bam')):
        print('already exists, go further')
        with open(all_log, "a") as log:
            log.write('already exists, go further'+ '\n')
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
        with open(all_log, "a") as log:
            log.write('done'+ '\n')

    # ______
    print('picard SortSam')
    with open(all_log, "a") as log:
        log.write('picard SortSam'+ '\n')
    if(os.path.exists(outdir + my_id + '/' + my_id + '_sorted.bam')):
        print('already exists, go further')
        with open(all_log, "a") as log:
            log.write('already exists, go further'+ '\n')
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
        with open(all_log, "a") as log:
            log.write('done'+ '\n')


    # ______
    print('picard addorreplace')
    with open(all_log, "a") as log:
        log.write('picard addorreplace'+ '\n')
    if(os.path.exists(outdir + my_id+'/'+my_id+'_formatted.bam')):
        print('already exists, go further')
        with open(all_log, "a") as log:
            log.write('already exists, go further'+ '\n')
    else:
        process = subprocess.Popen(['picard', 'AddOrReplaceReadGroups', 'I=' + outdir + my_id+'/'+my_id+'_sorted.bam',
                                    'O=' + outdir + my_id+'/'+my_id+'_formatted.bam', 'VALIDATION_STRINGENCY=LENIENT',
                                    'RGLB=lib1', 'RGPL=seq1', 'RGPU=unit1','RGSM=20', 'RGID=1'],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   universal_newlines=True)
        stderr, stdout = process.communicate()
        loggi(tmp_log, tmp_err, stdout, stderr, 'a')
        print('done')
        with open(all_log, "a") as log:
            log.write('done'+ '\n')

    # ______
    print('picard markduplicates')
    with open(all_log, "a") as log:
        log.write('picard markduplicates'+ '\n')
    if(os.path.exists(outdir + my_id + '/' + my_id + '_ready.bam')):
        print('already exists, go further')
        with open(all_log, "a") as log:
            log.write('already exists, go further'+ '\n')
    else:
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
        print('done')
        with open(all_log, "a") as log:
            log.write('done'+ '\n')

    # ______

    print('gatk baserecalibrator')
    with open(all_log, "a") as log:
        log.write('gatk baserecalibrator'+ '\n')
    if(os.path.exists(outdir + my_id + '/' + my_id +'.table')):
        print('already exists, go further')
        with open(all_log, "a") as log:
            log.write('already exists, go further'+ '\n')
    else:
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
        print('done')
        with open(all_log, "a") as log:
            log.write('done'+ '\n')
    # ______

    print('gatk applyBQSR')
    with open(all_log, "a") as log:
        log.write('gatk applyBQSR'+ '\n')
    if(os.path.exists(outdir + my_id + '/' + my_id +'_final.bam')):
        print('already exists, go further')
        with open(all_log, "a") as log:
            log.write('already exists, go further'+ '\n')
    else:
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
        print('done')
        with open(all_log, "a") as log:
            log.write('done'+ '\n')


    # ______
    print('gatk haplotypecaller')
    with open(all_log, "a") as log:
        log.write('gatk haplotypecaller'+ '\n')
    if(os.path.exists(outdir + my_id + '/' + my_id +'.vcf')):
        print('already exists, go further')
        with open(all_log, "a") as log:
            log.write('already exists, complete'+ '\n')
    else:
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
        print('done')
        with open(all_log, "a") as log:
            log.write('processing done with '+ my_id+ '\n')


def stats(my_id,clause ):
    if (clause == '_ready.bam' ):
        strf = outdir + my_id + '/' + my_id +'_ready.bam'
        statfilecov = outdir + my_id + '/'  + 'stats_nodup_cov.txt'
        statfile = outdir + my_id + '/'  +'stats_nodup_num.txt'

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
        subprocess.run(['samtools', 'coverage',
                                    strf,
                                    '-o',statfilecov])

        # __number of reads##############################
        subprocess.run(['samtools', 'view', '-c',
                                     strf,
                                    '-o',statfile])

def rm(my_id):
    rm_list = [os.path.join(outdir, my_id) + '/' + my_id + '_sortsam',
                    os.path.join(outdir, my_id) + '/' + my_id + '_sortsam.bai',
                    outdir + my_id + '/' + my_id + '_chop.bam',
                    outdir + my_id + '/' + my_id + '_sorted.bam',
                    outdir + my_id + '/' + my_id + '_formatted.bam',
                    outdir + my_id + '/' + my_id + '_ready.bam',
                    outdir + my_id + '/' + my_id + '.table',
                    outdir + my_id + '/' + my_id + '_final.bai',
                    outdir + my_id + '/' + my_id + '_final.bam']
    for i in rm_list:
        if(os.path.exists(i)):
            os.remove(i)




########my_id = 'BAM00030'#######################
def pipe_my_id(my_id):
    print('start SNP '+ my_id)
    #my_id = 'BAM00030'
    tmp_path = os.path.join(outdir,my_id)
    if (os.path.exists(tmp_path)):
        pass
    else:
        os.mkdir(tmp_path)
    if (os.path.exists(outdir + my_id + '/' + my_id + '.vcf')):
        print(my_id+ ' vcf file from gatk already exists, skip processing')
    else:
        process_bam(my_id)

    #print("Beginnig with " + my_id)
    #process_bam(my_id)
    #stats(my_id,'.bam')
    #stats(my_id,'_final.bam')
    #stats(my_id,'_ready.bam')
    stats(my_id,'.vcf')
    rm(my_id)

#path = "/media/ElissarDisk/ADASTRA/parameters/CONFIG.cfg"
# reading config ------------------------------------------
path = sys.argv[1]
'''
dicti = readConfig_SNP(path)
maindir = dicti['maindir']
indir = dicti['indir']

processed_ref = dicti['processed_ref']
ref_vcf = dicti['ref_vcf']
outdir = dicti['outdir']
logdir = dicti['logdir']
javapars = dicti['javapars']
met = dicti['metadata']
'''
config = configparser.ConfigParser()
config.read(path)
maindir = config["Directories"]["maindir"]
print(maindir)
indir = os.path.join(maindir,config["Directories"]["data_in"])
print(indir)
processed_ref = os.path.join(maindir,config["Files"]["ref_out1"])
ref_vcf = os.path.join(maindir,config["Files"]["ref_vcf"])
outdir = os.path.join(maindir,config["Directories"]["data_out"])
logdir = os.path.join(maindir,config["Directories"]["data_log"])
javapars = config["Parameters"]["javaparameters"]
met = os.path.join(maindir,config["Files"]["metadata"])

# pipe ----------------------------------------------------
my_id = sys.argv[2]

pipe_my_id(my_id)
