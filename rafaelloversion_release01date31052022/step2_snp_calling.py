import subprocess32 as subprocess
import os
import configparser
import sys as sys
import shutil
import pandas as pd

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
        log.write('samtools sort '+my_id + '\n')

    if(os.path.exists(os.path.join(outdir,my_id) + '/' + my_id + '_sortsam')):
        print('already exists, go further')
        with open(all_log, "a") as log:
            log.write(os.path.join(outdir,my_id) + '/' + my_id + '_sortsam'+' already exists, go further'+ '\n')
    else:
        try:
            process = subprocess.Popen(['samtools', 'sort', indir + my_id + '.bam',
                                      '-o',os.path.join(outdir,my_id) + '/' + my_id + '_sortsam'],
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE,
                                       universal_newlines=True
                                       )
            stderr, stdout = process.communicate()
            loggi(tmp_log, tmp_err, stdout, stderr, 'a')
            rc = process.returncode
            print('done')
            if (rc==0):
                with open(all_log, "a") as log:
                    log.write('samtools sort done with '+my_id + '\n')
            else:
                with open(all_log, "a") as log:
                    log.write('samtools sort failed with '+my_id + '\n')
                if (os.path.exists(os.path.join(outdir, my_id) + '/' + my_id + '_sortsam')):
                    os.remove(os.path.join(outdir, my_id) + '/' + my_id + '_sortsam')
                sys.exit(1)
        except Exception:
            with open(all_log, "a") as log:
                log.write('samtools sort failed with ' + my_id + '\n')
            if (os.path.exists(os.path.join(outdir, my_id) + '/' + my_id + '_sortsam')):
                os.remove(os.path.join(outdir, my_id) + '/' + my_id + '_sortsam')
            sys.exit(10)

    # ______
    print('samtools index')
    with open(all_log, "a") as log:
        log.write('samtools index '+ my_id +'\n')
    if(os.path.exists(os.path.join(outdir,my_id) + '/' + my_id + '_sortsam.bai')):
        print('already exists, go further')
        with open(all_log, "a") as log:
            log.write(os.path.join(outdir,my_id) + '/' + my_id + '_sortsam.bai'+' already exists, go further'+ '\n')
    else:
        try:
            process = subprocess.Popen(['samtools', 'index',os.path.join(outdir,my_id) + '/' + my_id + '_sortsam',
                                        os.path.join(outdir,my_id) + '/' + my_id + '_sortsam.bai'],
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE,
                                       universal_newlines=True
                                       )
            stderr, stdout = process.communicate()
            loggi(tmp_log, tmp_err, stdout, stderr, 'a')
            rc = process.returncode
            if (rc == 0):
                with open(all_log, "a") as log:
                    log.write('samtools index done with ' + my_id + '\n')
            else:
                with open(all_log, "a") as log:
                    log.write('samtools index failed with ' + my_id + '\n')
                if (os.path.exists(os.path.join(outdir, my_id) + '/' + my_id + '_sortsam.bai')):
                    os.remove(os.path.join(outdir, my_id) + '/' + my_id + '_sortsam.bai')
                sys.exit(2)
            print('done')
        except Exception:
            with open(all_log, "a") as log:
                log.write('samtools index failed with ' + my_id + '\n')
            if (os.path.exists(os.path.join(outdir, my_id) + '/' + my_id + '_sortsam.bai')):
                os.remove(os.path.join(outdir, my_id) + '/' + my_id + '_sortsam.bai')


    # ______
    print('samtools view -b')
    with open(all_log, "a") as log:
        log.write('samtools view -b '+my_id + '\n')
    if(os.path.exists(os.path.join(outdir,my_id) + '/' + my_id + '_chop.bam')):
        print('already exists, go further')
        with open(all_log, "a") as log:
            log.write(os.path.join(outdir,my_id) + '/' + my_id + '_chop.bam'+' already exists, go further'+ '\n')
    else:
        try:
            process = subprocess.Popen(['samtools', 'view', '-b', os.path.join(outdir,my_id) + '/' + my_id + '_sortsam',
                                        'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
                                        'chr11', 'chr12',
                                        'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21',
                                        'chr22', 'chrX', 'chrY',
                                        '-o' + os.path.join(outdir,my_id) + '/' + my_id + '_chop.bam'],
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE,
                                       universal_newlines=True)

            stderr, stdout = process.communicate()
            loggi(tmp_log, tmp_err, stdout, stderr, 'a')
            print('done')
            rc = process.returncode
            if (rc == 0):
                with open(all_log, "a") as log:
                    log.write('samtools view done with ' + my_id + '\n')
            else:
                with open(all_log, "a") as log:
                    log.write('samtools view failed with ' + my_id + '\n')
                if (os.path.exists(os.path.join(outdir, my_id) + '/' + my_id + '_chop.bam')):
                    os.remove(os.path.join(outdir, my_id) + '/' + my_id + '_chop.bam')
                sys.exit(3)
        except Exception:
            with open(all_log, "a") as log:
                log.write('samtools view failed with ' + my_id + '\n')
            if (os.path.exists(os.path.join(outdir, my_id) + '/' + my_id + '_chop.bam')):
                os.remove(os.path.join(outdir, my_id) + '/' + my_id + '_chop.bam')

    #удаление уже ненужных файлов
    if (os.path.exists(os.path.join(outdir,my_id) + '/' + my_id + '_sortsam')):
        os.remove(os.path.join(outdir,my_id) + '/' + my_id + '_sortsam')
    if (os.path.exists(os.path.join(outdir,my_id) + '/' + my_id + '_sortsam.bai')):
        os.remove(os.path.join(outdir,my_id) + '/' + my_id + '_sortsam.bai')

    # ______
    print('picard SortSam')
    with open(all_log, "a") as log:
        log.write('picard SortSam '+my_id + '\n')
    if(os.path.exists(os.path.join(outdir,my_id) + '/' + my_id + '_sorted.bam')):
        print('already exists, go further')
        with open(all_log, "a") as log:
            log.write(os.path.join(outdir,my_id) + '/' + my_id + '_sorted.bam'+'already exists, go further'+ '\n')
    else:
        try:
            process = subprocess.Popen([ 'picard', 'SortSam', 'I=' + os.path.join(outdir,my_id) + '/' + my_id + '_chop.bam',
                                        'O=' + os.path.join(outdir,my_id) + '/' + my_id + '_sorted.bam',
                                        'SORT_ORDER=coordinate','VALIDATION_STRINGENCY=LENIENT', 'TMP_DIR=',os.path.join(outdir,my_id)],
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE,
                                       universal_newlines=True)
            stderr, stdout = process.communicate()
            loggi(tmp_log, tmp_err, stdout, stderr, 'a')
            print('done')
            rc = process.returncode
            if (rc == 0):
                with open(all_log, "a") as log:
                    log.write('picard SortSam done with ' + my_id + '\n')
            else:
                with open(all_log, "a") as log:
                    log.write('picard SortSam failed with ' + my_id + '\n')
                if (os.path.exists(os.path.join(outdir, my_id) + '/' + my_id + '_sorted.bam')):
                    os.remove(os.path.join(outdir, my_id) + '/' + my_id + '_sorted.bam')
                sys.exit(4)
        except Exception:
            with open(all_log, "a") as log:
                log.write('picard SortSam failed with ' + my_id + '\n')
            if (os.path.exists(os.path.join(outdir, my_id) + '/' + my_id + '_sorted.bam')):
                os.remove(os.path.join(outdir, my_id) + '/' + my_id + '_sorted.bam')
    # удаление уже ненужных файлов
    if (os.path.exists(os.path.join(outdir,my_id) + '/' + my_id + '_chop.bam')):
        os.remove(os.path.join(outdir,my_id) + '/' + my_id + '_chop.bam')

    # ______
    print('picard addorreplace')
    with open(all_log, "a") as log:
        log.write('picard addorreplace '+my_id + '\n')
    if(os.path.exists(os.path.join(outdir,my_id)+'/'+my_id+'_formatted.bam')):
        print(os.path.join(outdir,my_id)+'/'+my_id+'_formatted.bam'+' already exists, go further')
        with open(all_log, "a") as log:
            log.write(os.path.join(outdir,my_id)+'/'+my_id+'_formatted.bam'+' already exists, go further'+ '\n')
    else:
        try:
            process = subprocess.Popen(['picard', 'AddOrReplaceReadGroups', 'I=' + os.path.join(outdir,my_id) +'/'+my_id+'_sorted.bam',
                                        'O=' + os.path.join(outdir,my_id) +'/'+my_id+'_formatted.bam', 'VALIDATION_STRINGENCY=LENIENT',
                                        'RGLB=lib1', 'RGPL=seq1', 'RGPU=unit1','RGSM=20', 'RGID=1','TMP_DIR=',os.path.join(outdir,my_id)],
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE,
                                       universal_newlines=True)
            stderr, stdout = process.communicate()
            loggi(tmp_log, tmp_err, stdout, stderr, 'a')
            print('done')
            rc = process.returncode
            if (rc == 0):
                with open(all_log, "a") as log:
                    log.write('picard addorreplace done with ' + my_id + '\n')
            else:
                with open(all_log, "a") as log:
                    log.write('picard addorreplace failed with ' + my_id + '\n')
                if (os.path.exists(os.path.join(outdir, my_id) + '/' + my_id + '_formatted.bam')):
                    os.remove(os.path.join(outdir, my_id) + '/' + my_id + '_formatted.bam')
                sys.exit(5)
        except Exception:
            with open(all_log, "a") as log:
                log.write('picard addorreplace failed with ' + my_id + '\n')
            if (os.path.exists(os.path.join(outdir, my_id) + '/' + my_id + '_formatted.bam')):
                os.remove(os.path.join(outdir, my_id) + '/' + my_id + '_formatted.bam')

    # удаление уже ненужных файлов
    if (os.path.exists(os.path.join(outdir,my_id) + '/' + my_id + '_sorted.bam')):
        os.remove(os.path.join(outdir,my_id) + '/' + my_id + '_sorted.bam')


    # ______
    print('picard markduplicates')
    with open(all_log, "a") as log:
        log.write('picard markduplicates '+my_id +'\n')
    if(os.path.exists(os.path.join(outdir,my_id) + '/' + my_id + '_ready.bam')):
        print('already exists, go further')
        with open(all_log, "a") as log:
            log.write(os.path.join(outdir,my_id) + '/' + my_id + '_ready.bam'+' already exists, go further'+ '\n')
    else:
        try:
            process = subprocess.Popen(['picard', 'MarkDuplicates',
                                        'I=' + os.path.join(outdir,my_id) + '/' + my_id + '_formatted.bam',
                                        'O=' + os.path.join(outdir,my_id) + '/' + my_id + '_ready.bam',
                                        'REMOVE_DUPLICATES=true','VALIDATION_STRINGENCY=LENIENT',
                                        'M=' + os.path.join(outdir,my_id) + '/' + my_id + '_metrics.txt','TMP_DIR=',os.path.join(outdir,my_id)],
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE,
                                       universal_newlines=True)
            stderr, stdout = process.communicate()
            loggi(tmp_log, tmp_err, stdout, stderr, 'a')
            print('done')
            rc = process.returncode
            if (rc == 0):
                with open(all_log, "a") as log:
                    log.write('picard markduplicates done with ' + my_id + '\n')
            else:
                with open(all_log, "a") as log:
                    log.write('picard markduplicates failed with ' + my_id + '\n')
                if (os.path.exists(os.path.join(outdir, my_id) + '/' + my_id + '_ready.bam')):
                    os.remove(os.path.join(outdir, my_id) + '/' + my_id + '_ready.bam')
                sys.exit(6)
        except Exception:
            with open(all_log, "a") as log:
                log.write('picard markduplicates failed with ' + my_id + '\n')
            if (os.path.exists(os.path.join(outdir, my_id) + '/' + my_id + '_ready.bam')):
                os.remove(os.path.join(outdir, my_id) + '/' + my_id + '_ready.bam')

    #удаление уже ненужных файлов
    if (os.path.exists(os.path.join(outdir,my_id) + '/' + my_id + '_formatted.bam')):
        os.remove(os.path.join(outdir,my_id) + '/' + my_id + '_formatted.bam')

    # ______

    print('gatk baserecalibrator')
    with open(all_log, "a") as log:
        log.write('gatk baserecalibrator '+my_id + '\n')
    if(os.path.exists(os.path.join(outdir,my_id) + '/' + my_id +'.table')):
        print('already exists, go further')
        with open(all_log, "a") as log:
            log.write(os.path.join(outdir,my_id) + '/' + my_id +'.table'+' already exists, go further'+ '\n')
    else:
        try:
            process = subprocess.Popen(['gatk', 'BaseRecalibrator',
                                        '--java-options',javapars,
                                        '-R', processed_ref,
                                        '-I', os.path.join(outdir,my_id) + '/' + my_id + '_ready.bam',
                                        '--known-sites', ref_vcf,
                                        '-O', os.path.join(outdir,my_id) + '/' + my_id +'.table' ],
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE,
                                       universal_newlines=True)
            stderr, stdout = process.communicate()
            loggi(tmp_log, tmp_err, stdout, stderr, 'a')
            print('done')
            rc = process.returncode
            if (rc == 0):
                with open(all_log, "a") as log:
                    log.write('gatk baserecalibrator done with ' + my_id + '\n')
            else:
                with open(all_log, "a") as log:
                    log.write('gatk baserecalibrator failed with ' + my_id + '\n')
                if (os.path.exists(os.path.join(outdir, my_id) + '/' + my_id + '.table')):
                    os.remove(os.path.join(outdir, my_id) + '/' + my_id + '.table')
                sys.exit(7)
        except Exception:
            with open(all_log, "a") as log:
                log.write('gatk baserecalibrator failed with ' + my_id + '\n')
            if (os.path.exists(os.path.join(outdir, my_id) + '/' + my_id + '.table')):
                os.remove(os.path.join(outdir, my_id) + '/' + my_id + '.table')
    # ______

    print('gatk applyBQSR')
    with open(all_log, "a") as log:
        log.write('gatk applyBQSR '+my_id + '\n')
    if(os.path.exists(os.path.join(outdir,my_id) + '/' + my_id +'_final.bam')):
        print('already exists, go further')
        with open(all_log, "a") as log:
            log.write(os.path.join(outdir,my_id) + '/' + my_id +'_final.bam'+' already exists, go further'+ '\n')
    else:
        try:
            process = subprocess.Popen(['gatk', 'ApplyBQSR',
                                        '-R',  processed_ref,
                                        '--java-options', javapars,
                                        '-I', os.path.join(outdir,my_id) + '/' + my_id + '_ready.bam',
                                        '--bqsr-recal-file', os.path.join(outdir,my_id) + '/' + my_id +'.table',
                                        '-O', os.path.join(outdir,my_id) + '/' + my_id +'_final.bam'],
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE,
                                       universal_newlines=True)
            stderr, stdout = process.communicate()
            loggi(tmp_log, tmp_err, stdout, stderr, 'a')
            print('done')
            rc = process.returncode
            if (rc == 0):
                with open(all_log, "a") as log:
                    log.write('gatk applyBQSR done with ' + my_id + '\n')
            else:
                with open(all_log, "a") as log:
                    log.write('gatk applyBQSR failed with ' + my_id + '\n')
                if (os.path.exists(os.path.join(outdir, my_id) + '/' + my_id + '_final.bam')):
                    os.remove(os.path.join(outdir, my_id) + '/' + my_id + '_final.bam')
                sys.exit(9)
        except Exception:
            with open(all_log, "a") as log:
                log.write('gatk applyBQSR failed with ' + my_id + '\n')
            if (os.path.exists(os.path.join(outdir, my_id) + '/' + my_id + '_final.bam')):
                os.remove(os.path.join(outdir, my_id) + '/' + my_id + '_final.bam')

    # удаление уже ненужных файлов
    if (os.path.exists(os.path.join(outdir,my_id) + '/' + my_id + '_ready.bam')):
        os.remove(os.path.join(outdir,my_id) + '/' + my_id + '_ready.bam')

    # ______
    print('gatk haplotypecaller')
    with open(all_log, "a") as log:
        log.write('gatk haplotypecaller '+ my_id +'\n')
    if(os.path.exists(os.path.join(outdir,my_id) + '/' + my_id +'.vcf')):
        print('already exists, go further')
        with open(all_log, "a") as log:
            log.write(os.path.join(outdir,my_id) + '/' + my_id +'.vcf'+' already exists, complete'+ '\n')
    else:
        try:
            process = subprocess.Popen(['gatk', 'HaplotypeCaller',
                                        '--java-options', javapars,
                                        '-R', processed_ref,
                                        '-I', os.path.join(outdir,my_id) + '/' + my_id +'_final.bam',
                                        '--dbsnp', ref_vcf,
                                        '-O', os.path.join(outdir,my_id) + '/' + my_id +'.vcf'],
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE,
                                       universal_newlines=True)
            stderr, stdout = process.communicate()
            loggi(tmp_log, tmp_err, stdout, stderr, 'a')
            print('done')
            rc = process.returncode
            if (rc == 0):
                with open(all_log, "a") as log:
                    log.write('gatk HaplotypeCaller done with ' + my_id + '\n')
                    log.write('processing done with ' + my_id + '\n')
            else:
                with open(all_log, "a") as log:
                    log.write('gatk HaplotypeCaller failed with ' + my_id + '\n')
                if (os.path.exists(os.path.join(outdir, my_id) + '/' + my_id + '.vcf')):
                    os.remove(os.path.join(outdir, my_id) + '/' + my_id + '.vcf')
                sys.exit(10)
        except Exception:
            with open(all_log, "a") as log:
                log.write('gatk HaplotypeCaller failed with ' + my_id + '\n')
            if (os.path.exists(os.path.join(outdir, my_id) + '/' + my_id + '.vcf')):
                os.remove(os.path.join(outdir, my_id) + '/' + my_id + '.vcf')


    if (os.path.exists(os.path.join(outdir,my_id) + '/' + my_id + '_final.bam')):
        os.remove(os.path.join(outdir,my_id) + '/' + my_id + '_final.bam')



def process_bam_trimmed(my_id):

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
    print('gatk haplotypecaller')
    with open(all_log, "a") as log:
        log.write('gatk haplotypecaller '+ my_id +'\n')
    if(os.path.exists(os.path.join(outdir,my_id) + '/' + my_id +'.vcf')):
        print('already exists, go further')
        with open(all_log, "a") as log:
            log.write(os.path.join(outdir,my_id) + '/' + my_id +'.vcf'+' already exists, complete'+ '\n')
    else:
        try:
            process = subprocess.Popen(['gatk', 'HaplotypeCaller',
                                        '--java-options', javapars,
                                        '-R', processed_ref,
                                        '-I', indir + my_id + '.bam',
                                        '--dbsnp', ref_vcf,
                                        '-O', os.path.join(outdir,my_id) + '/' + my_id +'.vcf'],
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE,
                                       universal_newlines=True)
            stderr, stdout = process.communicate()
            loggi(tmp_log, tmp_err, stdout, stderr, 'a')
            print('done')
            rc = process.returncode
            if (rc == 0):
                with open(all_log, "a") as log:
                    log.write('gatk HaplotypeCaller done with ' + my_id + '\n')
                    log.write('processing done with ' + my_id + '\n')
            else:
                with open(all_log, "a") as log:
                    log.write('gatk HaplotypeCaller failed with ' + my_id + '\n')
                if (os.path.exists(os.path.join(outdir, my_id) + '/' + my_id + '.vcf')):
                    os.remove(os.path.join(outdir, my_id) + '/' + my_id + '.vcf')
                rm(my_id)
                sys.exit(10)
        except Exception:
            with open(all_log, "a") as log:
                log.write('gatk HaplotypeCaller failed with ' + my_id + '\n')
            if (os.path.exists(os.path.join(outdir, my_id) + '/' + my_id + '.vcf')):
                os.remove(os.path.join(outdir, my_id) + '/' + my_id + '.vcf')
            rm(my_id)
            sys.exit(10)





def stats(my_id,clause ):
    if (clause == '_ready.bam' ):
        strf = os.path.join(outdir,my_id) + '/' + my_id +'_ready.bam'
        statfilecov = os.path.join(final_outdir,my_id) + '/'  + 'stats_nodup_cov.txt'
        statfile = os.path.join(final_outdir,my_id) + '/'  +'stats_nodup_num.txt'

    elif (clause == '_final.bam' ):
        strf = os.path.join(outdir,my_id) + '/' + my_id +'_final.bam'
        statfilecov = os.path.join(final_outdir,my_id) + '/'  + 'stats_final_cov.txt'
        statfile = os.path.join(final_outdir,my_id) + '/'  +'stats_final_num.txt'

    elif (clause == '.bam'):
        strf = indir + my_id + '.bam'
        statfilecov = os.path.join(final_outdir,my_id) + '/' + 'stats_start_cov.txt'
        statfile = os.path.join(final_outdir,my_id) + '/' + 'stats_start_num.txt'
    elif (clause == '.vcf'): #vcf
        strf = os.path.join(final_outdir,my_id) + '/' + my_id +'.vcf'
        statfile = os.path.join(final_outdir,my_id) + '/' +'vcf_stats.txt'

    ##########start#########################################

    if (clause == '.vcf'): #vcf_number of peaks
        file_rs = open(strf, "r")
        line_rs = file_rs.readline()
        n = 0
        while line_rs.startswith("##"):
            n += 1
            line_rs = file_rs.readline()
        file_rs.close()

        vcf = pd.read_csv(strf,
                             sep='\t', skiprows=n)

        with open(statfile, "w") as g:
            g.write('number of peaks '+ str(vcf.shape[0]))


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
               os.path.join(outdir, my_id) + '/' + my_id + '_chop.bam',
               os.path.join(outdir, my_id) + '/' + my_id + '_sorted.bam',
               os.path.join(outdir, my_id) + '/' + my_id + '_formatted.bam',
               os.path.join(outdir, my_id) + '/' + my_id + '_metrics.txt',
               os.path.join(outdir, my_id) + '/' + my_id + '_ready.bam',
               os.path.join(outdir, my_id) + '/' + my_id + '.table',
               os.path.join(outdir, my_id) + '/' + my_id + '_final.bai',
               os.path.join(outdir, my_id) + '/' + my_id + '_final.bam',
               os.path.join(outdir, my_id) + '/' + my_id + '.vcf',
               os.path.join(outdir, my_id) + '/' + my_id + '.vcf.idx'
               ]
    for i in rm_list:
        if(os.path.exists(i) and i!=os.path.join(final_outdir, my_id) + '/' + my_id + '.vcf'
                and i!=os.path.join(final_outdir, my_id) + '/' + my_id + '.vcf.idx'):
            os.remove(i)
    try:
        os.rmdir(os.path.join(outdir, my_id))
        print("Directory '% s' has been removed successfully" % os.path.join(outdir, my_id))
    except OSError as error:
        print(error)
        print("Directory '% s' can not be removed" % os.path.join(outdir, my_id))


def cp(my_id):
    if (os.path.exists(os.path.join(outdir, my_id) + '/' + my_id + '.vcf')
       and (os.path.join(outdir, my_id) + '/' + my_id + '.vcf' != os.path.join(final_outdir, my_id) + '/' + my_id + '.vcf')):
        shutil.copy(os.path.join(outdir, my_id) + '/' + my_id + '.vcf',
                    os.path.join(final_outdir, my_id) + '/' + my_id + '.vcf')
    elif (os.path.join(outdir, my_id) + '/' + my_id + '.vcf' == os.path.join(final_outdir, my_id) + '/' + my_id + '.vcf'):
        print(my_id + '.vcf is already in the directory')
    else:
        print(my_id + '.vcf does not exist to copy')

    if (os.path.exists(os.path.join(outdir, my_id) + '/' + my_id + '.vcf.idx')
       and (os.path.join(outdir, my_id) + '/' + my_id + '.vcf.idx' == os.path.join(final_outdir, my_id) + '/' + my_id + '.vcf.idx')):
        shutil.copy(os.path.join(outdir, my_id) + '/' + my_id + '.vcf.idx',
                    os.path.join(final_outdir, my_id) + '/' + my_id + '.vcf.idx')



########my_id = 'BAM00030'#######################
def pipe_my_id(my_id):
    print('start SNP '+ my_id)
    #my_id = 'BAM00030'
    tmp_path = os.path.join(outdir,my_id)
    #делаем директории
    if (os.path.exists(tmp_path)):
        pass
    else:
        os.mkdir(tmp_path)
    tmp_path = os.path.join(final_outdir,my_id)
    if (os.path.exists(tmp_path)):
        pass
    else:
        os.mkdir(tmp_path)
    #делаем обработку
    if(os.path.exists(os.path.join(final_outdir,my_id) + '/' + my_id + '.vcf')):
        print(my_id + ' vcf file from gatk already exists in final directory, skip processing')
        stats(my_id, '.vcf')
        rm(my_id)
    elif (os.path.exists(os.path.join(outdir, my_id) + '/' + my_id + '.vcf') ):
        print(my_id+ ' vcf file from gatk already exists, skip processing, start to copy to final directory')
        cp(my_id)
        stats(my_id, '.vcf')
        rm(my_id)
    else:
        print("Beginning processing with " + my_id)

        process_bam_trimmed(my_id)
        cp(my_id)
        stats(my_id, '.vcf')
        rm(my_id)
    #print("Beginnig with " + my_id)
    #process_bam(my_id)
    #stats(my_id,'.bam')
    #stats(my_id,'_final.bam')
    #stats(my_id,'_ready.bam')




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
outdir = config["Directories"]["temp_data_out"]
print(outdir)
final_outdir = config["Directories"]["final_data_out"]
print(final_outdir)
logdir = os.path.join(maindir,config["Directories"]["data_log"])
javapars = config["Parameters"]["javaparameters"]
met = os.path.join(maindir,config["Files"]["metadata"])
picard = config["Soft"]["picard"]
# pipe ----------------------------------------------------
my_id = sys.argv[2]

pipe_my_id(my_id)
