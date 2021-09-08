import configparser
import os

def createConfig(path):
    """
    Create a config file


    processed_ref = maindir + 'processed_ref/genome-norm.fasta'
    ref_vcf = maindir + 'reference/00-common_all.vcf.gz'


    """
    config = configparser.ConfigParser()
    config.add_section("Directories")
    config.set("Directories", "maindir", "/media/ElissarDisk/ADASTRA/")
    config.set("Directories", "bam", "BAM/")
    config.set("Directories", "data_in", "data/")
    config.set("Directories", "data_out", "processed_data/")
    config.set("Directories", "data_log","logs/data_processing/")
    config.set("Directories", "ref_log", "logs/reference_processing/")

    config.add_section("Files")
    config.set("Files", "metadata", "parameters/metadata.tsv")
    config.set("Files", "ref_in", "reference/hg19.fa")
    config.set("Files", "ref_vcf", "reference/00-common_all.vcf.gz")
    config.set("Files", "ref_out1", "processed_ref/genome-norm.fasta")
    config.set("Files", "ref_out2", "processed_ref/genome-norm.dict")


    config.add_section("Parameters")
    config.set("Parameters", "JavaParameters", "-Xmx12G -XX:ParallelGCThreads=4")

    with open(path, "w") as config_file:
        config.write(config_file)


def readConfig_SNP(path,dir):
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
    data_out = config.get("Directories", "data_out")
    data_log = config.get("Directories", "data_log")

    javapars = config.get("Parameters", "JavaParameters")


    if (dir == 'maindir'):
        return(maindir)
    elif (dir == 'indir'):
        return (maindir+indir)
    elif (dir == 'processed_ref'):
        return (maindir+processed_ref)
    elif (dir == 'ref_vcf'):
        return (maindir+ref_vcf)
    elif (dir == 'outdir'):
        return(maindir+data_out)
    elif (dir == 'logdir'):
        return(maindir+data_log)
    elif (dir == 'javapars'):
        return(javapars)

def readConfig_ref(path,dir):
    config = configparser.ConfigParser()
    config.read(path)

    """
    maindir = '/media/ElissarDisk/ADASTRA/'
    
    logdir = maindir + 'logs/reference_processing/'
    
    indir = maindir + 'reference/hg19.fa'
    out1 = outdir +'genome-norm.fasta'
    out2 = outdir + 'genome-norm.dict'
    loglog = logdir + "stdout"
    logerr = logdir + "stderr"
    """
    maindir = config.get("Directories", "maindir")
    logdir = config.get("Directories", "ref_log")
    indir = config.get("Files", "ref_in")
    out1 = config.get("Files", "ref_out1")
    out2 = config.get("Files", "ref_out2")

    if (dir == 'maindir'):
        return(maindir)
    elif (dir == 'logdir'):
        return(maindir+logdir)
    elif (dir == 'indir'):
        return(maindir+indir)
    elif (dir == 'out1'):
        return(maindir+out1)
    elif (dir == 'out2'):
        return(maindir+out2)
    elif (dir == 'loglog'):
        return(maindir+logdir+'stdout')
    elif (dir == 'logerr'):
        return(maindir+logdir+'stderr')

path = "C:/Users/Aryuna/Desktop/IB/viadastra/config.cfg"
createConfig(path)

"""
maindir = '/media/ElissarDisk/ADASTRA/'
indir = maindir + 'data/'
processed_ref = maindir + 'processed_ref/genome-norm.fasta'
ref_vcf = maindir + 'reference/00-common_all.vcf.gz'

outdir = maindir + 'processed_data/'
logdir = maindir + 'logs/data_processing/'
javapars = '-Xmx12G -XX:ParallelGCThreads=4'
"""
maindir = readConfig_SNP(path,'maindir')
print(maindir)
indir = readConfig_SNP(path,'indir')
print(indir)
processed_ref = readConfig_SNP(path,'processed_ref')
print(processed_ref)
ref_vcf = readConfig_SNP(path,'ref_vcf')
print(ref_vcf)
outdir = readConfig_SNP(path,'outdir')
print(outdir)
logdir = readConfig_SNP(path,'logdir')
print(logdir)
javapars = readConfig_SNP(path,'javapars')
print(javapars)