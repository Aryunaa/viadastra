import configparser
import os
#import sys


def createConfig(path):

    config = configparser.ConfigParser()
    config.add_section("Directories")
    config.set("Directories", "maindir", "/media/ElissarDisk/ADASTRA/")
    config.set("Directories", "bam", "BAM/")
    config.set("Directories", "data_in", "data/")
    config.set("Directories", "data_out", "processed_data/")
    config.set("Directories", "data_log", "logs/data_processing/")
    config.set("Directories", "ref_log", "logs/reference_processing/")

    config.add_section("Files")
    config.set("Files", "metadata", "parameters/metadata.tsv")
    config.set("Files", "ref_in", "reference/hg19.fa")
    config.set("Files", "ref_vcf", "reference/00-common_all.vcf.gz")
    config.set("Files", "ref_out1", "processed_ref/genome-norm.fasta")
    config.set("Files", "ref_out2", "processed_ref/genome-norm.dict")
    config.set("Files", "exception_list", "logs/exceptions_list")
    config.set("Files", "processing_list", "logs/processing_list")

    config.add_section("Parameters")
    config.set("Parameters", "JavaParameters", "-Xmx12G -XX:ParallelGCThreads=4")

    with open(path, "w") as config_file:
        config.write(config_file)


def readConfig_SNP(path):
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
    metadata = config.get("Files","metadata")
    processing_list = config.get("Files","processing_list")
    javapars = config.get("Parameters", "JavaParameters")

    confdict = {'maindir': maindir, 'indir': os.path.join(maindir,indir),
                'processed_ref': os.path.join(maindir,processed_ref),
                'outdir': os.path.join(maindir,data_out),'ref_vcf': os.path.join(maindir,ref_vcf),
                'logdir': os.path.join(maindir,data_log), 'metadata':os.path.join(maindir,metadata),
                'processing_list' : os.path.join(maindir,processing_list),
                'javapars': javapars}
    return confdict


def readConfig_ref(path):
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
    indir = config.get("Files", "ref_in")
    logdir = config.get("Directories", "ref_log")
    out1 = config.get("Files", "ref_out1")
    out2 = config.get("Files", "ref_out2")

    confdict = {'maindir': maindir, 'indir': os.path.join(maindir,indir),
                'logdir': os.path.join(maindir, logdir),
                'out1': os.path.join(maindir,out1),'out2': os.path.join(maindir,out2),
                'loglog':os.path.join(maindir,logdir,'stdout'), 'logerr':os.path.join(maindir,logdir,'stderr')}

    return confdict

'''
def readConfig(path, dir):
    config = configparser.ConfigParser()
    config.read(path)




    config.set("Files", "ref_out2", "processed_ref/genome-norm.dict")
    config.set("Files", "exception_list", "logs/exception_list")
    config.set("Files", "processing_list", "logs/processing_list")

    config.set("Parameters", "JavaParameters", "-Xmx12G -XX:ParallelGCThreads=4")

    maindir = config.get("Directories", "maindir")
    bam = config.get("Directories", "bam")
    data_in = config.get("Directories", "data_in")
    data_out = config.get("Directories","data_out")
    data_log = config.get("Directories", "data_log")
    ref_log = config.get("Directories", "ref_log")

    ref_in = config.get("Files", "ref_in")
    ref_vcf = config.get("Files", "ref_vcf")
    ref_out1 = config.get("Files", "ref_out1")
    ref_out2 = config.get 
    metadata = config.get("Files", "metadata")
'''


#path = "C:/Users/Aryuna/Desktop/IB/viadastra_pretty/config.cfg"
#path = sys.argv[1]
#createConfig(path)

"""
maindir = '/media/ElissarDisk/ADASTRA/'
indir = maindir + 'data/'
processed_ref = maindir + 'processed_ref/genome-norm.fasta'
ref_vcf = maindir + 'reference/00-common_all.vcf.gz'

outdir = maindir + 'processed_data/'
logdir = maindir + 'logs/data_processing/'
javapars = '-Xmx12G -XX:ParallelGCThreads=4'
"""
'''
dicti = readConfig_SNP(path)
maindir =dicti['maindir']
print(maindir)
indir = dicti['indir']
print(indir)
processed_ref = dicti['processed_ref']
print(processed_ref)
ref_vcf = dicti['ref_vcf']
print(ref_vcf)
outdir = dicti['outdir']
print(outdir)
logdir = dicti['logdir']
print(logdir)
javapars = dicti['javapars']
print(javapars)
'''