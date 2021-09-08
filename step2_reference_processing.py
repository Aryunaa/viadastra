import subprocess32 as subprocess
import os


'''
maindir = '/media/ElissarDisk/ADASTRA/'
outdir = maindir + 'processed_ref/'
logdir = maindir + 'logs/reference_processing/'

indir = maindir + 'reference/hg19.fa'
out1 = outdir +'genome-norm.fasta'
out2 = outdir + 'genome-norm.dict'
loglog = logdir + "stdout"
logerr = logdir + "stderr"
'''


def readConfig_ref(path, dir):
    config = configparser.ConfigParser()
    config.read(path)
    maindir = config.get("Directories", "maindir")
    logdir = config.get("Directories", "ref_log")
    indir = config.get("Files", "ref_in")
    out1 = config.get("Files", "ref_out1")
    out2 = config.get("Files", "ref_out2")

    if (dir == 'maindir'):
        return (maindir)
    elif (dir == 'logdir'):
        return (maindir + logdir)
    elif (dir == 'indir'):
        return (maindir + indir)
    elif (dir == 'out1'):
        return (maindir + out1)
    elif (dir == 'out2'):
        return (maindir + out2)
    elif (dir == 'loglog'):
        return (maindir + logdir + 'stdout')
    elif (dir == 'logerr'):
        return (maindir + logdir + 'stderr')

readConfig_ref(path,dir)
maindir =readConfig_ref(path,'maindir')
print(maindir)
indir = readConfig_ref(path,'indir')
print(indir)
logdir = readConfig_ref(path,'logdir')
out1 = readConfig_ref(path,'out1')
out2 = readConfig_ref(path,'out2')
loglog = readConfig_ref(path,'stdout')
logerr = readConfig_ref(path,'stderr')



#___________
process = subprocess.Popen(['picard', 'NormalizeFasta', 'I='+indir, 'O=' + out1],
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE,
                     universal_newlines=True)
stderr, stdout = process.communicate()
stdout.split('\n')
stderr.split('\n')


with open(loglog, "w") as log:
    log.write(stdout)
with open(logerr, "w") as err:
    err.write(stderr)

#----------
process = subprocess.Popen(['samtools', 'faidx', out1],
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE,
                     universal_newlines=True)
stderr, stdout = process.communicate()
stdout.split('\n')
stderr.split('\n')

with open(loglog, "a") as log:
    log.write(stdout)
with open(logerr, "a") as err:
    err.write(stderr)
#___________
process = subprocess.Popen(['picard', 'CreateSequenceDictionary', 'R='+out1, 'O='+out2],
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE,
                     universal_newlines=True)
stderr, stdout = process.communicate()
stdout.split('\n')
stderr.split('\n')

with open(loglog, "a") as log:
    log.write(stdout)
with open(logerr, "a") as err:
    err.write(stderr)
#___________

