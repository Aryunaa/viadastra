import subprocess32 as subprocess
import sys
import os
import configparser


path = sys.argv[1]

config = configparser.ConfigParser()
config.read(path)
maindir = config["Directories"]["maindir"]
print(maindir)
indir = os.path.join(maindir,config["Files"]["ref_in"])
print(indir)
logdir = os.path.join(maindir,config["Directories"]["ref_log"])
print(logdir)
out1 = os.path.join(maindir,config["Files"]["ref_out1"])
print(out1)
out2 = os.path.join(maindir, config["Files"]["ref_out2"])
print(out2)
loglog = os.path.join(logdir,'stdout')
print(loglog)
logerr = os.path.join(logdir,'stderr')
print(logerr)
ref_vcf = os.path.join(maindir,config["Files"]["ref_vcf"])

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
process = subprocess.Popen(['gatk', 'IndexFeatureFile', '-I',ref_vcf],
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


