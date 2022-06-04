import subprocess32 as subprocess
import os
import configparser
import sys as sys
import pathlib
import argparse


parser = argparse.ArgumentParser(description='script that launches snp calling')

parser.add_argument("-j","--jobs", help="Number of jobs", default=4, required=False)
parser.add_argument("-m","--memfree", help="Memfree parameter for gnu parallel", default='40G', required=False)
parser.add_argument("-c","--config", help="Path for cfg file", required=True)
args = parser.parse_args()

jobs = args.jobs
path = args.config
memfree = args.memfree

# reading config --------------------------
config = configparser.ConfigParser()
config.read(path)
maindir = config["Directories"]["maindir"]
print(maindir)
indir = os.path.join(maindir,config["Directories"]["data_in"])
print(indir)
processed_ref = os.path.join(maindir,config["Files"]["ref_out1"])
ref_vcf = os.path.join(maindir,config["Files"]["ref_vcf"])
outdir = os.path.join(maindir,config["Directories"]["temp_data_out"])
logdir = os.path.join(maindir,config["Directories"]["data_log"])
mainlogs = os.path.join(maindir,config["Directories"]["mainlogs"])
javapars = config["Parameters"]["javaparameters"]
met = os.path.join(maindir,config["Files"]["metadata"])
processing_list = os.path.join(maindir, config["Files"]["processing_list"])

dir = pathlib.Path(__file__).parent.absolute()
script = os.path.join(dir,'step2_snp_calling.py')

all_log = os.path.join(os.path.join(maindir, mainlogs), 'whole_log')
with open(all_log, "w") as log:
    log.write('STARTING! all'+ '\n')
try:
    process = subprocess.Popen(['parallel', '--memfree',memfree,'--retry-failed','--joblog',
                                os.path.join(os.path.join(maindir, mainlogs), 'parallel_log'),'-j', jobs,'python', script ,path ,'::::',processing_list],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   universal_newlines=True
                                   )
    stderr, stdout = process.communicate()
    with open(all_log, "a") as log:
        log.write(stdout)
    with open(all_log, "a") as err:
        err.write(stderr)
    with open(all_log, "a") as log:
        log.write("script has been performed successfully")

except KeyboardInterrupt:
    with open(all_log, "a") as log:
        log.write('exception KeyboardInterrupt emerged'+ '\n')
    sys.exit(10)
except Exception:
    with open(all_log, "a") as log:
        log.write(str(Exception) + '\n')
        log.write('Some exception emerged'+ '\n')
    sys.exit(10)
