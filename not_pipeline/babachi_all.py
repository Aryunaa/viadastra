import subprocess32 as subprocess
import os
import configparser
import sys as sys
import pathlib
import pandas as pd
import shlex

babachidir = '/home/ariuna/tcga_atacseq/babachi_all/lihc'

def babachi(configpath):
    config = configparser.ConfigParser()
    config.read(configpath)
    maindir = config["Directories"]["maindir"]
    print(maindir)
    filt_bed_rs = os.path.join(maindir, config["Directories"]["rssnps"])
    mainlogs = os.path.join(maindir, config["Directories"]["mainlogs"])
    all_log = os.path.join(mainlogs, 'babachilogs')
    met = '../metas/atacseqmeta_lihc.tsv'
    metadata = pd.read_csv(met, sep='\t')
    tmp_path = babachidir
    with open(all_log, "w") as log:
        log.write('babachi' + '\n')

    process = subprocess.Popen(
        shlex.split('find -type f -name "*.vcf" | parallel -j ' + jobs + ' babachi filter -a ' + trs, ' -O ', filt_bed),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True)
    stderr, stdout = process.communicate()
    with open(all_log, "a") as log:
        log.write(stdout)
        log.write(stderr)
    myids = list(metadata['ID'])



    print('babachi done')