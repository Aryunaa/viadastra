import pandas as pd
import os
import pysam
import sys
import pybedtools
from pybedtools import BedTool
import subprocess32 as subprocess

processed_data = '/home/ariuna/rafaello/bedfiles/babachi'
metapath = '/home/ariuna/rafaello/viadastra/additional/metadata_yes_no.tsv'
metadata = pd.read_csv(metapath, sep = '\t')
fit_path_atacseq = '/home/ariuna/rafaello/nbfit_try/grouped_bads/atacseq'
fit_path_chipseq = '/home/ariuna/rafaello/nbfit_try/grouped_bads/chipseq'

scorefiles = '/home/ariuna/rafaello/bedfiles/scorefiles'

listi = os.listdir(bedfiles)