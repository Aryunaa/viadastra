import pandas as pd
import os
import sys
import configparser
import pysam
import subprocess32 as subprocess
import pathlib
import argparse
from argparse import RawTextHelpFormatter
#from steps.step1_mk_links import *



class SmartFormatter(argparse.HelpFormatter):

    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)

parser = argparse.ArgumentParser(description='viadastra pipeline, perfoms snp-calling, filtrating and works with babachi and mixalime',
                                 formatter_class=SmartFormatter)
parser.add_argument("-c","--config", help="Path for cfg file, which contains all parameters", required=True)
helpm = "R|Define a processing step please.\n1 creates soft links and directories\n2.1 for reference processing"\
        + "\n2.2.1 snp calling for one file\n2.2 snp calling for all files "
print(helpm)
parser.add_argument("-s","--step",choices=['1', '2.1', '2.2.1', '2.2', 'e'],help=helpm, required=True)
parser.add_argument("-j","--jobs", help="Number of jobs", default='4', required=False)
parser.add_argument("-m","--memfree", help="Memfree parameter for gnu parallel", default='40G', required=False)
args = parser.parse_args()
jobs = args.jobs
path = args.config
memfree = args.memfree

