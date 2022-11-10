import pandas as pd
import os
import sys
import configparser
import pysam
import subprocess32 as subprocess
import pathlib
import argparse
from argparse import RawTextHelpFormatter
from steps.step1_mk_links import *
from steps.step2_snp_all import *
from steps.step2_filter import *
import requests


class SmartFormatter(argparse.HelpFormatter):

    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)


def tgmessage(message,chatID, apiToken):
    apiURL = f'https://api.telegram.org/bot{apiToken}/sendMessage'

    try:
        response = requests.post(apiURL, json={'chat_id': chatID, 'text': message})
        print(response.text)
    except Exception as e:
        print(e)

parser = argparse.ArgumentParser(description='viadastra pipeline, perfoms snp-calling, filtrating and works with babachi and mixalime',
                                 formatter_class=SmartFormatter)

parser.add_argument("-j","--jobs", help="Number of jobs", default='4', required=False)
parser.add_argument("-m","--memfree", help="Memfree parameter for gnu parallel", default='40G', required=False)
parser.add_argument("-tt",'--apitoken', help= "Telegram bot api token", default='1992203014:AAGXCU5ta31M-R10axejbBtxRJd0L1PNOow')
parser.add_argument("-ti",'--chatid', help= "Telegram bot chat id", default='639261746')
parser.add_argument("-a", '--allele_reads_tr',
                    help= "Allelic reads threshold. Input SNPs will be filtered by ref_read_count >= x and alt_read_count >= x",
                    default=5,required=False)

requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument("-c","--config", help="Path for cfg file, which contains all parameters", required=True)
helpm = "R|Define a processing step please.\n1 creates soft links and directories\n2 for reference processing"\
        + "\n3 snp calling for one file\n4 snp calling for all files "
requiredNamed.add_argument("-s","--step",choices=['1', '2', '3', '4', '5','6'],help=helpm, required=True)
#parser.parse_args(['-h'])
args = parser.parse_args()

jobs = args.jobs
memfree = args.memfree
path = args.config
step = args.step
chatid = args.chatid
apitoken = args.apitoken

print(step)
print(path)

if step=='1':
    print('starting to create soft links and directories')
    mklinks(path)

elif (step=='2'):
    print('reference processing')
    print('it is not done yet')
elif (step=='3'):
    print('snp calling for one file')
    print('it is not done yet')
    tgmessage("it is not done yet",chatid,apitoken)
elif (step=='4'):
    print('snp calling for all files')
    ret = call_all(jobs,path,memfree)
    print(ret)
    tgmessage("snp calling was finished",chatid,apitoken)
elif(step=='5'):
    print('filtrate,pull,sort,babachi')
    threshold=args.allele_reads_tr
    filter(path,threshold,jobs)
    pullsort(path)
    print("done successfully")

elif(step=='6'):
    print('annotate,group by bads')
elif(step=='7'):
    print('mixalime')





