import subprocess32 as subprocess
import os
import pandas as pd

print('start process')

metadata_path = '/home/ariuna/rafaello/viadastra/additional/metadata_v5.tsv'
#processing_list_path = '/home/ariuna/rafaello/viadastra/additional/processing_list'
processed_data = '/home/ariuna/rafaello/bedfiles'

tmp_path = os.path.join(processed_data,'babachi')
if (os.path.exists(tmp_path)):
    pass
else:
    os.mkdir(tmp_path)
metadata = pd.read_csv(metadata_path,sep='\t')






pull_chip = metadata[metadata['BADgroup']=='chipseq']

intersect = list(pull_chip['ID'])
paths_rs = []
for my_id in intersect:
    if(os.path.exists(os.path.join(processed_data, my_id+'.snps.bed'))):
        paths_rs.append(os.path.join(processed_data, my_id+'.snps.bed'))

print(paths_rs)

header_list = ['#CHROM', 'POS1','POS2', 'ID', 'REF', 'ALT', 'REF_COUNT', 'ALT_COUNT']
pulled_chips_rs = pd.concat([pd.read_csv(f,sep='\t',names=header_list) for f in paths_rs])
pulled_chips_rs.to_csv(os.path.join(tmp_path,'ppulled_chipseq.tsv'),mode='w', header=False,index=False,sep='\t')
print('chipseq pulled')
process = subprocess.run(['bedtools', 'sort','-i',
                            os.path.join(tmp_path,'ppulled_chipseq.tsv'),'>',os.path.join(tmp_path,'pulled_chipseq.tsv')],
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           universal_newlines=True
                           )
print('chipseq sorted')
print(pulled_chips_rs.shape[0])
#_______________________________________________________________________________
pull_atac = metadata[metadata['BADgroup']=='atacseq']

intersect = list(pull_atac['ID'])
paths_rs = []
for my_id in intersect:
    if(os.path.exists(os.path.join(processed_data, my_id+'.snps.bed'))):
        paths_rs.append(os.path.join(processed_data, my_id+'.snps.bed'))
print(paths_rs)
pulled_atac_rs = pd.concat([pd.read_csv(f,sep='\t',names=header_list) for f in paths_rs])
pulled_atac_rs.to_csv(os.path.join(tmp_path,'ppulled_atacseq.tsv'),mode='w', header=False,index=False,sep='\t')
print('atacseq pulled')
process = subprocess.run(['bedtools', 'sort','-i',
                            os.path.join(tmp_path,'ppulled_atacseq.tsv'),'>',os.path.join(tmp_path,'pulled_atacseq.tsv')],
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           universal_newlines=True
                           )
print('atacseq sorted')
print(pulled_atac_rs.shape[0])

#______________________
pull_oht = metadata[ (metadata['Treatment']=='OHT') & (metadata['BADgroup']=='atacseq')]
pull_nooht = metadata[ (metadata['Treatment']!='OHT') & (metadata['BADgroup']=='atacseq')]

intersect_oht = list(pull_oht['ID'])
paths_rs_oht = []
for my_id in intersect_oht:
    if(os.path.exists(os.path.join(processed_data, my_id+'.snps.bed'))):
        paths_rs_oht.append(os.path.join(processed_data, my_id+'.snps.bed'))
print(paths_rs_oht)
pulled_atac_rs_oht = pd.concat([pd.read_csv(f,sep='\t',names=header_list) for f in paths_rs_oht])
pulled_atac_rs_oht.to_csv(os.path.join(tmp_path,'ppulled_atacseq_oht.tsv'),mode='w', header=False,index=False,sep='\t')
print('atacseqoht pulled')
process = subprocess.run(['bedtools', 'sort','-i',
                            os.path.join(tmp_path,'ppulled_atacseq_oht.tsv'),'>',os.path.join(tmp_path,'pulled_atacseq_oht.tsv')],
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           universal_newlines=True
                           )
print('atacseqoht sorted')
print(pulled_atac_rs_oht.shape[0])


intersect_nooht = list(pull_nooht['ID'])
paths_rs_nooht = []
for my_id in intersect_nooht:
    if(os.path.exists(os.path.join(processed_data, my_id+'.snps.bed'))):
        paths_rs_nooht.append(os.path.join(processed_data, my_id+'.snps.bed'))
print(paths_rs_nooht)
pulled_atac_rs_nooht = pd.concat([pd.read_csv(f,sep='\t',names=header_list) for f in paths_rs_nooht])
pulled_atac_rs_nooht.to_csv(os.path.join(tmp_path,'ppulled_atacseq_nooht.tsv'),mode='w', header=False,index=False,sep='\t')
print('atacseqnooht pulled')
process = subprocess.run(['bedtools', 'sort','-i',
                            os.path.join(tmp_path,'ppulled_atacseq_nooht.tsv'),'>',os.path.join(tmp_path,'pulled_atacseq_nooht.tsv')],
                           universal_newlines=True
                           )
print('atacseqnooht sorted')
print(pulled_atac_rs_oht.shape[0])




#________________________________________

'''

process = subprocess.run(['babachi', os.path.join(tmp_path, 'pulled_chipseq.tsv'),
                              '-O', tmp_path])
process = subprocess.run(['babachi', os.path.join(tmp_path, 'pulled_atacseq.tsv'),
                              '-O', tmp_path])

#babachi visualize pulled_chipseq.tsv -b pulled_chipseq.badmap.bed
process = subprocess.run(['babachi','visualize', os.path.join(tmp_path, 'pulled_chipseq.tsv'),
                              '-b', os.path.join(tmp_path, 'pulled_chipseq.badmap.bed')
                              ])

process = subprocess.run(['babachi','visualize', os.path.join(tmp_path, 'pulled_atacseq.tsv'),
                              '-b', os.path.join(tmp_path, 'pulled_atacseq.badmap.bed')
                              ])
'''

print('done')