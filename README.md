# viadastra steps:

**1)** python step1_mk_links.py /media/ElissarDisk/ADASTRA/parameters/CONFIG.cfg

**2)** python step2_reference_processing.py /media/ElissarDisk/ADASTRA/parameters/CONFIG.cfg ; bash tgsender_no_logs.sh ; echo end
 

**snp calling for one file:**

python step2_snp_calling.py /media/ElissarDisk/ADASTRA/parameters/CONFIG.cfg BAM00030

python step2_snp_calling.py /media/ElissarDisk/ADASTRA/parameters/CONFIG.cfg BAM00030 ; bash tgsender_no_logs.sh ; echo end

**snp calling for lots files:**

python step2_snp_all.py 8 /media/ElissarDisk/ADASTRA/parameters/CONFIG.cfg

python step2_snp_all.py 8 /media/ElissarDisk/ADASTRA/parameters/CONFIG.cfg ; bash tgsender_no_logs.sh ; echo end

**3) vcf_filtrate:**


python step3_filtrate_vcf_get_stats.py 5 20 6 /media/ElissarDisk/ADASTRA/parameters/CONFIG.cfg ; bash tgsender_no_logs.sh ; echo end

python step3_fastfiltrate_vcf_get_stats.py 5 /media/ElissarDisk/ADASTRA/parameters/CONFIG.cfg ; bash tgsender_no_logs.sh ; echo end


python AS_figure_ref_bias.py /media/ElissarDisk/ADASTRA/processed_data/ref_alt_count.tsv /media/ElissarDisk/ADASTRA/processed_data


**4) making bad maps**


python step4_pullby_bads_and_sort.py /media/ElissarDisk/ADASTRA/parameters/CONFIG.cfg ; bash tgsender_no_logs.sh ; echo end

python step4_to_babachi.py  /media/ElissarDisk/ADASTRA/parameters/CONFIG.cfg ; bash tgsender_no_logs.sh ; echo end

babachi /media/ElissarDisk/ADASTRA/processed_data/pulled_chipseq_tobabachi.tsv --visualize

babachi /media/ElissarDisk/ADASTRA/processed_data/pulled_atacseq_tobabachi.tsv --visualize

python step5_group_by_bad.py 5 /media/ElissarDisk/ADASTRA/parameters/CONFIG.cfg ; bash tgsender_no_logs.sh ; echo end

