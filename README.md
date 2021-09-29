# viadastra steps:

**1)** python step1_mk_links.py /media/ElissarDisk/ADASTRA/parameters/CONFIG.cfg

**2)** python step2_reference_processing.py /media/ElissarDisk/ADASTRA/parameters/CONFIG.cfg ; bash tgsender_no_logs.sh ; echo end

**3)** 

**snp calling for one file:**

python step2_snp_calling.py /media/ElissarDisk/ADASTRA/parameters/CONFIG.cfg BAM00030

python step2_snp_calling.py /media/ElissarDisk/ADASTRA/parameters/CONFIG.cfg BAM00030 ; bash tgsender_no_logs.sh ; echo end

**snp calling for lots files:**

python step2_snp_all.py 4 /media/ElissarDisk/ADASTRA/parameters/CONFIG.cfg

python step2_snp_all.py 4 /media/ElissarDisk/ADASTRA/parameters/CONFIG.cfg ; bash tgsender_no_logs.sh ; echo end

**vcf_filtrate:**

python step2_filtrate_vcf_get_stats.py 5 /media/ElissarDisk/ADASTRA/parameters/CONFIG.cfg ; bash tgsender_no_logs.sh ; echo end
