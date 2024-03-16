#viadastra
Tutorial is here: 
python3 main.py -h 


1. create config file and meta\
meta must contain columns:
`['BAM', 'ID','path']`
2. manage env
```php
conda env create -f viadastra/viadastra_env.yml
```

2. make soft links 
```php
python3 main.py -c viadastra/configs/config.cfg -s 1
```
2. snp calling for all files
```php
python3 main.py -c viadastra/configs/config.cfg -s 4
```
or 
```php
parallel -j 4 python viadastra/steps/step2_snp_calling.py viadastra/configs/config.cfg :::: logs/processing_list ;  curl -s -X POST https://api.telegram.org/bot1992203014:AAGXCU5ta31M-R10axejbBtxRJd0L1PNOow/sendMessage -d chat_id=639261746 -d text="snp calling выполнен"
```

3. filtrate pull and sort 


