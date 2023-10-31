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
3. filtrate pull and sort 


