# CCLE_DataParser

Data parser for CCLE

     * filename  source/link
     * baseline expression: CCLE_expression.csv, source link: https://depmap.org/portal/download/ (19Q3: log2(RPKM+1))
     * baseline expression: CCLE_expression.csv, source link: https://depmap.org/portal/download/ (20Q1: log2(TPM+1))
     * mutation: CCLE_mutations.csv, source link: https://depmap.org/portal/download/
     * copy number variation: CCLE_gene_cn.csv, source link: https://depmap.org/portal/download/
     * cell line: sample_info.csv, source link: https://depmap.org/portal/download/
     * compound: , source link:
     * drug sensitivity: CCLE_NP24.2009_Drug_data_2015.02.24.csv, source link:
     * model list: https://cellmodelpassports.sanger.ac.uk/downloads

# How to use it?

Download files listed above and put it in a folder (e.g., /data/DR/db/CCLE/):

```{python}
DB_PATH = '/data/DR/db/CCLE/'
DB_FILE = {'MODEL':'model_list_20200204.csv',
           'EXP_tpm':'CCLE_expression.csv', # log2(TPM+1)
           'EXP_rpkm': 'CCLE_expression_log2RPKM.csv', # log2(RPKM+1)
           'CNV':'CCLE_gene_cn.csv',
           'MUT':'CCLE_mutations.csv',
           'RESP':'CCLE_NP24.2009_Drug_data_2015.02.24.csv',
           'CELL': 'sample_info.csv',
           'DRUG':''}
```

type $python useCCLE.py to generate the following files in a tidy format ready for analysis
