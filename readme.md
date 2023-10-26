# Project summary

DNA methylation QTL (methQTL) mapping of the sex chromosome of a chicken population based on 124 individuals. DNA methylation was measured with MeDIP-seq. Method for producing the samples can be found [here](https://doi.org/10.1038/s41559-020-01310-1). The code supplied will calculate and map QTL, filter based lod-score threshold (based on permutation tests), remove outliers and produce a table (`data/mQTL_table1.0.xlsx`, `data/f8_mqtl_table_chrz.csv`) and boxplots (`img/mQTL_plots1.0.pdf`)


```
project
|-- bin
|   |-- methtable.R
|   |-- outlier_check.R
|   `-- qtlscan.R
|-- data
|   |-- f8_mqtl_table_chrz.csv
|   |-- f8_mqtl_table_chrz_unfiltered.Rds
|   |-- f8_pmap_galgal6.csv
|   `-- mQTL_table1.0.xlsx
|-- img
|   `-- mQTL_plots1.0.pdf
|-- readme.md
`-- tmp
    |-- methQTL_toRemove.txt
    `-- tmp.pdf
```