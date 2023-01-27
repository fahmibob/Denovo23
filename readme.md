# Denovo assembly for SRA Illumina data

#1 retrieve sra data
a. download the summary file of target species group in sra
b. filter the preferred species and biggest data using r script
c. get the sra run code using esearch by using the expriment code
d. get the sra run data by using sra code using prefetch

#2. preprocess data
a. run quality control using preprocess.nf (fastp and sortmerna)

#3. perform denovo assembly
