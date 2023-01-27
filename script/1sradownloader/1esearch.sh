#!/bin/bash
export NCBI_API_KEY=3427cb5832a349c2516ada0aafcccb3b4009
s=$1 #Input experiment
esearch -db sra -query $s | efetch -format runinfo | tail -n+2 >> out.csv
