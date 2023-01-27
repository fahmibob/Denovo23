#!/bin/bash
preFetch () {
  echo "Retriving file from $1"
  echo "File is saved at $2"
  prefetch $1 -o $2 -X $3
}

fastq-dump () {
  parallel-fastq-dump --sra-id $1 \
  --outdir $2 --threads $3 \
  --split-3 \
  --gzip
}

qc () {
  fastqc $1 --outdir $2 -threads $3
}

multi () {
  multiqc --file-list $1 -outdir $2
}

eSearch (){
  esearch -db sra -query $1 | efetch -format runinfo | tail -n+2 >> $2
}

fastp_single () {
  fastp -i $1 -o $2 -h $3 -j $4 \
  -q 20 -l 36 -r -x -w $5
}

fastp_paired () {
  fastp -i $1 -I $2 -o $3 -O $4 \
  -h $5 -j $6 \
  -q 20 -l 36 -r -x -w $7
}

pairSortRna () {
  if [ "$5" = nil ]; then
    rm -rf /home/fahmi/sortmerna/run/kvdb/
    sortmerna --ref ~/sortmerna/smr_v4.3_default_db.fasta \
    --reads $1 --reads $2 \
    --paired_in --fastx --other $3 -out2 \
    --threads $4
    cp /home/fahmi/sortmerna/run/out/aligned.log $6
  else
    rm -rf /home/fahmi/sortmerna/$5/kvdb/
    sortmerna --ref ~/sortmerna/smr_v4.3_default_db.fasta \
    --reads $1 --reads $2 \
    --paired_in --fastx --other $3 -out2 \
    --threads $4 --workdir /home/fahmi/sortmerna/$5
    cp /home/fahmi/sortmerna/$5/out/aligned.log $6
  fi

}

singleSortRna () {
  if [ "$4" = nil ]; then
    rm -rf /home/fahmi/sortmerna/run/kvdb/
    sortmerna --ref ~/sortmerna/smr_v4.3_default_db.fasta --reads $1 \
    --fastx --other $2\
    --threads $3
    cp /home/fahmi/sortmerna/run/out/aligned.log $5
  else
    rm -rf /home/fahmi/sortmerna/$4/kvdb/
    sortmerna --ref ~/sortmerna/smr_v4.3_default_db.fasta --reads $1 \
    --fastx --other $2\
    --threads $3 --workdir /home/fahmi/sortmerna/$4
    cp /home/fahmi/sortmerna/$4/out/aligned.log $5
  fi
}


trimmomatic_single () {
  phred=$(zcat $1 |\
    head -n 1000 | \
    awk '{if(NR%4==0) printf("%s",$0);}' |  od -A n -t u1 | \
    awk 'BEGIN{min=100;max=0;} \
        {for(i=1;i<=NF;i++) \
            {if($i>max) max=$i; \
                 if($i<min) min=$i;}}END \
            {if(max<=74 && min<59) \
                       print "-phred33"; \
             else \
             if(max>73 && min>=64) \
                       print "-phred64"; \
             else \
             if(min>=59 && min<64 && max>73) \
                       print "-phred64"; else print "Unknown score encoding!";}')

  java -jar trimmomatic-0.39.jar SE -threads $3 $phred \
  -trimlog $4 $1 $2 \
  ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 \
  TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
}


trimmomatic_double () {
  phred=$(zcat $1 |\
    head -n 1000 | \
    awk '{if(NR%4==0) printf("%s",$0);}' |  od -A n -t u1 | \
    awk 'BEGIN{min=100;max=0;} \
        {for(i=1;i<=NF;i++) \
            {if($i>max) max=$i; \
                 if($i<min) min=$i;}}END \
            {if(max<=74 && min<59) \
                       print "-phred33"; \
             else \
             if(max>73 && min>=64) \
                       print "-phred64"; \
             else \
             if(min>=59 && min<64 && max>73) \
                       print "-phred64"; else print "Unknown score encoding!";}')

  java -jar trimmomatic-0.39.jar PE -threads $7 $phred \
  -trimlog $8 $1 $2 $3 $4 $5 $6 \
  ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 \
  TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
}

# Check if the function exists (bash specific)
if declare -f "$1" > /dev/null
then
  # call arguments verbatim
  "$@"
else
  # Show a helpful error
  echo "'$1' is not a known function name" >&2
  exit 1
fi
