// enable dsl2
nextflow.enable.dsl=2
//specify output in parameter job and common option
params.query=false
params.sortme=false
params.forking=7
params.jobs="default"
params.filterbyFile=false
params.ref="$baseDir/tools/sortmerna/smr_v4.3_default_db.fasta"

workflow {
  if (params.query){
    Channel.fromFilePairs(params.query, flat:true)
    .set{pair_ch}

    if (params.filterbyFile){
      def filterList = new File(params.filterbyFile).collect {it}
      pair_ch=pair_ch.map { if (it[0].toString() in filterList){it}}
    }

    pair_ch=pair_ch.map {data -> tuple(data[1].parent.getName(), data[0], data[1], data[2] )}

    fastp(pair_ch)
    sortmerna(fastp.out)

  } else if (params.sortme) {
    Channel.fromFilePairs(params.sortme, flat:true)
    .set{sortme_ch}

    if (params.filterbyFile){
      def filterList = new File(params.filterbyFile).collect {it}
      sortme_ch=sortme_ch.map { if (it[0].toString() in filterList){it}}
    }

    sortme_ch=sortme_ch.map {data -> tuple(data[1].parent.getName(), data[0], data[1], data[2] )}

    sortmerna(sortme_ch)
  }
}

process fastp {
  publishDir "$baseDir/output/${params.jobs}/fastp/${speciesName}", mode: 'copy'

  errorStrategy 'ignore'
  maxForks 4

  input:
  tuple val(speciesName), val(sampleName), path(forward), path(reverse)
  output:
  tuple val(speciesName), val(sampleName), path("*_1.fastq.gz"), path("*_2.fastq.gz")

  script:
  """
  fastp -i $forward -I $reverse -o 'FP'$sampleName'_1.fastq.gz' \
  -O 'FP'$sampleName'_2.fastq.gz' -q 20 -l 36 -r -x -w 2
  """
}

process sortmerna {
  publishDir "$baseDir/output/${params.jobs}/sortmerna/${speciesName}", mode: 'move'

  errorStrategy 'ignore'
  maxForks 8

  input:
  tuple val(speciesName), val(sampleName), path(forward), path(reverse)

  output:
  tuple path("SMR${sampleName}_fwd*"), path("SMR${sampleName}_rev*")
  path("out/*.log")

  script:
  sortme="$baseDir/tools/sortmerna/sortmerna"
  """
  $sortme --ref $params.ref \
  --reads $forward --reads $reverse --paired_in --fastx --other SMR$sampleName \
  -out2 --threads 8 --workdir .
  """
}
