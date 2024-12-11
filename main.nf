#!/usr/bin/env nextflow

// using DSL-2
nextflow.enable.dsl=2

// params set in nextflow.config

// import modules
include { get_irods_bam  } from '../modules/get_irods_bam.nf'
include { get_local_bam  } from '../modules/get_local_bam.nf'
include { samtools_index } from '../modules/samtools_index.nf'
include { MOSDEPTH       } from './modules/nf-core/mosdepth/main'   

// annotate mutations
process annotate_mutations {
  tag "${meta.id}"
  label 'normal10gb'
  publishDir "${params.out_dir}/${meta.donor_id}/${meta.id}/genotyping/",
    mode: "copy"
    
  input:
  tuple val(meta), path(mutations)
  
  output:
  tuple val(meta), path("${meta.id}_annotated_mutations.tsv")
  
  script:
  """
  annotate_mutations.R \\
    --mutations ${mutations}
  mv annotated_mutations.tsv ${meta.id}_annotated_mutations.tsv
  """
}

// genotype mutations
process genotype_mutations {
  tag "${meta.id}_${chr}"
  label 'normal10gb'
  
  input:
  tuple val(meta), path(bam), path(bai), path(mutations), val(chr)
  
  output:
  tuple val(meta), path("${meta.id}_${chr}_genotyped_mutations.tsv")
        
  script:
  """
  genotype_mutations.R \\
    --donor_id ${meta.donor_id} \\
    --chr ${chr} \\
    --mutations ${mutations} \\
    --bam ${bam} \\
    --q ${params.q} \\
    --mask ${params.mask} \\
    --mq ${params.mq}
  mv genotyped_mutations.tsv ${meta.id}_${chr}_genotyped_mutations.tsv
  """
}

// concat mutations
process concat_mutations {
  tag "${meta.id}"
  label 'normal10gb'
  publishDir "${params.out_dir}/${meta.donor_id}/${meta.id}/genotyping/",
    mode: "copy"
  
  input:
  tuple val(meta), path(geno)
  
  output:
  tuple val(meta), path("${meta.id}_genotyped_mutations.tsv")

  script:
  """
  head -1 ${geno[0]} > ${meta.id}_genotyped_mutations.tsv
  for file in ${geno} ; do
    sed 1d \$file >> ${meta.id}_genotyped_mutations.tsv
  done
  """
}

// report
process report {
  tag "${meta.donor_id}"
  label 'normal10gb'
  publishDir "${params.out_dir}/${meta.donor_id}/",
    mode: "copy"
  
  input:
  path rmd
  tuple val(meta), val(ids), path(geno),
        path(summary_txt), path(global_txt), path(regions_txt)
  
  output:
  path "${meta.donor_id}_report.html"
  
  script:
  def ids_c = 'c("' + ids.join('", "') + '")'
  def geno_c = 'c("' + geno.join('", "') + '")'
  def summary_txt_c = 'c("' + summary_txt.join('", "') + '")'
  def global_txt_c = 'c("' + global_txt.join('", "') + '")'
  def regions_txt_c = 'c("' + regions_txt.join('", "') + '")'
  """
  #!/usr/bin/env Rscript

  # capture params
  my_params <-
    list(
      donor_id = "${meta.donor_id}",
      ids = ${ids_c},
      geno = ${geno_c},
      summary_txt = ${summary_txt_c},
      global_txt = ${global_txt_c},
      regions_txt = ${regions_txt_c})
  saveRDS(my_params, "params.rds")

  # render
  rmarkdown::render(
    "${rmd}",
    output_dir = "./",
    output_file = "${meta.donor_id}_report.html",
    params = my_params)
  """
}

workflow {

  // get input bams
  Channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true, sep: "\t")
    | map { row ->
            def meta = [donor_id: row.donor_id, id: row.id]
            [meta, file(row.bam)]
    }
    | set { ch_bam }

  // get input mutations
  Channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true, sep: "\t")
    | map { row ->
            def meta = [donor_id: row.donor_id, id: row.id]
            [meta, file(row.mutations, checkIfExists: true)]
    }
    | set { ch_mutations }

  // get rmd file
  rmd = file("${baseDir}/bin/report.Rmd")

  // stage bams
  if (params.location == "irods") {
    get_irods_bam(ch_bam)
    | set { ch_bam2 }
  } else {
    get_local_bam(ch_bam)
    | set { ch_bam2 }
  }

  // index bams
  samtools_index(ch_bam2)

  // initialize fasta file with meta map:
  fasta = params.fasta ? Channel.fromPath(params.fasta).map{ it -> [ [id:it.baseName], it ] }.collect() : Channel.empty()

  // get bait set
  bait_set = Channel.fromPath(params.bait_set, checkIfExists: true)

  // mosdepth
  MOSDEPTH(samtools_index.out.combine(bait_set), fasta)

  // create value channel of chromosomes
  Channel.of(1..22).concat(Channel.of('X', 'Y')) | set { chromosomes }
  if (!params.no_chr) {
    chromosomes = chromosomes | map { 'chr' + it }
  }

  // annotate, genotype and collect mutations
  annotate_mutations(ch_mutations)
  genotype_mutations(samtools_index.out.join(ch_mutations).combine(chromosomes))
  concat_mutations(genotype_mutations.out.groupTuple(size: 24))

  // generate report
  concat_mutations.out \
  .join(MOSDEPTH.out.summary_txt) \
  .join(MOSDEPTH.out.global_txt) \
  .join(MOSDEPTH.out.regions_txt) \
  | map { meta, geno, summary_txt, global_txt, regions_txt ->
          [meta.subMap(['donor_id']),
            meta.id, geno, summary_txt, global_txt, regions_txt]
  }
  | groupTuple()
  | set { report_input }
  report(rmd, report_input)
  
}