#!/usr/bin/env nextflow

// using DSL-2
nextflow.enable.dsl=2

// params set in nextflow.config

// import modules
include { get_irods_bam       } from '../modules/get_irods_bam.nf'
include { get_local_bam       } from '../modules/get_local_bam.nf'
include { samtools_index      } from '../modules/samtools_index.nf'
include { bedtools_bamtofastq } from '../modules/bedtools_bamtofastq.nf'
include { MOSDEPTH            } from './modules/nf-core/mosdepth/main'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

// get gene coordinates
process get_gene_coords {
  tag "${gene}"
  label 'normal'
  
  input:
  val(gene)

  output:
  tuple val(gene), path("${gene}.bed")

  script:
  def chr = params.no_chr ? "" : "chr"
  def host = params.genome_build == "GRCh37" ? "host = 'grch37.ensembl.org'," : ""
  """
  #!/usr/bin/env Rscript

  # libraries
  library(biomaRt)
  
  # use the ensembl database
  ensembl <- useEnsembl(biomart = "genes", $host
                        dataset = "hsapiens_gene_ensembl")

  # get coordinates
  gene_coordinates <- getBM(
    attributes = c("hgnc_symbol", "chromosome_name", "start_position",
                   "end_position", "strand"),
    filters = "hgnc_symbol",
    values = c("${gene}"),
    mart = ensembl) |>
    # filter out non-canonical chromosomes
    dplyr::filter(!grepl("HSCHR|PATCH|CTG", chromosome_name)) |>
    dplyr::transmute(chr = paste0("$chr", chromosome_name),
                     start = start_position, end = end_position,
                     gene = hgnc_symbol)

  # save
  readr::write_tsv(gene_coordinates, "${gene}.bed", col_names = FALSE)
  """
}

// get gene coverage
process get_gene_cov {
  tag "${meta.id}_${gene}"
  label 'normal10gb'
  
  input:
  tuple val(meta), path(per_base_bed), val(gene), path(gene_bed)
  
  output:
  tuple val(meta), val(gene), path("${meta.id}_${gene}_gene_cov.tsv")
  
  script:
  """
  module load bedtools2-2.29.0/python-3.10.10
  echo -e "chr\\tstart\\tend\\tcov\\tgene\\tid" > ${meta.id}_${gene}_gene_cov.tsv
  zcat $per_base_bed |
  bedtools intersect -a stdin -b $gene_bed -wa -wb |
  cut -f1-4,8 |
  awk -v OFS='\\t' '{print \$0, "${meta.id}"}' \\
  >> ${meta.id}_${gene}_gene_cov.tsv
  """
}

// concat gene coverage
process concat_gene_cov {
  tag "${meta.id}"
  label 'normal10gb'
  publishDir "${params.out_dir}/${meta.donor_id}/${meta.id}/genotyping/",
    mode: "copy"
  
  input:
  tuple val(meta), val(gene), path(gene_covs)
  
  output:
  tuple val(meta), path("${meta.id}_gene_cov.tsv")

  script:
  """
  head -1 ${gene_covs[0]} > ${meta.id}_gene_cov.tsv
  for file in ${gene_covs} ; do
    sed 1d \$file >> ${meta.id}_gene_cov.tsv
  done
  """
}

// annotate mutations
process annotate_mutations {
  tag "${meta.id}"
  label 'normal10gb'
    
  input:
  tuple val(meta), path(mutations)
  
  output:
  tuple val(meta), path("${meta.id}_annotated_mutations.tsv")
  
  script:
  """
  # annotate mutations
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
    --chr ${chr} \\
    --mutations ${mutations} \\
    --bam ${bam} \\
    --min_bq ${params.min_bq} \\
    --mask ${params.mask} \\
    --min_mq ${params.min_mq}
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
  label 'normal50gb'
  publishDir "${params.out_dir}/${meta.donor_id}/",
    mode: "copy"
  
  input:
  path rmd
  tuple val(meta), val(ids), path(geno), path(gene_covs), 
        path(summary_txt), path(global_txt), path(regions_txt)
  
  output:
  path "${meta.donor_id}_report.html"
  
  script:
  def c_ids = 'c("' + ids.join('", "') + '")'
  def c_gene_covs = 'c("' + gene_covs.join('", "') + '")'
  def c_geno = 'c("' + geno.join('", "') + '")'
  def c_summary_txt = 'c("' + summary_txt.join('", "') + '")'
  def c_global_txt = 'c("' + global_txt.join('", "') + '")'
  def c_regions_txt = 'c("' + regions_txt.join('", "') + '")'
  """
  #!/usr/bin/env Rscript

  # capture params
  my_params <-
    list(
      donor_id = "${meta.donor_id}",
      gene_covs = ${c_gene_covs},
      ids = ${c_ids},
      geno = ${c_geno},
      summary_txt = ${c_summary_txt},
      global_txt = ${c_global_txt},
      regions_txt = ${c_regions_txt})
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

  // validate input parameters
  validateParameters()

  // print summary of supplied parameters
  log.info paramsSummaryLog(workflow)

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

  // get fastq
  if ( params.bamtofastq ) {
    bedtools_bamtofastq(ch_bam2)
  }

  // initialize fasta file with meta map
  fasta = params.fasta ? Channel.fromPath(params.fasta).map{ it -> [ [id:it.baseName], it ] }.collect() : Channel.empty()

  // get bait set
  bait_set = Channel.fromPath(params.bait_set, checkIfExists: true)

  // get genome-wide coverage
  MOSDEPTH(samtools_index.out.combine(bait_set), fasta)

  // get gene coords
  Channel.of(params.genes.split(','))
    | set { ch_gene }
  get_gene_coords(ch_gene)

  // get gene coverage
  get_gene_cov(MOSDEPTH.out.per_base_bed.combine(get_gene_coords.out))
  concat_gene_cov(get_gene_cov.out.groupTuple())

  // create value channel of chromosomes
  Channel.of(1..22).concat(Channel.of('X', 'Y')) | set { chromosomes }
  if (!params.no_chr) {
    chromosomes = chromosomes | map { 'chr' + it }
  }

  // annotate, genotype and collect mutations
  annotate_mutations(ch_mutations)
  genotype_mutations(samtools_index.out.join(annotate_mutations.out).combine(chromosomes))
  concat_mutations(genotype_mutations.out.groupTuple(size: 24))

  // generate report
  concat_mutations.out \
  .join(concat_gene_cov.out)
  .join(MOSDEPTH.out.summary_txt) \
  .join(MOSDEPTH.out.global_txt) \
  .join(MOSDEPTH.out.regions_txt) \
  | map { meta, geno, gene_covs, summary_txt, global_txt, regions_txt ->
          [meta.subMap(['donor_id']),
            meta.id, geno, gene_covs, summary_txt, global_txt, regions_txt]
  }
  | groupTuple()
  | set { report_input }
  report(rmd, report_input)
  
}