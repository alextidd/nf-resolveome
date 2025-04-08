#!/usr/bin/env nextflow

// using DSL-2
nextflow.enable.dsl=2

// params set in nextflow.config

// import modules
include { get_irods_bam        } from './modules/local/get_irods_bam'
include { get_local_bam        } from './modules/local/get_local_bam'
include { samtools_index       } from './modules/local/samtools_index'
include { get_gene_coords      } from './modules/local/get_gene_coords'
include { get_gene_cov         } from './modules/local/get_gene_cov'
include { concat_gene_cov      } from './modules/local/concat_gene_cov'
include { annotate_mutations   } from './modules/local/annotate_mutations'
include { genotype_mutations; genotype_mutations as genotype_snps } from './modules/local/genotype_mutations'
include { concat_mutations; concat_mutations as concat_snps  } from './modules/local/concat_mutations'
include { concat_snps_per_cell } from './modules/local/concat_snps_per_cell'
include { generate_nr_nv       } from './modules/local/generate_nr_nv'
include { plot_baf             } from './modules/local/plot_baf'
include { report               } from './modules/local/report'
include { MOSDEPTH; MOSDEPTH as MOSDEPTH_VDJ } from './modules/nf-core/mosdepth/main'
include { plot_vdj_cov         } from './modules/local/plot_vdj_cov'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
include { concat_snps } from './modules/local/concat_snps/main.nf'

workflow {

  // validate input parameters
  validateParameters()

  // print summary of supplied parameters
  log.info paramsSummaryLog(workflow)

  // get input bams
  Channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true)
    | map { row ->
            def meta = [donor_id: row.donor_id, id: row.id]
            [meta, file(row.bam)]
    }
    | set { ch_bam }

  // get input mutations
  Channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true)
    | map { row ->
            def meta = [donor_id: row.donor_id, id: row.id]
            [meta, "mutations", file(row.mutations, checkIfExists: true)]
    }
    | set { ch_mutations }
  
  // get input SNPs
  Channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true)
    | map { row ->
            def meta = [donor_id: row.donor_id, id: row.id]
            [meta, "snps", file(row.snps, checkIfExists: true)]
    }
    | set { ch_snps }
  
  // get refcds file
  refcds = file(params.refcds, checkIfExists: true)

  // get rmd file
  rmd = file("${baseDir}/bin/report.Rmd")

  // initialize fasta file with meta map
  fasta = params.fasta ? Channel.fromPath(params.fasta).map{ it -> [ [id:it.baseName], it ] }.collect() : Channel.empty()

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

  // get genome-wide and immune panel coverage
  bait_set_hyb = Channel.fromPath(params.bait_set_hyb, checkIfExists: true)
  MOSDEPTH(samtools_index.out.combine(bait_set_hyb), fasta)

  // get VDJ coverage
  bait_set_vdj = Channel.fromPath(params.bait_set_vdj, checkIfExists: true)
  MOSDEPTH_VDJ(samtools_index.out.combine(bait_set_vdj), fasta)

  // plot VDJ coverage
  bait_set_vdj2 = file(params.bait_set_vdj, checkIfExists: true)
  plot_vdj_cov(MOSDEPTH_VDJ.out.regions_bed, bait_set_vdj2)

  // get gene coords
  ch_gene = Channel.fromPath(params.genes, checkIfExists: true).splitText().map{it -> it.trim()}
  get_gene_coords(ch_gene)

  // get gene coverage
  get_gene_cov(MOSDEPTH.out.per_base_bed.combine(get_gene_coords.out))
  concat_gene_cov(get_gene_cov.out.groupTuple())

  // genotype mutations
  genotype_mutations(samtools_index.out.join(ch_mutations))

  // generate mutation summaries
  genotype_mutations.out
    | map { meta, set, geno -> [meta.subMap('donor_id'), set, geno] }
    | groupTuple(by: [0, 1])
    | set { ch_all_genos }
  concat_mutations(ch_all_genos)
  annotate_mutations(concat_mutations.out, refcds)

  // genotype SNPs in chunks of 100,000
  ch_snps_split = ch_snps.splitText(by: 100000, file: true, keepHeader: true)
  ch_bams_x_snps = samtools_index.out.combine(ch_snps_split, by: 0)
  genotype_snps(ch_bams_x_snps)
  concat_snps_per_cell(genotype_snps.out.groupTuple(by: [0, 1]))
  concat_snps_per_cell.out
    | map { meta, set, geno -> [meta.subMap('donor_id'), set, geno] }
    | groupTuple(by: [0, 1])
    | set { ch_all_snps }
  concat_snps(ch_all_snps)

  // plot BAF from genotyped SNPs
  plot_baf(concat_snps.out)

  // generate report
  plot_baf.out.out \
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