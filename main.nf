#!/usr/bin/env nextflow

// using DSL-2
nextflow.enable.dsl=2

// params set in nextflow.config

// import modules
include { get_irods_bam        } from './modules/local/get_irods_bam'
include { get_local_bam        } from './modules/local/get_local_bam'
include { samtools_index       } from './modules/local/samtools_index'
include { concat_gene_cov      } from './modules/local/concat_gene_cov'
include { annotate_mutations   } from './modules/local/annotate_mutations'
include { bedtools_intersect_snps } from './modules/local/bedtools_intersect_snps'
include { genotype_mutations; genotype_mutations as genotype_snps } from './modules/local/genotype_mutations'
include { concat_mutations; concat_mutations as concat_snps } from './modules/local/concat_mutations'
include { concat_snps_per_cell } from './modules/local/concat_snps_per_cell'
include { generate_nr_nv       } from './modules/local/generate_nr_nv'
include { plot_baf             } from './modules/local/plot_baf'
include { knit_plots as knit_baf; knit_plots as knit_vdj } from './modules/local/knit_plots'
include { merge_plots as merge_baf; merge_plots as merge_vdj } from './modules/local/merge_plots'
include { report               } from './modules/local/report'
include { MOSDEPTH; MOSDEPTH as MOSDEPTH_VDJ } from './modules/nf-core/mosdepth/main'
include { plot_vdj_cov         } from './modules/local/plot_vdj_cov'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

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
            [meta, "snps", row.snps]
    }
    | filter { it[2] != "NA" }
    | map { meta, set, snps -> [meta, set, file(snps, checkIfExists: true)] }
    | set { ch_snps }

  // get refcds file
  refcds = file(params.refcds, checkIfExists: true)

  // get report rmd file
  report_rmd = file("${baseDir}/bin/report.Rmd")

  // get parent rmd file
  parent_rmd = file("${baseDir}/bin/parent.Rmd")

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

  // knit VDJ coverage plots
  knit_vdj(plot_vdj_cov.out)
  knit_vdj.out
    | map { meta, vdj_plots, rmds  ->
            return [meta.subMap(['donor_id']), vdj_plots, rmds]
    }
    | groupTuple()
    | map { meta, vdj_plots, rmds ->
            def flat_vdj_plots = vdj_plots.flatten()
            return [meta, flat_vdj_plots, rmds]
    }
    | set { merge_vdj_input }
  merge_vdj(parent_rmd, merge_vdj_input, "VDJ coverage", "vdj_cov")

  // genotype mutations
  genotype_mutations(samtools_index.out.join(ch_mutations))

  // generate mutation summaries
  genotype_mutations.out
    | map { meta, set, geno -> [meta.subMap('donor_id'), set, geno] }
    | groupTuple(by: [0, 1])
    | set { ch_all_genos }
  concat_mutations(ch_all_genos)
  annotate_mutations(concat_mutations.out, refcds)

  // if seq_type = dnahyb, subset snps to those in the panel
  bait_set_hyb2 = file(params.bait_set_hyb, checkIfExists: true)
  if (params.seq_type == "dnahyb") {
    ch_snps2 = bedtools_intersect_snps(ch_snps, bait_set_hyb2)
  } else {
    ch_snps2 = ch_snps
  }

  // genotype SNPs in chunks of 100,000
  ch_snps_split = ch_snps2.splitText(by: 100000, file: true, keepHeader: true)
  ch_bams_x_snps = samtools_index.out.combine(ch_snps_split, by: 0)
  genotype_snps(ch_bams_x_snps)
  concat_snps_per_cell(genotype_snps.out.groupTuple(by: [0, 1]))
  concat_snps_per_cell.out
    | map { meta, set, geno -> [meta.subMap('donor_id'), set, geno] }
    | groupTuple(by: [0, 1])
    | set { ch_all_snps }
  concat_snps(ch_all_snps)

  // plot BAF from genotyped SNPs
  plot_baf(concat_snps_per_cell.out)

  // knit and merge BAF plots
  knit_baf(plot_baf.out)
  knit_baf.out
    | map { meta, baf_plots, rmds  ->
            return [meta.subMap(['donor_id']), baf_plots, rmds]
    }
    | groupTuple()
    | map { meta, baf_plots, rmds ->
            def flat_baf_plots = baf_plots.flatten()
            return [meta, flat_baf_plots, rmds]
    }
    | set { merge_baf_input }
  merge_baf(parent_rmd, merge_baf_input, "BAF plots", "genotyping/snps")
  
  // generate report
  if (params.knit_report) {
    plot_baf.out \
      .join(plot_vdj_cov.out) \
      .join(MOSDEPTH.out.summary_txt) \
      .join(MOSDEPTH.out.global_txt) \
      .join(MOSDEPTH.out.regions_txt) \
      | map { meta, baf_plots, vdj_plots, summary_txt, global_txt, regions_txt ->
              [meta.subMap(['donor_id']),
                meta.id, baf_plots, vdj_plots, summary_txt, global_txt, regions_txt]
      }
      | groupTuple()
      | set { report_input }
    report(report_rmd, report_input)
  }
  
}