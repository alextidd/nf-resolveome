#!/usr/bin/env nextflow

// using DSL-2
nextflow.enable.dsl=2

// params set in nextflow.config

// import modules
include { get_irods_bam       } from './modules/local/get_irods_bam'
include { get_local_bam       } from './modules/local/get_local_bam'
include { samtools_index      } from './modules/local/samtools_index'
include { get_gene_coords     } from './modules/local/get_gene_coords'
include { get_gene_cov        } from './modules/local/get_gene_cov'
include { concat_gene_cov     } from './modules/local/concat_gene_cov'
include { annotate_mutations  } from './modules/local/annotate_mutations'
include { genotype_mutations; genotype_mutations as genotype_snps } from './modules/local/genotype_mutations'
include { concat_snps         } from './modules/local/concat_snps'
include { concat_mutations    } from './modules/local/concat_mutations'
include { generate_nr_nv      } from './modules/local/generate_nr_nv'
include { plot_heatmap        } from './modules/local/plot_heatmap'
include { plot_baf            } from './modules/local/plot_baf'
include { report              } from './modules/local/report'
include { MOSDEPTH            } from './modules/nf-core/mosdepth/main'
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
            [meta, "snps", file(row.snps, checkIfExists: true)]
    }
    | set { ch_snps }

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

  // genotype mutations
  genotype_mutations(samtools_index.out.join(ch_mutations))

  // generate mutation summaries
  genotype_mutations.out
    | map { meta, set, geno -> [meta.subMap('donor_id'), set, geno] }
    | groupTuple(by: [0, 1])
    | set { ch_all_genos }
  concat_mutations(ch_all_genos)
  annotate_mutations(concat_mutations.out)

  // // plot heatmap from genotyped mutations
  // generate_nr_nv(ch_all_genos)
  // plot_heatmap(generate_nr_nv.out.out)

  // genotype SNPs in chunks of 100,000
  ch_snps_split = ch_snps.splitText(by: 100000, file: true, keepHeader: true)
  ch_bams_x_snps = samtools_index.out.combine(ch_snps_split, by: 0)
  genotype_snps(ch_bams_x_snps)
  concat_snps(genotype_snps.out.groupTuple())

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