#!/usr/bin/env nextflow

// params set in nextflow.config

// import modules
include { get_irods_bam        } from './modules/local/get_irods_bam'
include { get_local_bam        } from './modules/local/get_local_bam'
include { samtools_index       } from './modules/local/samtools_index'
include { concat_gene_cov      } from './modules/local/concat_gene_cov'
include { annotate_mutations   } from './modules/local/annotate_mutations'
include { bedtools_intersect_snps } from './modules/local/bedtools_intersect_snps'
include { bedtools_intersect_vdj_bed } from './modules/local/bedtools_intersect_vdj_bed'
include { genotype_mutations; genotype_mutations as genotype_snps } from './modules/local/genotype_mutations'
include { concat_mutations; concat_mutations as concat_snps } from './modules/local/concat_mutations'
include { concat_snps_per_cell } from './modules/local/concat_snps_per_cell'
include { generate_nr_nv       } from './modules/local/generate_nr_nv'
include { plot_baf             } from './modules/local/plot_baf'
include { merge_pdf as merge_pdf_vdj; merge_pdf as merge_pdf_baf } from './modules/local/merge_pdf'
include { knit_qc_report       } from './modules/local/knit_qc_report'
include { MOSDEPTH; MOSDEPTH as MOSDEPTH_VDJ } from './modules/nf-core/mosdepth/main'
include { plot_vdj_cov         } from './modules/local/plot_vdj_cov'

workflow {

  // get input bams
  channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true)
    | map { row ->
            def meta = [donor_id: row.donor_id, id: row.id]
            [meta, file(row.bam, checkIfExists: params.location == "local")]
    }
    | set { ch_bam }

  // get input mutations
  if (params.run_mutations) {
    channel
      .fromPath(params.samplesheet)
      .splitCsv(header: true)
      | map { row ->
              def meta = [donor_id: row.donor_id, id: row.id]
              [meta, "mutations", file(row.mutations, checkIfExists: true)] }
      | set { ch_mutations }
  }
  
  // get input SNPs
  if (params.run_snps) {
    channel
      .fromPath(params.samplesheet)
      .splitCsv(header: true)
      | map { row ->
              def meta = [donor_id: row.donor_id, id: row.id]
              [meta, "snps", row.snps]
      }
      | filter { it[2] != "NA" }
      | map { meta, set, snps -> [meta, set, file(snps, checkIfExists: true)] }
      | set { ch_snps }
  }

  // get refcds file
  refcds = file(params.refcds, checkIfExists: true)

  // initialise fasta file with meta map
  fasta = params.fasta ? channel.fromPath(params.fasta).map{ it -> [ [id:it.baseName], it ] }.collect() : channel.empty()

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

  // get hyb panel coverage
  dnahyb_bed = channel.fromPath(params.dnahyb_bed, checkIfExists: true)
  dnahyb_bed2 = file(params.dnahyb_bed, checkIfExists: true)
  MOSDEPTH(samtools_index.out.combine(dnahyb_bed), fasta)

  // get VDJ coverage
  // if seq_type = dnahyb, restrict the VDJ bed to regions captured by the hyb panel
  vdj_bed = channel.fromPath(params.vdj_bed, checkIfExists: true)
  if (params.seq_type == "dnahyb") {
    vdj_bed = bedtools_intersect_vdj_bed(vdj_bed, dnahyb_bed2)
  }
  vdj_bed = vdj_bed.first()
  MOSDEPTH_VDJ(samtools_index.out.combine(vdj_bed), fasta)

  // plot VDJ coverage
  plot_vdj_cov(MOSDEPTH_VDJ.out.regions_bed, vdj_bed)

  // merge VDJ coverage plots into PDF
  merge_pdf_vdj(
    plot_vdj_cov.out
      .map { meta, plots -> [meta.subMap(['donor_id']), plots] }
      .groupTuple()
      .map { meta, plots -> [meta, plots.flatten()] },
    "vdj_coverage_plots",
    "vdj_cov"
  )

  // genotype mutations
  if (params.run_mutations) {
    genotype_mutations(samtools_index.out.join(ch_mutations))
    genotype_mutations.out
      | map { meta, set, geno -> [meta.subMap('donor_id'), set, geno] }
      | groupTuple(by: [0, 1])
      | set { ch_all_genos }
    concat_mutations(ch_all_genos)
    annotate_mutations(concat_mutations.out, refcds, params.no_chr_prefix)
  }

  // genotype SNPs and plot BAF
  if (params.run_snps) {

    // if seq_type = dnahyb, subset SNPs to those in the panel
    if (params.seq_type == "dnahyb") {
      ch_snps2 = bedtools_intersect_snps(ch_snps, dnahyb_bed2)
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
    plot_baf(concat_snps_per_cell.out, refcds)

    // merge BAF plots into PDF
    merge_pdf_baf(
      plot_baf.out
        .map { meta, plots -> [meta.subMap(['donor_id']), plots] }
        .groupTuple()
        .map { meta, plots -> [meta, plots.flatten()] },
      "baf_plots",
      "genotyping/snps"
    )

  }

  // knit QC report
  // (coverage metrics are always available; SNP-derived metrics are only
  // included when run_snps is set)
  if (params.knit_qc_report) {

    qc_report_rmd = file("${baseDir}/bin/qc_report.Rmd", checkIfExists: true)

    // well/plate per id, e.g. well "A1", "B10" (only needed for the QC report)
    channel
      .fromPath(params.samplesheet)
      .splitCsv(header: true)
      | map { row -> [[donor_id: row.donor_id, id: row.id], row.well, row.plate] }
      | set { ch_wells }

    ch_qc_mosdepth =
      MOSDEPTH.out.summary_txt
        | join(MOSDEPTH.out.global_txt)
        | join(MOSDEPTH.out.regions_txt)
        | map { meta, summary_txt, global_txt, regions_txt ->
                [[donor_id: meta.donor_id, id: meta.id], summary_txt, global_txt, regions_txt] }
        | join(ch_wells)
        | map { meta, summary_txt, global_txt, regions_txt, well, plate ->
                [meta.subMap(['donor_id']),
                  meta.id, well, plate, summary_txt, global_txt, regions_txt]
        }
        | groupTuple()

    if (params.run_snps) {
      ch_qc_report_input =
        ch_qc_mosdepth
          | join(concat_snps.out)
          | map { meta, ids, wells, plates, summary_txt, global_txt, regions_txt, set, geno ->
                  [meta, ids, wells, plates, summary_txt, global_txt, regions_txt, true, geno] }
    } else {
      ch_qc_report_input =
        ch_qc_mosdepth
          | map { meta, ids, wells, plates, summary_txt, global_txt, regions_txt ->
                  [meta, ids, wells, plates, summary_txt, global_txt, regions_txt, false, []] }
    }
    knit_qc_report(qc_report_rmd, ch_qc_report_input, params.seq_type)

  }

}