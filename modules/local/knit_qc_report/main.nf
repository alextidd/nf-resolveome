// report
process knit_qc_report {
  tag { "${meta.donor_id}" }
  queue "week"
  memory { 200.GB * task.attempt }
  publishDir { "${params.out_dir}/${meta.donor_id}/qc/" },
    mode: "copy"
  
  input:
  path rmd
  tuple val(meta), val(ids), val(wells), val(plates),
        path(summary_txt), path(global_txt), path(regions_txt),
        val(run_snps), path(geno)
  val(seq_type)

  output:
  path "${meta.donor_id}_qc_report.html"
  path "${meta.donor_id}_metrics_per_cell.tsv"

  script:
  def c_ids = 'c("' + ids.join('", "') + '")'
  def c_wells = 'c("' + wells.join('", "') + '")'
  def c_plates = 'c("' + plates.join('", "') + '")'
  def r_run_snps = run_snps ? 'TRUE' : 'FALSE'
  """
  #!/usr/bin/env Rscript

  # capture params
  my_params <-
    list(
      donor_id = "${meta.donor_id}",
      ids = ${c_ids},
      wells = ${c_wells},
      plates = ${c_plates},
      seq_type = "${seq_type}",
      run_snps = ${r_run_snps})
  saveRDS(my_params, "params.rds")

  # render
  rmarkdown::render(
    "${rmd}",
    output_dir = "./",
    output_file = "${meta.donor_id}_qc_report.html",
    params = my_params)
  """
}