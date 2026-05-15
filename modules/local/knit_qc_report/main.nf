// report
process knit_qc_report {
  tag { "${meta.donor_id}" }
  queue "week"
  memory { 200.GB * task.attempt }
  publishDir { "${params.out_dir}/${meta.donor_id}/" },
    mode: "copy"
  
  input:
  path rmd
  tuple val(meta), val(ids),
        path(summary_txt), path(global_txt), path(regions_txt),
        val(set), path(geno)
  val(seq_type)

  output:
  path "${meta.donor_id}_qc_report.html"
  
  script:
  def c_ids = 'c("' + ids.join('", "') + '")'
  """
  #!/usr/bin/env Rscript

  # capture params
  my_params <-
    list(
      donor_id = "${meta.donor_id}",
      ids = ${c_ids},
      seq_type = "${seq_type}")
  saveRDS(my_params, "params.rds")

  # render
  rmarkdown::render(
    "${rmd}",
    output_dir = "./",
    output_file = "${meta.donor_id}_qc_report.html",
    params = my_params)
  """
}