// report
process knit_genotyping_report {
  tag { "${meta.donor_id}" }
  queue "week"
  memory { 200.GB * task.attempt }
  publishDir { "${params.out_dir}/${meta.donor_id}/genotyping/mutations/" },
    mode: "copy"

  input:
  path rmd
  tuple val(meta), val(ids), val(wells), val(plates), val(cell_types), path(geno)
  val(min_alt_vaf)
  val(min_alt_depth)

  output:
  path "${meta.donor_id}_genotyping_report.html"

  script:
  def c_ids = 'c("' + ids.join('", "') + '")'
  def c_wells = 'c("' + wells.join('", "') + '")'
  def c_plates = 'c("' + plates.join('", "') + '")'
  def c_cell_types = 'c("' + cell_types.join('", "') + '")'
  """
  #!/usr/bin/env Rscript

  # capture params
  my_params <-
    list(
      donor_id = "${meta.donor_id}",
      ids = ${c_ids},
      wells = ${c_wells},
      plates = ${c_plates},
      cell_types = ${c_cell_types},
      min_alt_vaf = ${min_alt_vaf},
      min_alt_depth = ${min_alt_depth})
  saveRDS(my_params, "params.rds")

  # render
  rmarkdown::render(
    "${rmd}",
    output_dir = "./",
    output_file = "${meta.donor_id}_genotyping_report.html",
    params = my_params)
  """
}
