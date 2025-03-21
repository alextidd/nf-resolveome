// report
process report {
  tag "${meta.donor_id}"
  queue "week"
  memory { 200.GB * task.attempt }
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