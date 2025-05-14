// merge BAF plots
process merge_baf_plots {
  tag "${meta.donor_id}"
  label 'normal20gb'
  publishDir "${params.out_dir}/${meta.donor_id}/genotyping/snps/",
    mode: "copy"
  
  input:
  path rmd
  tuple val(meta), path(baf_plots), path(baf_rmds)
  
  output:
  path "${meta.donor_id}_baf_plots.html"
  
  script:
  def c_baf_rmds = 'c("' + baf_rmds.join('", "') + '")'
  """
  #!/usr/bin/env Rscript

  rmarkdown::render(
    "${rmd}",
    output_dir = "./",
    output_file = "${meta.donor_id}_baf_plots.html",
    params = list(baf_rmds = ${c_baf_rmds}))
  """
}