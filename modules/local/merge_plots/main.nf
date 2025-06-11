// merge plots
process merge_plots {
  tag "${meta.donor_id}"
  label 'normal20gb'
  publishDir "${params.out_dir}/${meta.donor_id}/${publish_subdir}/",
    mode: "copy"
  
  input:
  path rmd
  tuple val(meta), path(plots), path(rmds)
  val title
  val publish_subdir
  
  output:
  path "${meta.donor_id}_plots.html"
  
  script:
  def c_rmds = 'c("' + rmds.join('", "') + '")'
  """
  #!/usr/bin/env Rscript
  
  rmarkdown::render(
    "${rmd}",
    output_dir = "./",
    output_file = "${meta.donor_id}_plots.html",
    params = list(donor_id = "${meta.donor_id}", title = "${title}", rmds = ${c_rmds}))
  """
}