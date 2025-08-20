// plot BAF
process plot_baf {
  tag "${meta.id}"
  label 'normal10gb'
  publishDir "${params.out_dir}/${meta.donor_id}/cells/${meta.id}/genotyping/${set}/",
    mode: "copy"
  
  input:
  tuple val(meta), val(set), path(geno)
  path refcds

  output:
  tuple val(meta), path("${meta.id}_*_plot.png")

  script:
  def arg_baf_chrs = params.baf_chrs ? "--baf_chrs " + params.baf_chrs : ""
  def arg_baf_genes = params.baf_genes ? "--baf_genes " + params.baf_genes : ""
  """
  plot_baf.R \\
    --geno ${geno} \\
    --id ${meta.id} \\
    --refcds ${refcds} \\
    ${arg_baf_chrs} \\
    ${arg_baf_genes}
  """
}