// plot BAF
process plot_baf {
  tag "${meta.id}"
  label 'normal10gb'
  publishDir "${params.out_dir}/${meta.donor_id}/${meta.id}/genotyping/${set}/",
    mode: "copy"
  
  input:
  tuple val(meta), val(set), path(geno)

  output:
  tuple val(meta), path(geno), emit: out
  path("${meta.id}_*_baf_plot.png"), emit: plot

  script:
  """
  plot_baf.R \\
    --geno ${geno} \\
    --id ${meta.id}
  """
}