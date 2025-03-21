// plot heatmap
process plot_heatmap {
  tag "${meta.donor_id}"
  label 'normal10gb'
  publishDir "${params.out_dir}/${meta.donor_id}/genotyping/${meta.set}/",
    mode: "copy"
  
  input:
  tuple val(meta), path(geno)

  output:
  path("${meta.donor_id}_mut_VAF_heatmap.pdf")

  script:
  """
  plot_heatmap.R \\
    --geno ${geno} \\
    --id ${meta.id}
  """
}