// plot VDJ coverage
process plot_vdj_cov {
  tag "${meta.id}"
  label 'normal10gb'
  publishDir "${params.out_dir}/${meta.donor_id}/${meta.id}/vdj_cov/",
    mode: "copy"
  
  input:
  tuple val(meta), path(regions_bed)
  path(bait_set_vdj)

  output:
  path("${meta.id}_*_mean_cov.png"), emit: plot

  script:
  """
  plot_vdj_cov.R \\
    --regions_bed ${regions_bed} \\
    --bait_set ${bait_set_vdj} \\
    --id ${meta.id}
  """
}