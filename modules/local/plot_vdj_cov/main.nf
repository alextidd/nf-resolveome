// plot VDJ coverage
process plot_vdj_cov {
  tag "${meta.id}"
  label 'normal10gb'
  publishDir "${params.out_dir}/${meta.donor_id}/cells/${meta.id}/vdj_cov/",
    mode: "copy"
  
  input:
  tuple val(meta), path(regions_bed)
  path(bait_set_vdj)

  output:
  tuple val(meta), path("${meta.id}_*_mean_cov.png")

  script:
  """
  plot_vdj_cov.R \\
    --regions_bed ${regions_bed} \\
    --bait_set ${bait_set_vdj} \\
    --id ${meta.id}
  """
}