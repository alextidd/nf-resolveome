// summarise coverage of the IG regions (CSR, IGKC, IGLC) per cell
process summarise_ig_cov {
  tag { "${meta.donor_id}" }
  label 'normal10gb'
  publishDir { "${params.out_dir}/${meta.donor_id}/vdj/" },
    mode: "copy"

  input:
  tuple val(meta), val(ids), path(regions_beds), path(summary_txts)
  val(seq_type)
  val(csr_max_norm_cov)
  val(light_chain_max_norm_cov)

  output:
  tuple val(meta), path("${meta.donor_id}_ig_cov.tsv")

  script:
  """
  summarise_ig_cov.R \\
    --ids ${ids.join(',')} \\
    --donor_id ${meta.donor_id} \\
    --seq_type ${seq_type} \\
    --csr_max_norm_cov ${csr_max_norm_cov} \\
    --light_chain_max_norm_cov ${light_chain_max_norm_cov}
  """
}
