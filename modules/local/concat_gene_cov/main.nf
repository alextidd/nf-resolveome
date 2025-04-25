// concat gene coverage
process concat_gene_cov {
  tag "${meta.id}"
  label 'normal10gb'
  publishDir "${params.out_dir}/${meta.donor_id}/${meta.id}/mosdepth/",
    mode: "copy"
  
  input:
  tuple val(meta), val(gene), path(gene_covs)
  
  output:
  tuple val(meta), path("${meta.id}_gene_cov.tsv")

  script:
  """
  head -1 ${gene_covs[0]} > ${meta.id}_gene_cov.tsv
  for file in ${gene_covs} ; do
    sed 1d \$file >> ${meta.id}_gene_cov.tsv
  done
  """
}