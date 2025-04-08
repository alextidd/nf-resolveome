// concat snps per cell
process concat_snps_per_cell {
  tag "${meta.id}"
  label 'normal10gb'
  publishDir "${params.out_dir}/${meta.donor_id}/${meta.id}/genotyping/${set}/",
    mode: "copy"
  
  input:
  tuple val(meta), val(set), path(genos, stageAs: "?/*")
  
  output:
  tuple val(meta), val(set), path("${meta.id}_genotyped_snps.tsv")

  script:
  """
  head -1 ${genos[0]} > ${meta.id}_genotyped_snps.tsv
  for file in ${genos} ; do
    sed 1d \$file >> ${meta.id}_genotyped_snps.tsv
  done
  """
}