// concat mutations
process concat_mutations {
  tag "${meta.donor_id}"
  label 'normal10gb'
  publishDir "${params.out_dir}/${meta.donor_id}/genotyping/${set}/",
    mode: "copy"
  
  input:
  tuple val(meta), val(set), path(genos)
  
  output:
  tuple val(meta), val(set), path("${meta.donor_id}_genotyped_mutations.tsv")
  
  script:
  """
  head -1 ${genos[0]} > ${meta.donor_id}_genotyped_mutations.tsv
  for file in ${genos} ; do
    sed 1d \$file >> ${meta.donor_id}_genotyped_mutations.tsv
  done
  """
}