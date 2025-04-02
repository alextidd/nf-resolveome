// annotate mutations
process annotate_mutations {
  tag "${meta.donor_id}"
  label 'normal10gb'
  publishDir "${params.out_dir}/${meta.donor_id}/genotyping/${set}/",
    mode: "copy"
    
  input:
  tuple val(meta), val(set), path(mutations)
  path(refcds)

  output:
  tuple val(meta), val(set), path("${meta.donor_id}_annotated_mutations.tsv")
  
  script:
  """
  annotate_mutations.R \\
    --mutations ${mutations} \\
    --donor_id ${meta.donor_id} \\
    --refcds ${refcds}
  """
}