// genotype mutations
process genotype_mutations {
  tag "${meta.id}"
  label 'normal10gb'
  
  input:
  tuple val(meta), path(bam), path(bai), val(set), path(mutations)
  
  output:
  tuple val(meta), val(set), path("${meta.id}_genotyped_mutations.tsv")
        
  script:
  """
  genotype_mutations.R \\
    --mutations ${mutations} \\
    --bam ${bam} \\
    --min_bq ${params.min_bq} \\
    --mask ${params.mask} \\
    --min_mq ${params.min_mq} \\
    --id ${meta.id}
  """
}