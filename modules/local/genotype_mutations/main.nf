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
  export R_LIBS_USER="/software/team268/at31/R/R-tmp-4.4/"
  export R_HOME=/software/isg/languages/R/4.4.0/exec/lib/R
  export LD_LIBRARY_PATH="/software/isg/languages/R/4.4.0/exec/lib/R/lib:/software/isg/languages/R/4.4.0/exec/lib:$LD_LIBRARY_PATH"

  genotype_mutations.R \\
    --mutations ${mutations} \\
    --bam ${bam} \\
    --min_bq ${params.min_bq} \\
    --mask ${params.mask} \\
    --min_mq ${params.min_mq} \\
    --id ${meta.id}
  """
}