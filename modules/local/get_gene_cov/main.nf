// get gene coverage
process get_gene_cov {
  tag "${meta.id}_${gene}"
  label 'normal10gb'
  
  input:
  tuple val(meta), path(per_base_bed), val(gene), path(gene_bed)
  
  output:
  tuple val(meta), val(gene), path("${meta.id}_${gene}_gene_cov.tsv")
  
  script:
  """
  module load bedtools2-2.29.0/python-3.10.10
  echo -e "chr\\tstart\\tend\\tcov\\tgene\\tid" > ${meta.id}_${gene}_gene_cov.tsv
  zcat $per_base_bed |
  bedtools intersect -a stdin -b $gene_bed -wa -wb |
  cut -f1-4,8 |
  awk -v OFS='\\t' '{print \$0, "${meta.id}"}' \\
  >> ${meta.id}_${gene}_gene_cov.tsv
  """
}