// get gene coordinates
process get_gene_coords {
  tag "${gene}"
  label 'normal'
  
  input:
  val(gene)

  output:
  tuple val(gene), path("${gene}.bed")

  script:
  def chr = params.no_chr ? "" : "chr"
  def host = params.genome_build == "GRCh37" ? "host = 'grch37.ensembl.org'," : ""
  """
  #!/usr/bin/env Rscript

  # libraries
  library(biomaRt)
  
  # use the ensembl database
  ensembl <- useEnsembl(biomart = "genes", $host
                        dataset = "hsapiens_gene_ensembl")

  # get coordinates
  gene_coordinates <- getBM(
    attributes = c("hgnc_symbol", "chromosome_name", "start_position",
                   "end_position", "strand"),
    filters = "hgnc_symbol",
    values = c("${gene}"),
    mart = ensembl) |>
    # filter out non-canonical chromosomes
    dplyr::filter(!grepl("HSCHR|PATCH|CTG", chromosome_name)) |>
    dplyr::transmute(chr = paste0("$chr", chromosome_name),
                     start = start_position, end = end_position,
                     gene = hgnc_symbol)

  # save
  readr::write_tsv(gene_coordinates, "${gene}.bed", col_names = FALSE)
  """
}