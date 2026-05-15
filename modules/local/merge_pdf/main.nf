// merge plots into a PDF
process merge_pdf {
  tag { "${meta.donor_id}" }
  label 'normal'
  publishDir { "${params.out_dir}/${meta.donor_id}/${publish_subdir}/" },
    mode: "copy"

  input:
  tuple val(meta), path(plots)
  val title
  val publish_subdir

  output:
  path "${meta.donor_id}_${title}.pdf"

  script:
  """
  pdfunite \$(ls *.pdf | sort) ${meta.donor_id}_${title}.pdf
  """
}
