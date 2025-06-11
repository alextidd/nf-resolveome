// knit plots
process knit_plots {
  tag "${meta.id}"
  label 'normal'
  
  input:
  tuple val(meta), path(plots)
  
  output:
  tuple val(meta), path(plots), path("${meta.id}_plots.Rmd")
  
  script:
  def c_plots = 'c("' + plots.join('", "') + '")'
  """
  cat <<'EOF' > ${meta.id}_plots.Rmd
  ```{r echo = FALSE, results = 'asis'}
  png_files <- setNames(${c_plots}, tools::file_path_sans_ext(basename(${c_plots})))
  cat("## ${meta.id}\\n\\n")
  for (i in seq_along(png_files)) {
    cat("### ", names(png_files)[i], "\\n\\n")
    cat("![](", png_files[i], ")", "\\n\\n", sep = "")
  }
  ```
  EOF
  """
}