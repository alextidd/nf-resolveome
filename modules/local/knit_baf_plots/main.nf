// knit BAF plots
process knit_baf_plots {
  tag "${meta.id}"
  label 'normal'
  
  input:
  tuple val(meta), path(baf_plots)
  
  output:
  tuple val(meta), path(baf_plots), path("${meta.id}_baf_plots.Rmd")
  
  script:
  def c_baf_plots = 'c("' + baf_plots.join('", "') + '")'
  """
  cat <<'EOF' > ${meta.id}_baf_plots.Rmd
  ```{r echo = FALSE, results = 'asis'}
  png_files <- setNames(${c_baf_plots}, tools::file_path_sans_ext(basename(${c_baf_plots})))
  cat("## ${meta.id}\\n\\n")
  for (i in seq_along(png_files)) {
    cat("### ", names(png_files)[i], "\\n\\n")
    cat("![](", png_files[i], ")", "\\n\\n", sep = "")
  }
  ```
  EOF
  """
}