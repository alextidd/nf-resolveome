// bedtools intersect of snps and bait set
process bedtools_intersect_snps {
  tag "${meta.id}"
  label 'normal'

  input:
  tuple val(meta), val(set), path(snps)
  path(bait_set_hyb)

  output:
  tuple val(meta), val(set), path("snps_intersected.tsv")

  script:
  """
  module load bedtools2-2.29.0/python-3.10.10

  awk 'NR==1 {
    # identify columns and save header names
    for (i=1; i<=NF; i++) {
      header[i] = \$i
      if (\$i == "chr") chr_col = i
      else if (\$i == "pos") pos_col = i
      else other_cols[i] = 1
    }

    # print new header: chr, start, pos, then other columns
    printf "chr\\tstart\\tpos"
    for (i=1; i<=NF; i++) {
      if (i != chr_col && i != pos_col && (i in other_cols)) {
        printf "\\t%s", header[i]
      }
    }
    printf "\\n"
    next
  }
  {
    # print new BED-style row
    chr = \$chr_col
    pos = \$pos_col
    printf "%s\\t%d\\t%d", chr, pos-1, pos
    for (i=1; i<=NF; i++) {
      if (i != chr_col && i != pos_col && (i in other_cols)) {
        printf "\\t%s", \$i
      }
    }
  printf "\\n"
  }' ${snps} > snps.bed

  head -1 snps.bed > snps_intersected.tsv.tmp

  bedtools intersect -a <(sed 1d snps.bed) -b ${bait_set_hyb} -wa -wb \\
    >> snps_intersected.tsv.tmp
  
  cut -f1,3- snps_intersected.tsv.tmp > snps_intersected.tsv
  """
}