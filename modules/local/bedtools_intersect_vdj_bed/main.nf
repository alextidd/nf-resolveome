// bedtools intersect of vdj bed and dnahyb bait set
process bedtools_intersect_vdj_bed {
  label 'normal'

  input:
  path(vdj_bed)
  path(dnahyb_bed)

  output:
  path("vdj_intersected.bed")

  script:
  """
  module load bedtools2-2.29.0/python-3.10.10

  bedtools intersect -a ${vdj_bed} -b ${dnahyb_bed} -wa -u > vdj_intersected.bed
  """
}
