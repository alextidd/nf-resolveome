// bedtools intersect of a named-region bed and the dnahyb bait set
process bedtools_intersect_vdj_bed {
  label 'normal'

  input:
  path(regions_bed)
  path(dnahyb_bed)
  val(prefix)

  output:
  path("${prefix}_intersected.bed")

  script:
  """
  module load bedtools2-2.29.0/python-3.10.10

  bedtools intersect -a ${regions_bed} -b ${dnahyb_bed} -wa -u > ${prefix}_intersected.bed
  """
}
