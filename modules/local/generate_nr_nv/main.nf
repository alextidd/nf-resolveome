// generate NR and NV matrices
process generate_nr_nv {
  tag "${meta.donor_id}"
  label 'normal10gb'
  publishDir "${params.out_dir}/${meta.donor_id}/genotyping/${set}/",
    mode: "copy"
  
  input:
  tuple val(meta), val(set), val(ids), path(genos)
  
  output:
  path(geno), emit: geno
  tuple val(meta),
        path("${meta.donor_id}_NV_genotyped_mutations.tsv"),
        path("${meta.donor_id}_NR_genotyped_mutations.tsv"), emit: out
  
  script:
  """
  # get mut ids as row names
  awk -F"\\t" 'NR == 1 {print ""} ; NR > 1 && \$21 == "nanoseq_mutations" {print \$1"-"\$2"-"\$3"-"\$4}' ${genos[0]} \
  > tmp_1

  # for all files get the NR and NV values
  genos=(${genos.join(' ')})
  ids=(${ids.join(' ')})
  for i in "\${!genos[@]}" ; do
    echo \${ids[i]}
    awk -F"\\t" -v id=\${ids[i]} 'NR == 1 {print id} ; NR > 1 {print \$5}' \${genos[i]} \
    > tmp_\${ids[i]}_nv
    awk -F"\\t" -v id=\${ids[i]} 'NR == 1 {print id} ; NR > 1 {print \$7}' \${genos[i]} \
    > tmp_\${ids[i]}_nr
  done

  # paste together
  paste tmp_1 tmp_*_nv > ${meta.donor_id}_NV_genotyped_mutations.tsv
  paste tmp_1 tmp_*_nr > ${meta.donor_id}_NR_genotyped_mutations.tsv
  """
}