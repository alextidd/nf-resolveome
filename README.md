# nf-resolveome

## Overview

`nf-resolveome` is a Nextflow pipeline designed for performing genotyping, 
VDJ, and LOH analyses for the ResolveDNA/ResolveOME protocol.

1. If `--location irods`, retrieve the BAM from iRODS.
2. Index the BAM.
3. Run MOSDEPTH to calculate genome-wide coverage.
4. Run MOSDEPTH to calculate coverage of the bait set.
5. Run MOSDEPTH to calculate and plot coverage of the the VDJ regions.
6. Optionally, run MOSDEPTH over the IG regions to call class switch recombination and light chain usage.
7. Genotype and optionally annotate mutations.
8. Genotype SNPs and generate BAF plots.
9. Optionally, knit the report.

## Usage

### Samplesheet

First, prepare a comma-delimited samplesheet with your input data. It should 
look like this:

| id                          | donor_id | bam                           | mutations                     | snps                      |
|----------------------------|----------|--------------------------------|-------------------------------|---------------------------|
| plate3_wellA12_dna_run49882 | PD63118  | plate3_wellA12_dna_run49882.bam | PD63118_nanoseq_mutations.tsv | PD63118_caveman_snps.tsv |
| plate3_wellA2_dna_run49882  | PD63118  | plate3_wellA2_dna_run49882.bam  | PD63118_nanoseq_mutations.tsv | PD63118_caveman_snps.tsv |
| plate3_wellA3_dna_run49882  | PD63118  | plate3_wellA3_dna_run49882.bam  | PD63118_nanoseq_mutations.tsv | PD63118_caveman_snps.tsv |
| plate3_wellA4_dna_run49882  | PD63118  | plate3_wellA4_dna_run49882.bam  | PD63118_nanoseq_mutations.tsv | PD63118_caveman_snps.tsv |                                            |

Each row represents a BAM file (`bam`) corresponding to a unique ID (`id`), with 
associated donor ID (`donor_id`), mutations to genotype (`mutations`), 
heterozygous SNPs to plot on the BAF plot (`snps`).

#### Mutations and SNPs

The mutation and SNP files should be tab-delimited and look like this:

```
chr     pos     ref     alt
1       2488104 A       C
1       2488105 T       A
1       2488105 T       C
1       2488106 G       A
1       2488106 G       T
1       2488146 A       AC
```

#### Bait set

The `--dnahyb_bed` should be a bed file containing all targeted regions, like
this:

```
1       1718762 1718882
1       1720479 1720719
1       1721814 1722054
1       1724656 1724776
1       1735818 1736058
1       1737885 1738005
1       1747187 1747307
```

#### VDJ regions

The `--vdj_bed` should be a bed file with named regions whose coverage is
of interest to discern VDJ recombination, like this:

```
2       87565634        87566158        IGKV3OR2-268
2       89156674        89157196        IGKC
2       89160080        89160117        IGKJ5
2       89160398        89160434        IGKJ4
2       89160733        89160770        IGKJ3
2       89161037        89161074        IGKJ2
2       89161398        89161435        IGKJ1
2       89184913        89185669        IGKV4-1
2       89196748        89197300        IGKV5-2
2       89214597        89214894        IGKV7-3
```

#### IG regions

The `--ig_bed` should be a bed file of named regions whose coverage reports
immunoglobulin locus rearrangement, like this:

```
chr14   105836764       105856218       CSR
chr2    88857160        88857683        IGKC
chr22   22895374        22895834        IGLC1
chr22   22900975        22901437        IGLC2
```

The region names are what the pipeline keys on: `CSR` (the IGHD-IGHM span,
which contains Smu, and is deleted when a B cell class switches), `IGKC`
(deleted by the kappa-deleting element in IGL-using cells) and `IGLC1` to
`IGLC7` (the light chain control). Coverage of each is normalised to the
overall mean coverage of the cell, and a region falling below the relevant
threshold is called as deleted. Both calls detect *biallelic* deletion only: a
cell that has rearranged one allele retains about half coverage and is called
`FALSE`. Regions outside the hyb panel are dropped from the bed and reported as
`NA` rather than as deleted.

Note that coverage here excludes duplicates, so the values are not on the same
scale as Picard `CollectHsMetrics` `MEAN_BAIT_COVERAGE`, which counts them.

### Run

Now, you can run the pipeline using, for example:

```bash
nextflow run nf-resolveome \
    --samplesheet samplesheet.csv \
    --dnahyb_bed data/immune_panel.bed \
    --vdj_bed data/ig_tcr_genes.bed \
    --fasta /lustre/scratch124/casm/team78pipelines/canpipe/live/ref/Homo_sapiens/GRCh37d5/genome.fa \
    --location local \
    --baf_chrs 1,9
```

To get a full list of parameters, use the `--help` flag:

```bash
$ nextflow run nf-resolveome --help
```

## Parameters

### Input/output options

- `--location`: Are the BAMs saved locally or on iRODs?  (accepted: irods, local) [default: local]
- `--samplesheet`: Comma-separated file containing the columns 'id', 'donor_id', 'bam', 'mutations', and 'snps'. Also requires 'well' (e.g. 'A1', 'B10') and 'plate' columns if `--knit_qc_report` is set, to plot metrics spatially by well, faceted by plate. If `--knit_genotyping_report` is set, the optional 'well', 'plate' and 'cell_type' columns are used to label the cells in the report. 
- `--min_bq`: Minimum base quality for genotyping. [default: 30] 
- `--min_mq`: Minimum mapping quality for genotyping. [default: 30] 
- `--mask`: Mask for genotyping. [default: 3844] 
- `--out_dir`: Output directory. [default: out/] 
- `--annotate_mutations`: Annotate the genes and impacts of mutations using dndscv? 
- `--dnahyb_bed`: A bed file of the bait set used for hybridisation. 
- `--vdj_bed`: A bed file of the VDJ regions of interest. 
- `--ig_bed`: A bed file of named IG regions (CSR, IGKC, IGLC1-7) whose coverage calls class switch recombination and light chain usage. If unset, the IG coverage step is skipped. 
- `--csr_max_norm_cov`: Coverage of the CSR region below this fraction of the cell's overall mean coverage is called as class switch recombination. [default: 0.02]
- `--light_chain_max_norm_cov`: Coverage of IGKC below this fraction of the cell's overall mean coverage is called as deletion of the IGK locus, i.e. an IGL-using cell. [default: 0.02]
- `--baf_chrs`: Any chromosomes of interest to zoom in on when making the BAF plots, for higher breakpoint resolution, delimited with a comma (e.g. 1,2,3). 
- `--knit_qc_report`: Knit the report? [default: false]
- `--knit_genotyping_report`: Knit the genotyping report? Requires `--run_mutations`. [default: false]
- `--min_alt_vaf`: Minimum alt VAF for a mutation to be reported in the genotyping report. [default: 0.2]
- `--min_alt_depth`: Minimum alt depth for a mutation to be reported in the genotyping report. [default: 3]

### Reference files

- `--genome_build`: Genome build.  (accepted: GRCh37, GRCh38)
- `--no_chr_prefix`: Is there a 'chr' prefix on the chromosome names? [default: false]
- `--refcds`: Path to RefCDS Rda object from the dndscv package.
- `--fasta`: Fasta file for the genome build. 

If you wish to run the pipeline using a different genome build, you will have 
to change these. Alternative RefCDS Rda objects can be downloaded 
[here](https://github.com/im3sanger/dndscv_data/tree/master/data).

**N.B.** Make sure the contig names (e.g. chr1 vs 1) are consistent across the
BAM, the mutations files, the SNPs files, and the reference files.

## Output

An example of an output directory, for one cell of patient `PD64513`:

```
.
└── PD64513
    ├── cells
    │   └── PD64513_plate003_wellA10_CD4.Naive_dnahyb
    │       ├── genotyping
    │       │   └── snps
    │       │       └── ..._snps_baf_plot.pdf
    │       ├── mosdepth
    │       │   ├── ....mosdepth.global.dist.txt
    │       │   ├── ....mosdepth.region.dist.txt
    │       │   ├── ....mosdepth.summary.txt
    │       │   ├── ....regions.bed.gz
    │       │   └── ....regions.bed.gz.csi
    │       └── vdj
    │           ├── ....regions.bed.gz          # VDJ regions
    │           ├── ....ig.regions.bed.gz       # IG regions
    │           └── ..._BCR_chr14_mean_cov.pdf
    ├── genotyping
    │   ├── mutations
    │   │   ├── PD64513_annotated_mutations.tsv
    │   │   ├── PD64513_genotyped_mutations.tsv
    │   │   └── PD64513_genotyping_report.html
    │   └── snps
    │       ├── PD64513_baf_plots.pdf
    │       └── PD64513_genotyped_snps.tsv
    ├── qc
    │   ├── PD64513_metrics_per_cell.tsv
    │   └── PD64513_qc_report.html
    └── vdj
        ├── PD64513_ig_cov.tsv
        └── PD64513_vdj_coverage_plots.pdf
```

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.