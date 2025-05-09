# nf-resolveome

## Overview

`nf-resolveome` is a Nextflow pipeline designed for performing genotyping, 
VDJ, and LOH analyses for the ResolveDNA/ResolveOME protocol.

1. If `--location irods`, retrieve the BAM from iRODS.
2. Index the BAM.
3. Run MOSDEPTH to calculate genome-wide coverage.
4. Run MOSDEPTH to calculate coverage of the bait set.
5. Run MOSDEPTH to calculate and plot coverage of the the VDJ regions.
6. Genotype and optionally annotate mutations.
7. Genotype SNPs and generate BAF plots.
8. Optionally, knit the report.

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

The `--bait_set_hyb` should be a bed file containing all targeted regions, like
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

The `--bait_set_vdj` should be a bed file with named regions whose coverage is
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

### Run

Now, you can run the pipeline using, for example:

```bash
nextflow run nf-resolveome \
    --samplesheet samplesheet.csv \
    --bait_set_hyb data/immune_panel.bed \
    --bait_set_vdj data/ig_tcr_genes.bed \
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
- `--samplesheet`: Comma-separated file containing the columns 'id', 'donor_id', 'bam', 'mutations', and 'snps'. 
- `--min_bq`: Minimum base quality for genotyping. [default: 30] 
- `--min_mq`: Minimum mapping quality for genotyping. [default: 30] 
- `--mask`: Mask for genotyping. [default: 3844] 
- `--out_dir`: Output directory. [default: out/] 
- `--annotate_mutations`: Annotate the genes and impacts of mutations using dndscv? 
- `--bait_set_hyb`: A bed file of the bait set used for hybridisation. 
- `--bait_set_vdj`: A bed file of the VDJ regions of interest. 
- `--baf_chrs`: Any chromosomes of interest to zoom in on when making the BAF plots, for higher breakpoint resolution, delimited with a comma (e.g. 1,2,3). 
- `--knit_report`: Knit the report? [default: false]

### Reference files

- `--genome_build`: Genome build.  (accepted: GRCh37, GRCh38) [default: GRCh37] 
- `--no_chr`: Is there a 'chr' prefix on the chromosome names? [default: false]
- `--refcds`: Path to RefCDS Rda object from the dndscv package.
- `--fasta`: Fasta file for the genome build. 

If you wish to run the pipeline using a different genome build, you will have 
to change these. Alternative RefCDS Rda objects can be downloaded 
[here](https://github.com/im3sanger/dndscv_data/tree/master/data).

**N.B.** Make sure the contig names (e.g. chr1 vs 1) are consistent across the
BAM, the mutations files, the SNPs files, and the reference files.

## Output

An example of an output directory, for the sample `plate3_wellA12_dna_run49882`
in patient `PD63118`:

```
.
└── PD63118
    ├── genotyping
    │   ├── mutations
    │   │   ├── PD63118_annotated_mutations.tsv
    │   │   └── PD63118_genotyped_mutations.tsv
    │   └── snps
    │       └── PD63118_genotyped_snps.tsv
    └── plate3_wellA12_dna_run49882
        ├── genotyping
        │   └── snps
        │       ├── plate3_wellA12_dna_run49882_caveman_snps_baf_chr1_plot.png
        │       ├── plate3_wellA12_dna_run49882_caveman_snps_baf_chr9_plot.png
        │       ├── plate3_wellA12_dna_run49882_caveman_snps_baf_plot.png
        │       └── plate3_wellA12_dna_run49882_genotyped_snps.tsv
        ├── mosdepth
        │   ├── plate3_wellA12_dna_run49882_gene_cov.tsv
        │   ├── plate3_wellA12_dna_run49882.mosdepth.global.dist.txt
        │   ├── plate3_wellA12_dna_run49882.mosdepth.region.dist.txt
        │   ├── plate3_wellA12_dna_run49882.mosdepth.summary.txt
        │   ├── plate3_wellA12_dna_run49882.per-base.bed.gz
        │   ├── plate3_wellA12_dna_run49882.per-base.bed.gz.csi
        │   ├── plate3_wellA12_dna_run49882.regions.bed.gz
        │   ├── plate3_wellA12_dna_run49882.regions.bed.gz.csi
        │   └── versions.yml
        └── vdj_cov
            ├── plate3_wellA12_dna_run49882_chr14_BCR_mean_cov.png
            ├── plate3_wellA12_dna_run49882_chr14_TCR_mean_cov.png
            ├── plate3_wellA12_dna_run49882_chr22_BCR_mean_cov.png
            ├── plate3_wellA12_dna_run49882_chr2_BCR_mean_cov.png
            └── plate3_wellA12_dna_run49882_chr7_TCR_mean_cov.png
```

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.