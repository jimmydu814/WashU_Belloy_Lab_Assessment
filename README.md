# TREM2 Rare Variant Association Analysis Using SKAT-O

## Project Overview

This project evaluates whether rare functional variants in the **TREM2**
gene are associated with Alzheimer's disease using the optimized
Sequence Kernel Association Test (SKAT-O).

The analysis workflow consists of:

1. Input data inspection
2. Variant quality control using PLINK 2
3. Functional annotation using Ensembl VEP
4. Classification of variants into:
   - Loss-of-function variants
   - Missense variants
5. SKAT-O association testing
6. Interpretation and reporting

The supplied genotype data include variants within TREM2 and 5 kb
beyond the gene boundaries.

## Phenotype Coding

The covariate file uses:

- `DX = 0`: Control
- `DX = 1`: Alzheimer's disease case

The `ID` column contains subject identifiers.

All supplied covariates will be included in the final SKAT-O model.

## Software

The analysis is being performed in WSL/Linux using:

- PLINK 2
- bcftools
- Bash
- Ensembl VEP
- R
- SKAT R package

## Input Data Inspection

Before quality control, both the VCF and covariate files were inspected.

### VCF checks

The following were checked:

- Number of samples
- Number of variants
- Sample identifiers
- Chromosomes represented
- Presence of genotype (`GT`) information
- Multiallelic variants
- Existing VCF filter values
- Variant coordinates and alleles
- Genome build information

The VCF header did not explicitly specify the reference genome build.

Script: 01_inspect_vcf.sh

Commend used:

```bash
bash scripts/01_inspect_vcf.sh data/WES_TREM2_5k_boundaries_raw.vcf.gz | tee results/vcf_inspection.txt
```

### Covariate checks

The covariate file was inspected for:

- Column names
- Subject count
- Phenotype values
- Duplicate subject IDs
- Missing values
- Concordance between covariate IDs and VCF sample IDs

VCF and covariate sample identifiers were found to be consistent.

Script: 01b_inspect_covariates.sh

Commend used:

```bash
bash scripts/01b_inspect_covaraites.sh data/covariates.txt | tee covariates_inspection.txt
```
## Variant Quality Control

Variant quality control was performed using PLINK 2.

The required filtering criteria were:

- Genotyping rate >= 90%
- Hardy-Weinberg equilibrium in controls: P > 1e-5
- Minor allele count (MAC) > 1
- Minor allele frequency (MAF) < 1%

Filters were applied sequentially so that the number of variants removed
at each stage could be documented.

Script: 02_plink_qc.sh

Commend used:
```bash
bash scripts/02_plink_qc.sh data/WES_TREM2_5k_boundaries_raw.vcf.gz data/covariates.txt results/qc
```

### QC Results

| QC Step | Before | After | Removed |
|---|---:|---:|---:|
| Genotyping rate >= 90% | 279 | 163 | 116 |
| HWE in controls, P > 1e-5 | 163 | 163 | 0 |
| MAC > 1 | 163 | 34 | 129 |
| MAF < 1% | 34 | 33 | 1 |

A total of **246 variants were removed**, leaving **33 variants** after
quality control.

### Export Final QC Variants to VCF

```bash
plink2 \
    --pfile results/qc/04_final \
    --export vcf bgz \
    --out results/qc/04_final
```

## Variant Annotation with Ensembl VEP

After PLINK QC, 33 variants remained in the final dataset. The filtered PLINK dataset was exported to VCF format and annotated using Ensembl Variant Effect Predictor (VEP) release 116 with the GRCh38 cache.
VEP commend:
```bash
~/ensembl-vep/vep \
    --input_file results/qc/04_final.vcf.gz \
    --output_file results/vep/vep_annotations.tsv \
    --cache \
    --offline \
    --assembly GRCh38 \
    --dir_cache /mnt/d/vep_cache \
    --symbol \
    --canonical \
    --mane \
    --biotype \
    --flag_pick_allele_gene \
    --tab \
    --fields "Uploaded_variation,Location,Allele,Gene,SYMBOL,Feature,Feature_type,BIOTYPE,Consequence,MANE_SELECT,CANONICAL,PICK,Protein_position,Amino_acids" \
    --force_overwrite
```

### Filtering Variant for SKAT-O
VEP annotations were filtered to keep only TREM2 transcript annotations. Variants were classified as LoF if any TREM2 transcript had a specified LoF consequence, or as missense if no LoF consequence was present but a missense annotation was found. Each genomic variant was counted once, with LoF given priority to avoid double counting.

Script: 03_filter_vep.sh

Commend used: 
```bash
bash scripts/03_filter_vep.sh \
    results/vep/vep_annotations.tsv \
    results/vep
``` 
### Extracting Variant Sets for SKAT-O
The LoF and missense variant lists were extracted from the final QC PLINK dataset and converted to additive genotype matrices for use in R.

Script: 04_extract_skat_variants.sh

Command used:

```bash
bash scripts/04_extract_skat_variants.sh \
    results/qc/04_final \
    results/vep/trem2_lof_variants.txt \
    results/vep/trem2_missense_variants.txt \
    results/skat
```

### SKAT-O Analysis

SKAT-O was run in R using the `SKAT` package for the LoF and missense variant sets. Alzheimer’s disease status (`DX`) was modeled as a binary outcome and adjusted for sex, APOE ε4 dosage, APOE ε2 dosage, and PC1–PC5.

Model used:

```text
DX ~ sex + final_APOE4d + final_APOE2d + PC1 + PC2 + PC3 + PC4 + PC5
```

Results(`results/skat/skato_results.tsv`):
| Variant category | Number of variants included | P-value |
|---|---:|---:|
| LoF | 3 | 1.987 × 10^-14 |
| Missense | 16 | 1.594 × 10^-12 |

Scripts: 05_run_skato.R

Commend used:

```bash
Rscript scripts/05_run_skato.R \
    data/covariates.txt \
    results/skat/trem2_lof_genotypes.raw \
    results/skat/trem2_missense_genotypes.raw \
    results/skat \
    2>&1 | tee results/skat/05_run_skato.log
```

### Scaling to Many Genes
To scale the analysis to many genes, I would package the workflow and software dependencies into a container using Docker. The current workflow starts with a VCF containing only TREM2 and its surrounding 5 kb region, so I would add an initial step to extract the region for each gene of interest from a larger VCF. I would also modify the downstream scripts so that the gene name is provided as an input rather than being hard-coded as TREM2.

For parallel processing, each gene could be analyzed as an independent job on a high-performance computing cluster using a job scheduler such as SLURM. A SLURM job array could run the same QC, VEP annotation, functional classification, and SKAT-O workflow for many genes in parallel, with each job producing its own results and log files.

## Code and Data Documentation

This GitHub repository reflects how I would document and share the analysis. The workflow is organized into separate scripts for data inspection, PLINK QC, VEP annotation/filtering, genotype extraction, and SKAT-O analysis, with the corresponding summary results and log files retained for reproducibility.

The raw genotype and covariate data are not included in the repository because they were provided specifically for this assessment and may contain sensitive or restricted information. Instead, the repository contains the code, documentation, compact result summaries, variant lists, software/version information, and a README describing how to reproduce the analysis when the input data are available. Intermediate result files are also excluded in the repository. See the README file in the results/ folder for details on which outputs are retained and which are excluded.