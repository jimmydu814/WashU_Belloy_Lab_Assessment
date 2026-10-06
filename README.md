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

scirpts: 01_inspect_vcf.sh and 01b_inspecct_covarates.sh

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

Commend run:

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

Script: 01_inspect_vcf.sh

Commend run:
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

Scirpt: 02_plink_qc.sh

Commend ran:
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
