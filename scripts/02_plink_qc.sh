#!/bin/bash

set -euo pipefail

# ============================================================
# TREM2 rare-variant quality control
#
# Phenotype coding:
#   DX = 0 -> Control
#   DX = 1 -> Alzheimer's disease case
#
# Usage:
# bash scripts/02_plink_qc.sh \
#     data/input.vcf.gz \
#     data/covariates.txt \
#     results/qc
# ============================================================

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input.vcf.gz> <covariates.txt> <output_dir>"
    exit 1
fi

VCF="$1"
COV="$2"
OUTDIR="$3"

# ------------------------------------------------------------
# Check inputs
# ------------------------------------------------------------

if [ ! -f "$VCF" ]; then
    echo "ERROR: VCF not found: $VCF"
    exit 1
fi

if [ ! -f "$COV" ]; then
    echo "ERROR: Covariate file not found: $COV"
    exit 1
fi

if ! command -v plink2 >/dev/null 2>&1; then
    echo "ERROR: plink2 is not installed or not in PATH."
    exit 1
fi

mkdir -p "$OUTDIR"


# ------------------------------------------------------------
# Function to count variants in a PLINK .pvar file
# ------------------------------------------------------------

count_variants() {
    awk '!/^#/ && NF {n++} END {print n+0}' "$1"
}


echo "========================================"
echo "TREM2 PLINK QC"
echo "========================================"
echo
echo "VCF:        $VCF"
echo "Covariates: $COV"
echo "Output:     $OUTDIR"
echo
echo "Phenotype coding:"
echo "  0 = Control"
echo "  1 = Alzheimer's disease case"
echo


# ============================================================
# STEP 0: Convert VCF to PLINK 2 format
# ============================================================

echo "----------------------------------------"
echo "STEP 0: Convert VCF to PLINK 2 format"
echo "----------------------------------------"

plink2 \
    --vcf "$VCF" \
    --make-pgen \
    --out "$OUTDIR/00_raw"

RAW_COUNT=$(count_variants "$OUTDIR/00_raw.pvar")

echo "Starting variants: $RAW_COUNT"
echo


# ============================================================
# STEP 1: Genotyping rate >= 90%
#
# --geno 0.1 removes variants with >10% missing genotypes.
# ============================================================

echo "----------------------------------------"
echo "STEP 1: Genotyping rate >= 90%"
echo "----------------------------------------"

plink2 \
    --pfile "$OUTDIR/00_raw" \
    --geno 0.1 \
    --make-pgen \
    --out "$OUTDIR/01_geno"

GENO_COUNT=$(count_variants "$OUTDIR/01_geno.pvar")
GENO_REMOVED=$((RAW_COUNT - GENO_COUNT))

echo "Before:  $RAW_COUNT"
echo "After:   $GENO_COUNT"
echo "Removed: $GENO_REMOVED"
echo


# ============================================================
# STEP 2A: Identify controls and cases from covariate file
#
# DX = 0 -> control
# DX = 1 -> case
# ============================================================

echo "----------------------------------------"
echo "STEP 2A: Identify controls and cases"
echo "----------------------------------------"

CONTROL_FILE="$OUTDIR/controls.keep"

awk '
NR == 1 {
    id_col = 0
    dx_col = 0

    for (i = 1; i <= NF; i++) {
        if ($i == "ID") {
            id_col = i
        }

        if ($i == "DX") {
            dx_col = i
        }
    }

    if (id_col == 0 || dx_col == 0) {
        print "ERROR: Could not find ID or DX column." > "/dev/stderr"
        exit 1
    }

    next
}

$dx_col == 0 {
    print $id_col
}
' "$COV" > "$CONTROL_FILE"


CONTROL_COUNT=$(awk '
NR == 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == "DX") {
            dx_col = i
        }
    }
    next
}

$dx_col == 0 {
    count++
}

END {
    print count + 0
}
' "$COV")


CASE_COUNT=$(awk '
NR == 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == "DX") {
            dx_col = i
        }
    }
    next
}

$dx_col == 1 {
    count++
}

END {
    print count + 0
}
' "$COV")


echo "Controls (DX=0): $CONTROL_COUNT"
echo "Cases    (DX=1): $CASE_COUNT"

if [ "$CONTROL_COUNT" -eq 0 ]; then
    echo "ERROR: No controls with DX=0 were found."
    exit 1
fi

if [ "$CASE_COUNT" -eq 0 ]; then
    echo "ERROR: No cases with DX=1 were found."
    exit 1
fi

echo


# ============================================================
# STEP 2B: HWE P > 1e-5 in controls only
#
# Controls are temporarily selected to determine which variants
# pass the HWE threshold.
# ============================================================

echo "----------------------------------------"
echo "STEP 2B: HWE P > 1e-5 in controls"
echo "----------------------------------------"

plink2 \
    --pfile "$OUTDIR/01_geno" \
    --keep "$CONTROL_FILE" \
    --hwe 1e-5 \
    --make-pgen \
    --out "$OUTDIR/02_hwe_controls"


# Get list of variants that passed HWE in controls

plink2 \
    --pfile "$OUTDIR/02_hwe_controls" \
    --write-snplist \
    --out "$OUTDIR/02_hwe_pass"


# Apply passing variants back to ALL subjects

plink2 \
    --pfile "$OUTDIR/01_geno" \
    --extract "$OUTDIR/02_hwe_pass.snplist" \
    --make-pgen \
    --out "$OUTDIR/02_hwe"

HWE_COUNT=$(count_variants "$OUTDIR/02_hwe.pvar")
HWE_REMOVED=$((GENO_COUNT - HWE_COUNT))

echo "Before:  $GENO_COUNT"
echo "After:   $HWE_COUNT"
echo "Removed: $HWE_REMOVED"
echo


# ============================================================
# STEP 3: Minor allele count > 1
#
# MAC > 1 means retain variants with MAC >= 2.
# ============================================================

echo "----------------------------------------"
echo "STEP 3: Minor allele count > 1"
echo "----------------------------------------"

plink2 \
    --pfile "$OUTDIR/02_hwe" \
    --mac 2 \
    --make-pgen \
    --out "$OUTDIR/03_mac"

MAC_COUNT=$(count_variants "$OUTDIR/03_mac.pvar")
MAC_REMOVED=$((HWE_COUNT - MAC_COUNT))

echo "Before:  $HWE_COUNT"
echo "After:   $MAC_COUNT"
echo "Removed: $MAC_REMOVED"
echo


# ============================================================
# STEP 4: Minor allele frequency < 1%
# ============================================================

echo "----------------------------------------"
echo "STEP 4: Minor allele frequency < 1%"
echo "----------------------------------------"

plink2 \
    --pfile "$OUTDIR/03_mac" \
    --max-maf 0.01 \
    --make-pgen \
    --out "$OUTDIR/04_final"

FINAL_COUNT=$(count_variants "$OUTDIR/04_final.pvar")
MAF_REMOVED=$((MAC_COUNT - FINAL_COUNT))

echo "Before:  $MAC_COUNT"
echo "After:   $FINAL_COUNT"
echo "Removed: $MAF_REMOVED"
echo


# ============================================================
# QC summary
# ============================================================

TOTAL_REMOVED=$((RAW_COUNT - FINAL_COUNT))

SUMMARY="$OUTDIR/qc_summary.tsv"

{
    echo -e "Step\tBefore\tAfter\tRemoved"
    echo -e "Genotyping_rate\t$RAW_COUNT\t$GENO_COUNT\t$GENO_REMOVED"
    echo -e "HWE_controls\t$GENO_COUNT\t$HWE_COUNT\t$HWE_REMOVED"
    echo -e "MAC_gt_1\t$HWE_COUNT\t$MAC_COUNT\t$MAC_REMOVED"
    echo -e "MAF_lt_0.01\t$MAC_COUNT\t$FINAL_COUNT\t$MAF_REMOVED"
} > "$SUMMARY"


echo "========================================"
echo "QC COMPLETE"
echo "========================================"
echo
echo "Starting variants: $RAW_COUNT"
echo "Final variants:    $FINAL_COUNT"
echo "Total removed:     $TOTAL_REMOVED"
echo
echo "QC summary:"
echo "$SUMMARY"
echo
echo "Final PLINK dataset:"
echo "$OUTDIR/04_final"