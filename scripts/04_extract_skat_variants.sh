#!/bin/bash

set -euo pipefail

# ============================================================
# Extract TREM2 LoF and missense variants for SKAT-O
#
# Inputs:
#   1. Final QC PLINK 2 dataset prefix
#   2. LoF variant ID list
#   3. Missense variant ID list
#   4. Output directory
#
# Usage:
# bash scripts/04_extract_skat_variants.sh \
#     results/qc/04_final \
#     results/vep/trem2_lof_variants.txt \
#     results/vep/trem2_missense_variants.txt \
#     results/skat
# ============================================================

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <plink_prefix> <lof_list> <missense_list> <output_dir>"
    exit 1
fi

PLINK="$1"
LOF_LIST="$2"
MISSENSE_LIST="$3"
OUTDIR="$4"


# ------------------------------------------------------------
# Check inputs
# ------------------------------------------------------------

if [ ! -f "${PLINK}.pgen" ] || \
   [ ! -f "${PLINK}.pvar" ] || \
   [ ! -f "${PLINK}.psam" ]; then
    echo "ERROR: PLINK dataset not found: $PLINK"
    exit 1
fi

if [ ! -f "$LOF_LIST" ]; then
    echo "ERROR: LoF variant list not found: $LOF_LIST"
    exit 1
fi

if [ ! -f "$MISSENSE_LIST" ]; then
    echo "ERROR: Missense variant list not found: $MISSENSE_LIST"
    exit 1
fi

if ! command -v plink2 >/dev/null 2>&1; then
    echo "ERROR: plink2 is not installed or not in PATH."
    exit 1
fi

mkdir -p "$OUTDIR"


# ------------------------------------------------------------
# Count expected variants
# ------------------------------------------------------------

LOF_EXPECTED=$(sort -u "$LOF_LIST" | grep -vc '^$' || true)
MISSENSE_EXPECTED=$(sort -u "$MISSENSE_LIST" | grep -vc '^$' || true)

echo "========================================"
echo "TREM2 SKAT-O VARIANT EXTRACTION"
echo "========================================"
echo
echo "Expected LoF variants:      $LOF_EXPECTED"
echo "Expected missense variants: $MISSENSE_EXPECTED"
echo


# ============================================================
# STEP 1: Extract LoF variants
# ============================================================

echo "----------------------------------------"
echo "STEP 1: Extract LoF variants"
echo "----------------------------------------"

plink2 \
    --pfile "$PLINK" \
    --extract "$LOF_LIST" \
    --make-pgen \
    --out "$OUTDIR/trem2_lof"

LOF_EXTRACTED=$(grep -vc '^#' "$OUTDIR/trem2_lof.pvar")

echo "LoF variants extracted: $LOF_EXTRACTED"

if [ "$LOF_EXTRACTED" -ne "$LOF_EXPECTED" ]; then
    echo "ERROR: Expected $LOF_EXPECTED LoF variants but extracted $LOF_EXTRACTED."
    echo "Check whether the variant IDs in the VEP list match the PLINK .pvar IDs."
    exit 1
fi

echo


# ============================================================
# STEP 2: Extract missense variants
# ============================================================

echo "----------------------------------------"
echo "STEP 2: Extract missense variants"
echo "----------------------------------------"

plink2 \
    --pfile "$PLINK" \
    --extract "$MISSENSE_LIST" \
    --make-pgen \
    --out "$OUTDIR/trem2_missense"

MISSENSE_EXTRACTED=$(grep -vc '^#' "$OUTDIR/trem2_missense.pvar")

echo "Missense variants extracted: $MISSENSE_EXTRACTED"

if [ "$MISSENSE_EXTRACTED" -ne "$MISSENSE_EXPECTED" ]; then
    echo "ERROR: Expected $MISSENSE_EXPECTED missense variants but extracted $MISSENSE_EXTRACTED."
    echo "Check whether the variant IDs in the VEP list match the PLINK .pvar IDs."
    exit 1
fi

echo


# ============================================================
# STEP 3: Export additive genotype matrices
#
# PLINK --export A produces:
#   rows    = subjects
#   columns = variants
#   values  = 0, 1, or 2 copies of the counted allele
# ============================================================

echo "----------------------------------------"
echo "STEP 3: Export genotype matrices"
echo "----------------------------------------"

plink2 \
    --pfile "$OUTDIR/trem2_lof" \
    --export A \
    --out "$OUTDIR/trem2_lof_genotypes"

plink2 \
    --pfile "$OUTDIR/trem2_missense" \
    --export A \
    --out "$OUTDIR/trem2_missense_genotypes"

echo


# ============================================================
# Final summary
# ============================================================

echo "========================================"
echo "EXTRACTION COMPLETE"
echo "========================================"
echo
echo "LoF variants:      $LOF_EXTRACTED"
echo "Missense variants: $MISSENSE_EXTRACTED"
echo
echo "PLINK datasets:"
echo "  $OUTDIR/trem2_lof"
echo "  $OUTDIR/trem2_missense"
echo
echo "Genotype matrices:"
echo "  $OUTDIR/trem2_lof_genotypes.raw"
echo "  $OUTDIR/trem2_missense_genotypes.raw"