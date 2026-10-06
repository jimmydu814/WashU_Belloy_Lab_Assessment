#!/bin/bash

set -euo pipefail

# Usage:
# bash scripts/01_inspect_vcf.sh data/input.vcf.gz

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input.vcf.gz>"
    exit 1
fi

VCF="$1"

# Check that the input file exists
if [ ! -f "$VCF" ]; then
    echo "ERROR: File not found: $VCF"
    exit 1
fi

# Check that bcftools is installed
if ! command -v bcftools >/dev/null 2>&1; then
    echo "ERROR: bcftools is not installed or not in PATH."
    exit 1
fi

echo "========================================"
echo "VCF INPUT INSPECTION"
echo "========================================"
echo

echo "Input file:"
echo "$VCF"
echo


echo "----------------------------------------"
echo "1. Genome build / reference information"
echo "----------------------------------------"

BUILD_INFO=$(bcftools view -h "$VCF" | grep -Ei 'reference=|assembly=|GRCh37|GRCh38|hg19|hg38' || true)

if [ -z "$BUILD_INFO" ]; then
    echo "No obvious genome build information found in VCF header."
else
    echo "$BUILD_INFO"
fi

echo


echo "----------------------------------------"
echo "2. Number of samples"
echo "----------------------------------------"

SAMPLE_COUNT=$(bcftools query -l "$VCF" | wc -l)

echo "Samples: $SAMPLE_COUNT"
echo


echo "----------------------------------------"
echo "3. First 10 sample IDs"
echo "----------------------------------------"

bcftools query -l "$VCF" | sed -n '1,10p'

echo


echo "----------------------------------------"
echo "4. Number of variant records"
echo "----------------------------------------"

VARIANT_COUNT=$(bcftools view -H "$VCF" | wc -l)

echo "Variants: $VARIANT_COUNT"
echo


echo "----------------------------------------"
echo "5. Chromosomes present"
echo "----------------------------------------"

bcftools query -f '%CHROM\n' "$VCF" | sort -u

echo


echo "----------------------------------------"
echo "6. Genotype (GT) field"
echo "----------------------------------------"

GT_INFO=$(bcftools view -h "$VCF" | grep '^##FORMAT=<ID=GT' || true)

if [ -z "$GT_INFO" ]; then
    echo "WARNING: GT field was not found in the VCF header."
else
    echo "$GT_INFO"
fi

echo


echo "----------------------------------------"
echo "7. Number of multiallelic variants"
echo "----------------------------------------"

MULTI_COUNT=$(bcftools view -H -m3 "$VCF" | wc -l)

echo "Multiallelic variants: $MULTI_COUNT"
echo


echo "----------------------------------------"
echo "8. FILTER values"
echo "----------------------------------------"

bcftools query -f '%FILTER\n' "$VCF" | sort | uniq -c

echo


echo "----------------------------------------"
echo "9. First 10 variants"
echo "----------------------------------------"

printf "CHROM\tPOS\tID\tREF\tALT\tFILTER\n"

bcftools query \
    -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER\n' \
    "$VCF" | head -10

echo

echo "========================================"
echo "Inspection complete"
echo "========================================"