#!/bin/bash

set -euo pipefail

# ============================================================
# Classify TREM2 variants from VEP annotations
#
# Classification rule:
#   LoF if ANY TREM2 transcript contains:
#     splice_acceptor_variant
#     splice_donor_variant
#     stop_gained
#     frameshift_variant
#
#   Otherwise Missense if ANY TREM2 transcript contains:
#     missense_variant
#
# Each genomic variant is classified only once.
#
# Usage:
# bash scripts/03_filter_vep.sh \
#     results/vep/vep_annotations.tsv \
#     results/vep
# ============================================================

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <vep_annotations.tsv> <output_dir>"
    exit 1
fi

VEP="$1"
OUTDIR="$2"

if [ ! -f "$VEP" ]; then
    echo "ERROR: VEP file not found: $VEP"
    exit 1
fi

mkdir -p "$OUTDIR"

ALL_TREM2="$OUTDIR/trem2_all_transcripts.tsv"
CLASS_SUMMARY="$OUTDIR/trem2_variant_classes.tsv"
LOF="$OUTDIR/trem2_lof_variants.txt"
MISSENSE="$OUTDIR/trem2_missense_variants.txt"


# ============================================================
# STEP 1: Keep all TREM2 transcript annotations
# ============================================================

echo "STEP 1: Selecting all TREM2 transcript annotations"

awk -F'\t' '
BEGIN {
    OFS="\t"
}

$0 ~ /^#Uploaded_variation/ {

    for (i=1; i<=NF; i++) {
        name=$i
        sub(/^#/, "", name)
        col[name]=i
    }

    print
    next
}

!/^#/ &&
$(col["SYMBOL"]) == "TREM2" &&
$(col["Feature_type"]) == "Transcript" {
    print
}
' "$VEP" > "$ALL_TREM2"


# ============================================================
# STEP 2: Classify each UNIQUE variant
# ============================================================

echo "STEP 2: Classifying unique TREM2 variants"

awk -F'\t' \
    -v lof_file="$LOF" \
    -v missense_file="$MISSENSE" \
    -v summary_file="$CLASS_SUMMARY" '

BEGIN {
    OFS="\t"
}

$0 ~ /^#Uploaded_variation/ {

    for (i=1; i<=NF; i++) {
        name=$i
        sub(/^#/, "", name)
        col[name]=i
    }

    next
}

!/^#/ {

    variant=$(col["Uploaded_variation"])
    consequence=$(col["Consequence"])

    variants[variant]=1

    if (consequence ~ /(^|[&,])(splice_acceptor_variant|splice_donor_variant|stop_gained|frameshift_variant)([&,]|$)/) {
        has_lof[variant]=1
    }

    if (consequence ~ /(^|[&,])missense_variant([&,]|$)/) {
        has_missense[variant]=1
    }
}

END {

    print "Variant", "Has_LoF", "Has_Missense", "Final_Class" > summary_file

    for (variant in variants) {

        lof = (variant in has_lof) ? 1 : 0
        miss = (variant in has_missense) ? 1 : 0

        if (lof) {
            class="LoF"
            print variant > lof_file
        }
        else if (miss) {
            class="Missense"
            print variant > missense_file
        }
        else {
            class="Neither"
        }

        print variant, lof, miss, class > summary_file
    }
}
' "$ALL_TREM2"


# ============================================================
# Count results
# ============================================================

TOTAL_COUNT=$(awk 'NR>1 {count++} END {print count+0}' "$CLASS_SUMMARY")
LOF_COUNT=$(wc -l < "$LOF" 2>/dev/null || echo 0)
MISSENSE_COUNT=$(wc -l < "$MISSENSE" 2>/dev/null || echo 0)

echo
echo "========================================"
echo "CLASSIFICATION COMPLETE"
echo "========================================"
echo
echo "Unique TREM2 variants: $TOTAL_COUNT"
echo "LoF variants:          $LOF_COUNT"
echo "Missense variants:     $MISSENSE_COUNT"
echo
echo "Outputs:"
echo "  $ALL_TREM2"
echo "  $CLASS_SUMMARY"
echo "  $LOF"
echo "  $MISSENSE"