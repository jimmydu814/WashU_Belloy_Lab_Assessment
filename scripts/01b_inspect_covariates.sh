#!/bin/bash

set -euo pipefail

# Usage:
# bash scripts/01b_inspect_covariates.sh data/covariates.txt
#
# Optional VCF comparison:
# bash scripts/01b_inspect_covariates.sh data/covariates.txt data/input.vcf.gz

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
    echo "Usage: $0 <covariates.txt> [input.vcf.gz]"
    exit 1
fi

COV="$1"
VCF="${2:-}"

# Check input file
if [ ! -f "$COV" ]; then
    echo "ERROR: Covariate file not found: $COV"
    exit 1
fi

echo "========================================"
echo "COVARIATE FILE INSPECTION"
echo "========================================"
echo

echo "Input file:"
echo "$COV"
echo


echo "----------------------------------------"
echo "1. Header / column names"
echo "----------------------------------------"

sed -n '1p' "$COV"

echo


echo "----------------------------------------"
echo "2. First 5 data rows"
echo "----------------------------------------"

sed -n '2,6p' "$COV"

echo


echo "----------------------------------------"
echo "3. Number of subjects"
echo "----------------------------------------"

# Subtract one for the header
TOTAL_LINES=$(wc -l < "$COV")
SUBJECT_COUNT=$((TOTAL_LINES - 1))

echo "Subjects: $SUBJECT_COUNT"

echo


echo "----------------------------------------"
echo "4. Number of columns"
echo "----------------------------------------"

COLUMN_COUNT=$(awk 'NR==1 {print NF}' "$COV")

echo "Columns: $COLUMN_COUNT"

echo


echo "----------------------------------------"
echo "5. Identify ID and DX columns"
echo "----------------------------------------"

awk '
NR==1 {
    id_col = 0
    dx_col = 0

    for (i=1; i<=NF; i++) {
        if ($i == "ID") {
            id_col = i
        }

        if ($i == "DX") {
            dx_col = i
        }
    }

    if (id_col == 0) {
        print "WARNING: ID column not found"
    } else {
        print "ID column:", id_col
    }

    if (dx_col == 0) {
        print "WARNING: DX column not found"
    } else {
        print "DX column:", dx_col
    }
}
' "$COV"

echo


echo "----------------------------------------"
echo "6. DX phenotype values"
echo "----------------------------------------"

awk '
NR==1 {
    for (i=1; i<=NF; i++) {
        if ($i == "DX") {
            dx_col = i
        }
    }
    next
}

dx_col {
    counts[$dx_col]++
}

END {
    for (value in counts) {
        print value, counts[value]
    }
}
' "$COV" | sort

echo


echo "----------------------------------------"
echo "7. First 10 subject IDs"
echo "----------------------------------------"

awk '
NR==1 {
    for (i=1; i<=NF; i++) {
        if ($i == "ID") {
            id_col = i
        }
    }
    next
}

id_col && NR <= 11 {
    print $id_col
}
' "$COV"

echo


echo "----------------------------------------"
echo "8. Duplicate subject IDs"
echo "----------------------------------------"

DUPLICATES=$(awk '
NR==1 {
    for (i=1; i<=NF; i++) {
        if ($i == "ID") {
            id_col = i
        }
    }
    next
}

id_col {
    count[$id_col]++
}

END {
    for (id in count) {
        if (count[id] > 1) {
            print id, count[id]
        }
    }
}
' "$COV")

if [ -z "$DUPLICATES" ]; then
    echo "No duplicate IDs detected."
else
    echo "$DUPLICATES"
fi

echo


echo "----------------------------------------"
echo "9. Missing-value indicators"
echo "----------------------------------------"

echo "Counts of common missing-value codes:"

awk '
NR > 1 {
    for (i=1; i<=NF; i++) {
        if ($i == "NA" ||
            $i == "N/A" ||
            $i == "." ||
            $i == "-9") {
            missing[$i]++
        }
    }
}

END {
    if (length(missing) == 0) {
        print "No common missing-value codes detected."
    } else {
        for (value in missing) {
            print value, missing[value]
        }
    }
}
' "$COV"

echo


# Optional sample-ID comparison with the VCF
if [ -n "$VCF" ]; then

    echo "----------------------------------------"
    echo "10. Compare covariate IDs with VCF IDs"
    echo "----------------------------------------"

    if [ ! -f "$VCF" ]; then
        echo "ERROR: VCF file not found: $VCF"
        exit 1
    fi

    if ! command -v bcftools >/dev/null 2>&1; then
        echo "ERROR: bcftools is required for VCF ID comparison."
        exit 1
    fi

    TMP_COV=$(mktemp)
    TMP_VCF=$(mktemp)

    trap 'rm -f "$TMP_COV" "$TMP_VCF"' EXIT

    # Extract covariate IDs
    awk '
    NR==1 {
        for (i=1; i<=NF; i++) {
            if ($i == "ID") {
                id_col = i
            }
        }
        next
    }

    id_col {
        print $id_col
    }
    ' "$COV" | sort -u > "$TMP_COV"

    # Extract VCF IDs
    bcftools query -l "$VCF" | sort -u > "$TMP_VCF"

    COV_IDS=$(wc -l < "$TMP_COV")
    VCF_IDS=$(wc -l < "$TMP_VCF")

    SHARED=$(comm -12 "$TMP_COV" "$TMP_VCF" | wc -l)
    COV_ONLY=$(comm -23 "$TMP_COV" "$TMP_VCF" | wc -l)
    VCF_ONLY=$(comm -13 "$TMP_COV" "$TMP_VCF" | wc -l)

    echo "Unique covariate IDs: $COV_IDS"
    echo "Unique VCF IDs:       $VCF_IDS"
    echo "Shared IDs:           $SHARED"
    echo "Covariate only:       $COV_ONLY"
    echo "VCF only:             $VCF_ONLY"

    echo
fi


echo "========================================"
echo "Inspection complete"
echo "========================================"