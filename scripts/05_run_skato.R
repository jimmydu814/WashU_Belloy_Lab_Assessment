#!/usr/bin/env Rscript

# ============================================================
# TREM2 SKAT-O analysis
#
# Usage:
# Rscript scripts/05_run_skato.R \
#     data/covariates.txt \
#     results/skat/trem2_lof_genotypes.raw \
#     results/skat/trem2_missense_genotypes.raw \
#     results/skat
# ============================================================

suppressPackageStartupMessages(library(SKAT))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
    stop(
        "Usage: Rscript 05_run_skato.R ",
        "<covariates.txt> <lof.raw> <missense.raw> <output_dir>"
    )
}

cov_file <- args[1]
lof_file <- args[2]
missense_file <- args[3]
outdir <- args[4]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 1. Read covariates
# ============================================================

cov <- read.table(
    cov_file,
    header = TRUE,
    stringsAsFactors = FALSE,
    colClasses = "character"
)

# Preserve IDs as characters, but convert other columns
# to numeric when appropriate.
for (name in setdiff(names(cov), "ID")) {
    cov[[name]] <- type.convert(
        cov[[name]],
        as.is = TRUE
    )
}

if (!all(c("ID", "DX") %in% names(cov))) {
    stop("Covariate file must contain ID and DX columns.")
}

cov$ID <- as.character(cov$ID)

# Confirm phenotype coding
if (!all(na.omit(unique(cov$DX)) %in% c(0, 1))) {
    stop("DX must be coded as 0 = control and 1 = case.")
}

# Check duplicate IDs
if (anyDuplicated(cov$ID)) {
    stop("Duplicate subject IDs detected in covariate file.")
}


# ============================================================
# 2. Identify covariates
#
# Use ALL supplied columns except ID and DX.
# ============================================================

covariate_cols <- setdiff(names(cov), c("ID", "DX"))

cat("Covariates included in null model:\n")
print(covariate_cols)

# Convert nonnumeric categorical covariates to factors
for (name in covariate_cols) {
    if (is.character(cov[[name]])) {
        cov[[name]] <- factor(cov[[name]])
    }
}


# ============================================================
# 3. Read PLINK .raw genotype file
# ============================================================

read_genotypes <- function(filename) {

    x <- read.table(
        filename,
        header = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE,
        comment.char = "",
        colClasses = "character"
    )

    if (!"IID" %in% names(x)) {
        stop(paste("IID column missing from", filename))
    }

    x$IID <- as.character(x$IID)

    metadata_cols <- c(
        "FID", "#FID", "IID",
        "PAT", "MAT", "SEX", "PHENOTYPE"
    )

    variant_cols <- setdiff(
        names(x),
        metadata_cols
    )

    if (length(variant_cols) == 0) {
        stop(paste("No genotype columns found in", filename))
    }

    list(
        data = x,
        variants = variant_cols
    )
}


lof_raw <- read_genotypes(lof_file)
missense_raw <- read_genotypes(missense_file)


# ============================================================
# 4. Determine common subjects
# ============================================================

common_ids <- Reduce(
    intersect,
    list(
        cov$ID,
        lof_raw$data$IID,
        missense_raw$data$IID
    )
)

cat("\nSubjects with genotype + covariate data:",
    length(common_ids), "\n")


# Reorder covariates
analysis_cov <- cov[
    match(common_ids, cov$ID),
    ,
    drop = FALSE
]


# Remove subjects missing phenotype or any covariate
model_cols <- c("DX", covariate_cols)

complete <- complete.cases(
    analysis_cov[, model_cols, drop = FALSE]
)

analysis_cov <- analysis_cov[complete, , drop = FALSE]

analysis_ids <- analysis_cov$ID

cat("Subjects after removing missing phenotype/covariates:",
    length(analysis_ids), "\n")


# ============================================================
# 5. Build genotype matrices
# ============================================================

make_Z <- function(raw_object, ids) {

    x <- raw_object$data
    variant_cols <- raw_object$variants

    idx <- match(ids, x$IID)

    if (any(is.na(idx))) {
        stop("Some analysis subjects are missing from genotype data.")
    }

    Z <- sapply(
        x[idx, variant_cols, drop = FALSE],
        function(v) suppressWarnings(as.numeric(v))
    )

    # Ensure matrix structure even if only one variant exists
    if (is.null(dim(Z))) {
        Z <- matrix(Z, ncol = 1)
        colnames(Z) <- variant_cols
    }

    rownames(Z) <- ids

    # Validate genotype values
    valid <- is.na(Z) | Z %in% c(0, 1, 2)

    if (!all(valid)) {
        stop("Unexpected genotype values detected.")
    }


    # --------------------------------------------------------
    # PLINK 2 --export A counts REF alleles by default.
    #
    # SKAT expects 0/1/2 to represent minor-allele counts.
    # Flip variants when the currently counted allele has
    # frequency > 0.5.
    # --------------------------------------------------------

    counted_freq <- colMeans(Z, na.rm = TRUE) / 2
    flip <- counted_freq > 0.5

    if (any(flip)) {
        Z[, flip] <- 2 - Z[, flip]
    }

    cat(
        "Variants converted to minor-allele coding:",
        sum(flip), "of", ncol(Z), "\n"
    )

    Z
}


Z_lof <- make_Z(lof_raw, analysis_ids)
Z_missense <- make_Z(missense_raw, analysis_ids)

cat("\nLoF variants:", ncol(Z_lof), "\n")
cat("Missense variants:", ncol(Z_missense), "\n")


# ============================================================
# 6. Build SKAT null model
#
# DX ~ all supplied covariates
# ============================================================

if (length(covariate_cols) > 0) {

    null_formula <- reformulate(
        covariate_cols,
        response = "DX"
    )

} else {

    null_formula <- DX ~ 1
}

cat("\nSKAT null model:\n")
print(null_formula)


null_model <- SKAT_Null_Model(
    null_formula,
    data = analysis_cov,
    out_type = "D"
)


# ============================================================
# 7. Run SKAT-O
#
# method      = SKATO
# method.bin  = Hybrid
# Beta(1,25)  = standard rare-variant weighting
# ============================================================

cat("\nRunning LoF SKAT-O...\n")

lof_result <- SKATBinary(
    Z_lof,
    null_model,
    method = "SKATO",
    method.bin = "Hybrid",
    weights.beta = c(1, 25),
    impute.method = "bestguess"
)


cat("Running missense SKAT-O...\n")

missense_result <- SKATBinary(
    Z_missense,
    null_model,
    method = "SKATO",
    method.bin = "Hybrid",
    weights.beta = c(1, 25),
    impute.method = "bestguess"
)


# ============================================================
# 8. Create results table
# ============================================================

get_tested_variants <- function(result, fallback) {

    if (!is.null(result$param$n.marker.test)) {
        return(result$param$n.marker.test)
    }

    fallback
}


results <- data.frame(
    Category = c("LoF", "Missense"),

    Input_variants = c(
        ncol(Z_lof),
        ncol(Z_missense)
    ),

    Tested_variants = c(
        get_tested_variants(lof_result, ncol(Z_lof)),
        get_tested_variants(missense_result, ncol(Z_missense))
    ),

    P_value = c(
        lof_result$p.value,
        missense_result$p.value
    ),

    Total_MAC = c(
        if (!is.null(lof_result$MAC))
            sum(lof_result$MAC) else NA,
        if (!is.null(missense_result$MAC))
            sum(missense_result$MAC) else NA
    ),

    Carriers = c(
        if (!is.null(lof_result$m))
            lof_result$m else NA,
        if (!is.null(missense_result$m))
            missense_result$m else NA
    )
)


# ============================================================
# 9. Save results
# ============================================================

write.table(
    results,
    file = file.path(outdir, "skato_results.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# Save the exact covariate model
capture.output(
    null_formula,
    file = file.path(outdir, "skato_model.txt")
)


# Save R/package versions for reproducibility
capture.output(
    sessionInfo(),
    file = file.path(outdir, "R_session_info.txt")
)


# ============================================================
# 10. Print results
# ============================================================

cat("\n========================================\n")
cat("SKAT-O RESULTS\n")
cat("========================================\n\n")

print(results)

cat("\nModel:\n")
print(null_formula)

cat("\nResults saved to:\n")
cat(file.path(outdir, "skato_results.tsv"), "\n")