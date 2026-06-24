# run_simulation.R
# ─────────────────────────────────────────────────────────────────────────────
# Runs one complete simulation replicate on the synthetic kidney-inspired
# dataset and saves results in a displayable format.
#
# Workflow:
#   1. Load the pre-generated synthetic dataset (data/synthetic_kidney.rds).
#      Run generate_data.R first if the file does not exist.
#   2. Extract count matrices for group 1 and group 2.
#   3. Run all comparison methods.
#   4. Compute BH-adjusted p-values (adj_p_val) for each method.
#   5. Compute FDP (false discovery proportion) and power at threshold q = 0.05.
#   6. Save:
#      - results/adj_pval_table.csv   — gene × method matrix of adj p-values
#      - results/fdp_power_table.csv  — summary FDP / power table
# ─────────────────────────────────────────────────────────────────────────────

# ── 0. Housekeeping ───────────────────────────────────────────────────────────
# Set working directory to data_code_paper/ before sourcing this script, or
# adjust the paths below.

source("comparison_methods.R")   # loads CQRM + all wrapper functions
dir.create("results", showWarnings = FALSE)

# ── Parameters ────────────────────────────────────────────────────────────────
# tau_vec  : quantile levels for all CQRM methods
# B        : null draws for permutation calibration
#            (use B = 10000 for production; B = 500 for a fast demo run)
# M        : jitter replicates per gene per quantile (CQRM / CQRM.best / CQRM.region)
# q_thresh : FDR threshold for FDP / power summary
tau_vec  <- c(0.2, 0.4, 0.6, 0.8)
B        <- 1000    # increase to 10000 for paper-quality results
M        <- 20      # increase to 50  for paper-quality results
q_thresh <- 0.05

# ── 1. Load synthetic data ────────────────────────────────────────────────────
if (!file.exists("data/synthetic_kidney.rds")) {
  cat("Synthetic data not found. Running generate_data.R...\n")
  source("generate_data.R")
}
dat        <- readRDS("data/synthetic_kidney.rds")
counts     <- dat$counts          # G × (2 * size) count matrix
group      <- dat$group           # 0 / 1 group vector
is_true_DE <- dat$is_true_DE      # logical vector, length G
size       <- dat$size            # samples per group
G          <- nrow(counts)
n_DE       <- sum(is_true_DE)

count1 <- counts[, group == 0]
count2 <- counts[, group == 1]

cat(sprintf("Loaded synthetic data: %d genes (%d DE, %d non-DE), %d samples/group\n",
            G, n_DE, G - n_DE, size))

# ── 2. Run all methods ────────────────────────────────────────────────────────
# Each function returns a numeric vector of length G.
# For methods that already return adj p-values (DESeq2, edgeR, limma, NOISeq)
# we store them directly. For methods returning raw p-values (Wilcoxon,
# ZIAQ, CQRM variants) we apply BH adjustment.

cat("Running DESeq2...\n")
pv_deseq2   <- mydeseq2(count1, count2)                # already adj

cat("Running Wilcoxon...\n")
pv_wilcoxon <- my.wilcoxon.deseq2(count1, count2)       # already adj (BH inside)

cat("Running edgeR...\n")
pv_edger    <- myedger(count1, count2)                   # already adj

cat("Running limma-voom...\n")
pv_limma    <- tryCatch(my.limmavoom(count1, count2),
                         error = function(e) { message(e); rep(NA, G) })  # adj

cat("Running NOISeq...\n")
pv_noiseq   <- tryCatch(my.NOISeq(count1, count2),
                         error = function(e) { message(e); rep(NA, G) })  # prob-based

cat("Running CQRM.approx...\n")
pv_approx   <- my.cqrm.approx(count1, count2, tau_vec = tau_vec, B = B)

cat("Running CQRM.best...\n")
pv_best     <- my.cqrm.best(count1, count2, tau_vec = tau_vec, B = B, M = M)

cat("Running CQRM.region...\n")
pv_region   <- my.cqrm.region(count1, count2, tau_vec = tau_vec, k = 2,
                               B = B, M = M)

cat("Running CQRM...\n")
pv_cqrm     <- my.cqrm(count1, count2, tau_vec = tau_vec, B = B, M = M)

# ── 3. BH-adjust raw p-values ─────────────────────────────────────────────────
adj_approx <- p.adjust(pv_approx, method = "BH")
adj_best   <- p.adjust(pv_best,   method = "BH")
adj_region <- p.adjust(pv_region, method = "BH")
adj_cqrm   <- p.adjust(pv_cqrm,  method = "BH")

# ── 4. Build adj p-value table (gene × method) ───────────────────────────────
adj_pval_table <- data.frame(
  Gene       = seq_len(G),
  True_DE    = is_true_DE,
  DESeq2     = pv_deseq2,
  Wilcoxon   = pv_wilcoxon,
  edgeR      = pv_edger,
  limma_voom = pv_limma,
  NOISeq     = pv_noiseq,
  CQRM_Approx  = adj_approx,
  CQRM_Best    = adj_best,
  CQRM_Region  = adj_region,
  CQRM         = adj_cqrm
)

write.csv(adj_pval_table, "results/adj_pval_table.csv", row.names = FALSE)
cat("Saved: results/adj_pval_table.csv\n")

# ── 5. FDP / Power summary at q_thresh ────────────────────────────────────────
compute_fdp_power <- function(adj_pv, true_de, q = 0.05) {
  adj_pv[is.na(adj_pv)] <- 1
  disc   <- which(adj_pv <= q)
  if (length(disc) == 0) return(c(FDP = NA, Power = 0))
  fdp   <- sum(!disc %in% which(true_de)) / length(disc)
  power <- sum(disc %in% which(true_de))  / sum(true_de)
  c(FDP = round(fdp, 4), Power = round(power, 4))
}

method_names <- c("DESeq2", "Wilcoxon", "edgeR", "limma-voom", "NOISeq",
                  "CQRM-Approx", "CQRM-Best", "CQRM-Region", "CQRM")
adj_list <- list(pv_deseq2, pv_wilcoxon, pv_edger, pv_limma, pv_noiseq,
                 adj_approx, adj_best, adj_region, adj_cqrm)

fdp_power <- do.call(rbind, lapply(seq_along(method_names), function(i) {
  res <- compute_fdp_power(adj_list[[i]], is_true_DE, q = q_thresh)
  data.frame(Method = method_names[i], FDP = res["FDP"], Power = res["Power"])
}))
rownames(fdp_power) <- NULL

write.csv(fdp_power, "results/fdp_power_table.csv", row.names = FALSE)
cat("Saved: results/fdp_power_table.csv\n")

# ── 6. Print summary tables ───────────────────────────────────────────────────
cat("\n=== FDP / Power at q =", q_thresh, "===\n")
print(fdp_power, digits = 4)

cat("\n=== Top 10 genes by CQRM adj p-value ===\n")
top_genes <- head(adj_pval_table[order(adj_pval_table$CQRM), ], 10)
print(top_genes[, c("Gene", "True_DE", "DESeq2", "Wilcoxon", "edgeR",
                    "CQRM_Approx", "CQRM_Best", "CQRM_Region", "CQRM")],
      digits = 4, row.names = FALSE)
