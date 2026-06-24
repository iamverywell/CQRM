# generate_data.R
# ─────────────────────────────────────────────────────────────────────────────
# Generates a synthetic kidney-inspired dataset for the CQRM simulation study.
#
# Strategy (mirrors synthetic_nonpara_server_1.R):
#   1. Load the kidney RNA-seq dataset from the SimSeq package.
#   2. Filter genes with > 80% zero counts.
#   3. Identify "true DE" genes using a consensus of DESeq2 and Wilcoxon at
#      a stringent threshold (adj p < 0.001 in both).
#   4. Build a synthetic dataset by sampling from the DE / non-DE gene pools:
#      - Sample n_gene * DE_prop genes from the DE pool.
#      - Sample n_gene * (1 - DE_prop) genes from the non-DE pool, then
#        permute their group labels so they are truly null.
#   5. Subsample to `size` samples per group.
#   6. Save the result to data/synthetic_kidney.rds.
#
# Parameters:
#   n_gene   — total number of genes in the synthetic dataset   (default 500)
#   DE_prop  — proportion of true DE genes                      (default 0.10)
#   size     — samples per group                                (default 50)
#   seed     — random seed for reproducibility                  (default 42)
#
# Required packages: SimSeq, DESeq2, edgeR, gtools
# ─────────────────────────────────────────────────────────────────────────────

library(SimSeq)
library(DESeq2)
library(gtools)

# ── Parameters ────────────────────────────────────────────────────────────────
n_gene  <- 500
DE_prop <- 0.10
size    <- 50    # samples per group drawn from kidney (kidney has 72 per group)
seed    <- 42
alpha_consensus <- 0.001   # threshold for consensus DE gene pool

# ── 1. Load and filter kidney data ────────────────────────────────────────────
data(kidney)
counts_full <- kidney$counts
group_full  <- as.numeric(kidney$treatment) - 1   # 0 / 1

# Remove genes where >80% of samples have zero counts
zero_frac <- apply(counts_full == 0, 1, sum) / ncol(counts_full)
counts_full <- counts_full[zero_frac < 0.8, ]
cat(sprintf("Genes after zero filter: %d\n", nrow(counts_full)))

count1_full <- counts_full[, group_full == 0]
count2_full <- counts_full[, group_full == 1]

# ── 2. Identify true DE genes (consensus: DESeq2 + Wilcoxon) ─────────────────
cat("Running DESeq2 on the full kidney dataset...\n")
cond_idx <- factor(group_full)
dds <- DESeqDataSetFromMatrix(counts_full, DataFrame(cond_idx), ~cond_idx)
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds)
padj_deseq2 <- res$padj
padj_deseq2[is.na(padj_deseq2)] <- 1

cat("Running Wilcoxon on DESeq2-normalized counts...\n")
norm_counts <- counts(dds, normalized = TRUE)
p_wilcoxon <- sapply(seq_len(nrow(norm_counts)), function(j) {
  g1 <- log2(norm_counts[j, group_full == 0] + 0.001)
  g2 <- log2(norm_counts[j, group_full == 1] + 0.001)
  if (mean(g1) == mean(g2)) 1 else
    wilcox.test(g1, g2, alternative = "two.sided")$p.value
})
padj_wilcoxon <- p.adjust(p_wilcoxon, method = "BH")

# Consensus: significant in both at the stringent threshold
DE_ind     <- (padj_deseq2 <= alpha_consensus) & (padj_wilcoxon <= alpha_consensus)
non_DE_ind <- !DE_ind
cat(sprintf("Consensus true DE genes: %d  |  non-DE genes: %d\n",
            sum(DE_ind), sum(non_DE_ind)))

# ── 3. Build synthetic dataset ────────────────────────────────────────────────
set.seed(seed)

# Subsample columns to `size` per group
col1 <- sample(seq_len(ncol(count1_full)), size)
col2 <- sample(seq_len(ncol(count2_full)), size)
counts_sub <- cbind(count1_full[, col1], count2_full[, col2])

# Sample DE genes (with replacement if pool is too small)
n_DE     <- round(n_gene * DE_prop)
n_non_DE <- n_gene - n_DE
idx_DE <- sample(which(DE_ind),
                 size    = n_DE,
                 replace = sum(DE_ind) < n_DE)

# Sample non-DE genes and permute group labels within each gene
idx_non_DE <- sample(which(non_DE_ind),
                     size    = n_non_DE,
                     replace = sum(non_DE_ind) < n_non_DE)
perm_group <- permute(c(rep(0, size), rep(1, size)))   # random permutation

counts_DE     <- counts_sub[idx_DE, ]
counts_non_DE <- counts_sub[idx_non_DE, perm_group + 1]  # permuted columns

# Stack: first n_DE rows are true DE, rest are null
syn_counts <- rbind(counts_DE, counts_non_DE)
rownames(syn_counts) <- seq_len(nrow(syn_counts))

# Group vector and true DE indicator
syn_group  <- c(rep(0, size), rep(1, size))
is_true_DE <- c(rep(TRUE, n_DE), rep(FALSE, n_non_DE))

# ── 4. Save ───────────────────────────────────────────────────────────────────
dir.create("data", showWarnings = FALSE)
saveRDS(
  list(
    counts     = syn_counts,
    group      = syn_group,
    is_true_DE = is_true_DE,
    n_DE       = n_DE,
    n_non_DE   = n_non_DE,
    size       = size,
    seed       = seed
  ),
  file = "data/synthetic_kidney.rds"
)
cat("Saved: data/synthetic_kidney.rds\n")
cat(sprintf("  Genes: %d  (DE: %d, non-DE: %d)\n", n_gene, n_DE, n_non_DE))
cat(sprintf("  Samples per group: %d\n", size))
