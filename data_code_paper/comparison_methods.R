# comparison_methods.R
# ─────────────────────────────────────────────────────────────────────────────
# Clean wrappers for all comparison methods used in the CQRM simulation study.
#
# Methods included:
#   mydeseq2           — DESeq2 (adj p-value)
#   my.wilcoxon.deseq2 — Wilcoxon on DESeq2-normalized log counts (BH adjusted)
#   myedger            — Robust edgeR GLM (BH adjusted)
#   my.limmavoom       — limma-voom (adj p-value)
#   my.NOISeq          — NOISeqBIO, TMM normalization (1 - prob, treated as q)
#   my.ZIAQ            — Zero-Inflated Adjusted Quantile regression (p-value)
#   CQRM wrappers      — my.cqrm, my.cqrm.best, my.cqrm.region, my.cqrm.approx
#
# CQRM method notes:
#   All CQRM wrappers precompute solve(subsigma) outside the B-loop for
#   efficiency. The gene loop still calls solve(subsigma) per gene per subset
#   because sigma changes per gene.
#
# Required packages:
#   MASS, lqmm, quantreg — for CQRM methods (or: library(CQRM))
#   DESeq2, edgeR, limma — comparison methods
#   NOISeq               — comparison method
#   quantreg, metap      — for ZIAQ (provide path to ZIAQ.R or install package)
#
# Usage:
#   source("comparison_methods.R")
#   pv <- mydeseq2(count1, count2)
# ─────────────────────────────────────────────────────────────────────────────

# ── Load CQRM (either as installed package or from source) ────────────────────
if (requireNamespace("CQRM", quietly = TRUE)) {
  library(CQRM)
} else {
  # Adjust path if running from a different working directory
  pkg_root <- file.path(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
  source(file.path(pkg_root, "R", "internals.R"))
  source(file.path(pkg_root, "R", "CQRM.R"))
}
library(MASS)
library(lqmm)
library(quantreg)


# ══════════════════════════════════════════════════════════════════════════════
# Comparison method wrappers
# ══════════════════════════════════════════════════════════════════════════════

# ── DESeq2 ────────────────────────────────────────────────────────────────────
# Returns: BH-adjusted p-values (padj); NA replaced by 1.
mydeseq2 <- function(count1, count2) {
  requireNamespace("DESeq2", quietly = TRUE) || stop("DESeq2 not available")
  n1 <- ncol(count1); n2 <- ncol(count2)
  cond_idx <- factor(c(rep(1, n1), rep(2, n2)))
  dat <- cbind(count1, count2)
  rownames(dat) <- seq_len(nrow(dat))
  dds <- DESeq2::DESeqDataSetFromMatrix(dat, S4Vectors::DataFrame(cond_idx), ~cond_idx)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  padj <- DESeq2::results(dds)$padj
  padj[is.na(padj)] <- 1
  padj
}


# ── Wilcoxon on DESeq2-normalized counts ──────────────────────────────────────
# Returns: BH-adjusted Wilcoxon p-values (uses DESeq2 normalization).
my.wilcoxon.deseq2 <- function(count1, count2) {
  requireNamespace("DESeq2", quietly = TRUE) || stop("DESeq2 not available")
  n1 <- ncol(count1); n2 <- ncol(count2)
  cond_idx <- factor(c(rep(1, n1), rep(2, n2)))
  dat <- cbind(count1, count2)
  rownames(dat) <- seq_len(nrow(dat))
  dds <- DESeq2::DESeqDataSetFromMatrix(dat, S4Vectors::DataFrame(cond_idx), ~cond_idx)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  norm <- DESeq2::counts(dds, normalized = TRUE)
  pv <- sapply(seq_len(nrow(norm)), function(j) {
    g1 <- log2(norm[j, cond_idx == 1] + 0.001)
    g2 <- log2(norm[j, cond_idx == 2] + 0.001)
    if (mean(g1) == mean(g2)) 1 else
      wilcox.test(g1, g2, alternative = "two.sided")$p.value
  })
  p.adjust(pv, method = "BH")
}


# ── edgeR (robust GLM) ────────────────────────────────────────────────────────
# Returns: BH-adjusted p-values.
myedger <- function(count1, count2) {
  requireNamespace("edgeR", quietly = TRUE) || stop("edgeR not available")
  count <- cbind(count1, count2)
  cond  <- c(rep(0, ncol(count1)), rep(1, ncol(count2)))
  y      <- edgeR::DGEList(counts = count, group = cond)
  y      <- edgeR::calcNormFactors(y)
  design <- model.matrix(~cond)
  y      <- edgeR::estimateGLMRobustDisp(y, design, prior.df = 20)
  dFit   <- edgeR::glmFit(y, design)
  dLrt   <- edgeR::glmLRT(dFit, coef = 2)
  p.adjust(dLrt$table$PValue, method = "BH")
}


# ── limma-voom ────────────────────────────────────────────────────────────────
# Returns: adj.P.Val from topTable; unfiltered genes get q = 1.
my.limmavoom <- function(count1, count2) {
  requireNamespace("edgeR", quietly = TRUE) || stop("edgeR not available")
  requireNamespace("limma", quietly = TRUE) || stop("limma not available")
  count <- cbind(count1, count2)
  rownames(count) <- seq_len(nrow(count))
  cond  <- c(rep(0, ncol(count1)), rep(1, ncol(count2)))
  design <- model.matrix(~cond)
  dgel   <- edgeR::DGEList(counts = count)
  dgel   <- edgeR::calcNormFactors(dgel)
  keep   <- edgeR::filterByExpr(dgel, design)
  dgel   <- dgel[keep, keep.lib.sizes = FALSE]
  v      <- limma::voom(dgel, design, plot = FALSE)
  fit    <- limma::eBayes(limma::lmFit(v, design))
  de     <- limma::topTable(fit, coef = ncol(design), n = nrow(dgel),
                             sort.by = "none")
  q_all  <- rep(1, nrow(count))
  gene_i <- as.numeric(rownames(de))
  q_all[gene_i] <- replace(de$adj.P.Val, is.na(de$adj.P.Val), 1)
  q_all
}


# ── NOISeq ────────────────────────────────────────────────────────────────────
# Returns: 1 - prob (smaller = more likely DE); NA replaced by 1.
my.NOISeq <- function(count1, count2) {
  requireNamespace("NOISeq", quietly = TRUE) || stop("NOISeq not available")
  count   <- cbind(count1, count2)
  cond    <- c(rep(0, ncol(count1)), rep(1, ncol(count2)))
  factors <- data.frame(Group = cond)
  mydata  <- NOISeq::readData(data = count, factors = factors)
  res     <- NOISeq::noiseqbio(mydata, norm = "tmm", factor = "Group", filter = 0)
  q       <- 1 - as.data.frame(res@results[[1]])$prob
  q[is.na(q)] <- 1
  q
}


# ── ZIAQ ──────────────────────────────────────────────────────────────────────
# Zero-Inflation Adjusted Quantile regression.
# Requires the ZIAQ.R helper (available at the repo for the original paper).
# Returns: raw p-values (apply p.adjust for FDR control).
#
# ziaq_path: path to ZIAQ.R if not already sourced
my.ZIAQ <- function(count1, count2,
                    tau_vec   = seq(0.1, 0.9, 0.1),
                    ziaq_path = NULL) {
  if (!exists("ziaq", mode = "function")) {
    if (is.null(ziaq_path)) stop("ZIAQ not loaded. Supply ziaq_path= to source ZIAQ.R")
    source(ziaq_path)
  }
  requireNamespace("quantreg", quietly = TRUE) || stop("quantreg not available")
  requireNamespace("metap",    quietly = TRUE) || stop("metap not available")
  n1 <- ncol(count1); n2 <- ncol(count2)
  group  <- factor(c(rep(1, n1), rep(2, n2)))
  dat    <- cbind(count1, count2)
  rownames(dat) <- seq_len(nrow(dat))
  colDat <- data.frame(group = group)
  res    <- ziaq(dat, colDat, formula = ~group, group = "group",
                 probs = tau_vec, log_i = TRUE, parallel = FALSE, no.core = 1)
  res$pvalue
}


# ══════════════════════════════════════════════════════════════════════════════
# CQRM method wrappers
# (solve(subsigma) precomputed outside the B-loop in all three matrix methods)
# ══════════════════════════════════════════════════════════════════════════════

# ── CQRM (all 2^K-1 subsets) ─────────────────────────────────────────────────
# Returns: calibrated p-values (apply p.adjust for FDR control).
#
# Parameters:
#   count1, count2 — gene × sample count matrices
#   tau_vec        — quantile levels, default c(0.2, 0.4, 0.6, 0.8)
#   B              — null draws for calibration, default 10000
#   seed           — random seed, default 12345
#   M              — jitter replicates per gene per quantile, default 50
my.cqrm <- function(count1, count2,
                    tau_vec = c(0.2, 0.4, 0.6, 0.8),
                    B = 10000, seed = 12345, M = 50) {
  CQRM(count1, count2, tau_vec = tau_vec, B = B, seed = seed, M = M)
}


# ── CQRM-Best (single-quantile MinP) ─────────────────────────────────────────
# Returns: calibrated p-values (apply p.adjust for FDR control).
#
# Parameters: same as my.cqrm; no matrix solve — scalar variance per quantile.
my.cqrm.best <- function(count1, count2,
                         tau_vec = c(0.2, 0.4, 0.6, 0.8),
                         B = 10000, seed = 12345, M = 50) {
  CQRM.best(count1, count2, tau_vec = tau_vec, B = B, seed = seed, M = M)
}


# ── CQRM-Region (contiguous-window MinP) ─────────────────────────────────────
# Returns: calibrated p-values (apply p.adjust for FDR control).
#
# Parameters:
#   k — window width (consecutive quantiles per test), default 2
my.cqrm.region <- function(count1, count2,
                            tau_vec = c(0.2, 0.4, 0.6, 0.8), k = 2,
                            B = 10000, seed = 12345, M = 50) {
  CQRM.region(count1, count2, tau_vec = tau_vec, k = k,
               B = B, seed = seed, M = M)
}


# ── CQRM-Approx (fast, no jittering) ─────────────────────────────────────────
# Returns: calibrated p-values (apply p.adjust for FDR control).
#
# Parameters:
#   zeta — floor for log transformation, default 1e-5; no M argument
my.cqrm.approx <- function(count1, count2,
                            tau_vec = c(0.2, 0.4, 0.6, 0.8),
                            B = 10000, seed = 12345, zeta = 1e-5) {
  CQRM.approx(count1, count2, tau_vec = tau_vec, B = B,
               seed = seed, zeta = zeta)
}
