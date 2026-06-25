# ════════════════════════════════════════════════════════════════════
# reference-values.R
# ────────────────────────────────────────────────────────────────────
# Reproduces every R reference value used in docs/validation-tables.md.
#
# A JOSS reviewer (or any user) can run this script independently to
# verify the numerical claims of the biostat-php test suite.
#
# Run:
#   Rscript tests/fixtures/reference-values.R
#
# Required packages: stats (base), epitools, car
# Optional packages: lme4, geepack, mice, ResourceSelection, pROC
# ════════════════════════════════════════════════════════════════════

options(scipen = 999, digits = 7)

cat("┌──────────────────────────────────────────────────────────────┐\n")
cat("│  biostat-php — R reference values for docs/validation-tables │\n")
cat("│  ", R.version.string, "\n", sep = "")
cat("└──────────────────────────────────────────────────────────────┘\n\n")

set.seed(42)

# ════════════════════════════════════════════════════════════════════
# Table 1 — Descriptive statistics
# ════════════════════════════════════════════════════════════════════
cat("── Table 1 — Descriptive ───────────────────────────────────────\n")
x <- c(2.3, 4.1, 5.7, 6.2, 7.8, 8.4, 9.1, 10.5, 11.2, 12.9)
cat(sprintf("  mean(x)               = %.6f\n", mean(x)))
cat(sprintf("  sd(x)                 = %.6f\n", sd(x)))
cat(sprintf("  median(x)             = %.6f\n", median(x)))
cat(sprintf("  quantile(x, 0.25)     = %.6f\n",
            quantile(x, 0.25, type = 7, names = FALSE)))
cat(sprintf("  quantile(x, 0.75)     = %.6f\n\n",
            quantile(x, 0.75, type = 7, names = FALSE)))

# ════════════════════════════════════════════════════════════════════
# Table 2 — 2 x 2 contingency
# ════════════════════════════════════════════════════════════════════
cat("── Table 2 — 2 x 2 contingency ────────────────────────────────\n")
m <- matrix(c(45, 18, 30, 22), nrow = 2, byrow = TRUE)
ct  <- chisq.test(m, correct = FALSE)
ctY <- chisq.test(m, correct = TRUE)
cat(sprintf("  chisq.test (no Yates) : X2 = %.4f, p = %.4f\n",
            ct$statistic, ct$p.value))
cat(sprintf("  chisq.test (Yates)    : X2 = %.4f, p = %.4f\n",
            ctY$statistic, ctY$p.value))
if (requireNamespace("epitools", quietly = TRUE)) {
  o <- epitools::oddsratio.wald(m)$measure
  cat(sprintf("  OR (Wald)             : %.4f [%.4f, %.4f]\n",
              o[2, 1], o[2, 2], o[2, 3]))
} else {
  cat("  (epitools not installed — skipping OR)\n")
}
# Zero-cell with Haldane-Anscombe
m0  <- matrix(c(10, 0, 5, 10), nrow = 2, byrow = TRUE)
m0c <- m0 + 0.5
or0 <- (m0c[1, 1] * m0c[2, 2]) / (m0c[1, 2] * m0c[2, 1])
se0 <- sqrt(sum(1 / m0c))
ci0 <- exp(log(or0) + c(-1, 1) * 1.96 * se0)
cat(sprintf("  OR (Haldane-Anscombe) : %.4f [%.4f, %.4f]   (manual)\n\n",
            or0, ci0[1], ci0[2]))

# ════════════════════════════════════════════════════════════════════
# Table 3 — Welch t-test
# ════════════════════════════════════════════════════════════════════
cat("── Table 3 — Welch t-test ─────────────────────────────────────\n")
a <- c(10, 12, 11, 13, 14, 12, 11)
b <- c(15, 17, 16, 18, 16, 17, 15)
tt <- t.test(a, b)
cat(sprintf("  t.test : t = %.4f, df = %.4f, p = %.3e\n\n",
            tt$statistic, tt$parameter, tt$p.value))

# ════════════════════════════════════════════════════════════════════
# Table 4 — One-way ANOVA
# ════════════════════════════════════════════════════════════════════
cat("── Table 4 — One-way ANOVA ────────────────────────────────────\n")
g1 <- c(10, 12, 11); g2 <- c(15, 17, 16); g3 <- c(20, 22, 21)
d <- stack(list(g1 = g1, g2 = g2, g3 = g3))
av <- summary(aov(values ~ ind, data = d))[[1]]
cat(sprintf("  aov   : F(%d, %d) = %.4f, p = %.3e\n\n",
            av[1, "Df"], av[2, "Df"], av[1, "F value"], av[1, "Pr(>F)"]))

# ════════════════════════════════════════════════════════════════════
# Table 5 — Correlation
# ════════════════════════════════════════════════════════════════════
cat("── Table 5 — Correlation ──────────────────────────────────────\n")
x <- c(10, 12, 14, 16, 18, 20, 22, 24, 26, 28)
y <- c(15, 18, 17, 22, 24, 26, 25, 29, 32, 30)
pr <- cor.test(x, y, method = "pearson")
sp <- suppressWarnings(cor.test(x, y, method = "spearman"))
cat(sprintf("  pearson  : r   = %.4f, p = %.3e\n", pr$estimate, pr$p.value))
cat(sprintf("  spearman : rho = %.4f, p = %.3e\n\n", sp$estimate, sp$p.value))

# ════════════════════════════════════════════════════════════════════
# Table 6 — Benjamini-Hochberg
# ════════════════════════════════════════════════════════════════════
cat("── Table 6 — Benjamini-Hochberg FDR ───────────────────────────\n")
p <- c(0.001, 0.008, 0.039, 0.041, 0.042, 0.060, 0.074, 0.205, 0.212, 0.216)
adj <- p.adjust(p, method = "BH")
for (i in seq_along(p)) {
  cat(sprintf("  raw = %.4f   adj = %.4f\n", p[i], adj[i]))
}
cat("\n")

# ════════════════════════════════════════════════════════════════════
# Table 7 — Logistic regression
# ════════════════════════════════════════════════════════════════════
cat("── Table 7 — Logistic regression ──────────────────────────────\n")
y_l <- c(0,0,0,0,0, 0,0,1,0,1, 0,1,0,1,1, 1,1,0,1,1, 1,0,1,1,1, 1,1,1,1,1)
x_l <- c(1,2,2,3,3, 4,4,4,5,5, 5,6,6,6,7, 7,7,8,8,8, 9,9,9,10,10,
         10,11,11,12,12)
mod <- glm(y_l ~ x_l, family = binomial(link = "logit"))
co  <- summary(mod)$coefficients
cat(sprintf("  (Intercept) = %.4f   SE = %.4f   p = %.4f\n",
            co["(Intercept)", 1], co["(Intercept)", 2], co["(Intercept)", 4]))
cat(sprintf("  x           = %.4f   SE = %.4f   p = %.4f\n",
            co["x_l", 1], co["x_l", 2], co["x_l", 4]))
cat(sprintf("  OR(x)       = exp(b) = %.4f\n", exp(co["x_l", 1])))
# AUC via Mann-Whitney
phat <- fitted(mod)
n1 <- sum(y_l == 1); n0 <- sum(y_l == 0)
auc <- (sum(outer(phat[y_l == 1], phat[y_l == 0], ">")) +
        0.5 * sum(outer(phat[y_l == 1], phat[y_l == 0], "=="))) / (n1 * n0)
cat(sprintf("  AUC         = %.4f\n", auc))
# Hosmer-Lemeshow (g=5)
if (requireNamespace("ResourceSelection", quietly = TRUE)) {
  hl <- ResourceSelection::hoslem.test(y_l, phat, g = 5)
  cat(sprintf("  HL (g=5)    : X2 = %.4f, df = %d, p = %.4f\n\n",
              hl$statistic, hl$parameter, hl$p.value))
} else {
  cat("\n")
}

# ════════════════════════════════════════════════════════════════════
# Table 8 — VIF
# ════════════════════════════════════════════════════════════════════
cat("── Table 8 — VIF (car::vif) ───────────────────────────────────\n")
if (requireNamespace("car", quietly = TRUE)) {
  set.seed(42)
  x1 <- rnorm(100); x2 <- rnorm(100); x3 <- rnorm(100)
  y_v <- 2 * x1 + x2 + rnorm(100, 0, 0.5)
  v_ok <- car::vif(lm(y_v ~ x1 + x2 + x3))
  cat(sprintf("  uncorrelated : x1 = %.3f   x2 = %.3f   x3 = %.3f\n",
              v_ok[1], v_ok[2], v_ok[3]))
  x4 <- x1 + rnorm(100, 0, 0.05)
  v_bad <- car::vif(lm(y_v ~ x1 + x4 + x2))
  cat(sprintf("  collinear    : x1 = %.1f    x4 = %.1f    x2 = %.3f\n\n",
              v_bad[1], v_bad[2], v_bad[3]))
} else {
  cat("  (car not installed — skipping)\n\n")
}

# ════════════════════════════════════════════════════════════════════
# Optional — GLMM (lme4) and GEE (geepack)
# ════════════════════════════════════════════════════════════════════
cat("── Optional — GLMM (lme4::glmer) ──────────────────────────────\n")
if (requireNamespace("lme4", quietly = TRUE)) {
  set.seed(17)
  n_c <- 30; n_per <- 8
  cluster <- rep(1:n_c, each = n_per)
  u <- rnorm(n_c, 0, 0.7)
  x <- runif(n_c * n_per, -2, 2)
  eta <- -0.6 + 0.5 * x + u[cluster]
  p <- 1 / (1 + exp(-eta))
  y <- rbinom(length(p), 1, p)
  d <- data.frame(y, x, cluster = factor(cluster))
  m <- suppressMessages(suppressWarnings(
        lme4::glmer(y ~ x + (1 | cluster), family = binomial, data = d)))
  s <- summary(m)
  cat(sprintf("  beta(x)     = %.4f   SE = %.4f\n",
              s$coefficients["x", 1], s$coefficients["x", 2]))
  s2u <- as.numeric(s$varcor$cluster)
  cat(sprintf("  sigma2_u    = %.4f\n", s2u))
  cat(sprintf("  ICC         = %.4f\n\n", s2u / (s2u + pi^2 / 3)))

  if (requireNamespace("geepack", quietly = TRUE)) {
    cat("── Optional — GEE (geepack::geeglm) ───────────────────────────\n")
    m_gee <- geepack::geeglm(y ~ x, id = cluster, data = d,
                             family = binomial, corstr = "exchangeable")
    s_gee <- summary(m_gee)
    cat(sprintf("  beta(x)     = %.4f   robust SE = %.4f\n",
                s_gee$coefficients["x", "Estimate"],
                s_gee$coefficients["x", "Std.err"]))
    cat(sprintf("  alpha       = %.4f\n\n",
                s_gee$corr["Estimate", "alpha"]))
  }
} else {
  cat("  (lme4 not installed — skipping GLMM and GEE)\n\n")
}

cat("─────────────────────────────────────────────────────────────────\n")
cat("All reference values printed.\n")
cat("Compare to docs/validation-tables.md and tests/*.php\n")
