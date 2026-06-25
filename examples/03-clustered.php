<?php

/**
 * ════════════════════════════════════════════════════════════════════
 * Example 3: clustered design — GLMM (PQL) and GEE (sandwich)
 * ────────────────────────────────────────────────────────────────────
 * Run:    php examples/03-clustered.php
 *
 * Simulates a small clustered binary outcome where the within-cluster
 * dependence is non-negligible, and compares the cluster-naive logistic
 * regression to GLMM and GEE.
 * ════════════════════════════════════════════════════════════════════
 */

declare(strict_types=1);

require_once __DIR__ . '/../autoload.php';

use TouilElhadj\BiostatPhp\BiostatAnalysis;

$stats = new BiostatAnalysis();

// ── Build a deterministic clustered dataset with random-like noise ──
// 30 clusters × 8 observations, true β_x ≈ 0.5, with random-intercepts
// drawn from a pseudo-Gaussian deterministic sequence so the example
// is reproducible without seed handling.
$n_clusters = 30;
$cluster = [];
$y       = [];
$X       = [];
$seed    = 17;

for ($c = 0; $c < $n_clusters; $c++) {
    // Deterministic pseudo-Gaussian via Box–Muller on (sin, cos)
    $u = cos($c * 1.31 + $seed) * 0.7;
    for ($j = 0; $j < 8; $j++) {
        $x   = ($j - 3.5) * 0.4 + sin($c * 0.9 + $j * 2.1) * 0.2;
        $eta = -0.6 + 0.5 * $x + $u;
        $p   = 1.0 / (1.0 + exp(-$eta));
        // Pseudo-Bernoulli: deterministic hash-based threshold
        $rand = (sin(($c + 1) * 911.0 + ($j + 1) * 137.0) + 1.0) / 2.0;
        $cluster[] = $c;
        $X[]       = [$x];
        $y[]       = $rand < $p ? 1 : 0;
    }
}

echo "─── Naive logistic (ignoring clustering) ─────────────────\n";
$xflat = array_map(fn($row) => $row[0], $X);
$naive = $stats->logisticRegression($y, $xflat);
printf("  β        = %.4f   p = %.4f\n", $naive['coef'], $naive['p']);

echo "\n─── GLMM (logistic, random intercept, PQL) ───────────────\n";
$glmm = $stats->glmmLogistic($y, $X, $cluster, ['x'], 50);
if (isset($glmm['error'])) {
    echo "  error: {$glmm['error']}\n";
} else {
    printf("  β (x)    = %.4f\n", $glmm['coef'][0] ?? $glmm['coef']);
    printf("  σ²_u     = %.4f\n", $glmm['sigma2_u'] ?? 0);
    printf("  ICC      = %.4f\n", $glmm['icc']);
    printf("  converged = %s\n", $glmm['converged'] ? 'yes' : 'no');
}

echo "\n─── GEE (exchangeable, Liang–Zeger sandwich) ─────────────\n";
$gee = $stats->geeLogistic($y, $X, $cluster, ['x']);
if (isset($gee['error'])) {
    echo "  error: {$gee['error']}\n";
} else {
    printf("  β (x)        = %.4f\n", $gee['coef'][0] ?? $gee['coef']);
    printf("  SE model     = %.4f\n", $gee['model_se'][0] ?? 0);
    printf("  SE robust    = %.4f\n", $gee['robust_se'][0] ?? 0);
    printf("  α (working)  = %.4f\n", $gee['alpha'] ?? 0);
}

echo "\nDone.\n";
