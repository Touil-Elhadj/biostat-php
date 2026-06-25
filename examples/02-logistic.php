<?php

/**
 * ════════════════════════════════════════════════════════════════════
 * Example 2: logistic regression (simple and multivariate)
 *            + Hosmer–Lemeshow goodness-of-fit
 *            + AUC from predictions
 * ────────────────────────────────────────────────────────────────────
 * Run:    php examples/02-logistic.php
 * ════════════════════════════════════════════════════════════════════
 */

declare(strict_types=1);

require_once __DIR__ . '/../autoload.php';

use TouilElhadj\BiostatPhp\BiostatAnalysis;

$stats = new BiostatAnalysis();

// ── Simple logistic: outcome y vs. single continuous predictor x ──
$y = [0,0,0,0,0, 0,0,1,0,1, 0,1,0,1,1, 1,1,0,1,1, 1,0,1,1,1, 1,1,1,1,1];
$x = [1,2,2,3,3, 4,4,4,5,5, 5,6,6,6,7, 7,7,8,8,8, 9,9,9,10,10, 10,11,11,12,12];

echo "─── Simple logistic regression ───────────────────────────\n";
$r = $stats->logisticRegression($y, $x);
printf("  β (slope)  = %.4f\n", $r['coef']);
printf("  OR         = %.3f   95%% CI [%.3f, %.3f]\n",
    $r['or'], $r['ci_low'], $r['ci_high']
);
printf("  p-value    = %.4f\n", $r['p']);

// Build predicted probabilities to run Hosmer-Lemeshow
$b0 = $r['intercept'] ?? 0.0;
$b1 = $r['coef'];
$p  = array_map(
    fn($xi) => 1.0 / (1.0 + exp(-($b0 + $b1 * $xi))),
    $x
);

echo "\n─── Hosmer–Lemeshow goodness-of-fit (5 groups) ───────────\n";
$hl = $stats->hosmerLemeshow($y, $p, 5);
printf("  chi²       = %.3f   df = %d   p = %.4f\n",
    $hl['chi2'], $hl['df'], $hl['p']
);
echo "  (a non-significant p means no evidence of misfit — good)\n";

// ── Multivariate logistic with 2 covariates ─────────────────────────
echo "\n─── Multivariate logistic regression ─────────────────────\n";
$X = [];
foreach ($x as $i => $xi) {
    $X[] = [$xi, $xi * 0.5 + 1.0];   // a roughly redundant 2nd covariate
}
$rm = $stats->logisticRegressionMulti($y, $X, ['x1', 'x2']);
printf("  AUC        = %.4f\n", $rm['auc']);
printf("  converged  = %s   iter = %d\n",
    $rm['converged'] ? 'yes' : 'no',
    $rm['iter']
);

echo "\nDone.\n";
