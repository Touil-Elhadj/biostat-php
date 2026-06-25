<?php

/**
 * ════════════════════════════════════════════════════════════════════
 * Example 1: descriptive + 2 × 2 tables + Welch t + correlation
 * ────────────────────────────────────────────────────────────────────
 * Run:    php examples/01-descriptive.php
 * ════════════════════════════════════════════════════════════════════
 */

declare(strict_types=1);

// In production:  require_once 'vendor/autoload.php';
require_once __DIR__ . '/../autoload.php';

use TouilElhadj\BiostatPhp\BiostatAnalysis;

$stats = new BiostatAnalysis();

echo "─── Descriptive statistics ───────────────────────────────\n";
$x = [2.3, 4.1, 5.7, 6.2, 7.8, 8.4, 9.1, 10.5, 11.2, 12.9];
printf("  mean       = %.4f\n", $stats->mean($x));
printf("  std        = %.4f\n", $stats->std($x));
printf("  median     = %.4f\n", $stats->median($x));
printf("  Q25 / Q75  = %.4f / %.4f\n",
    $stats->quantile($x, 0.25),
    $stats->quantile($x, 0.75)
);

echo "\n─── 2 × 2 table  (a=45, b=18, c=30, d=22) ────────────────\n";
$chi = $stats->chi2Test2x2(45, 18, 30, 22);
printf("  chi²       = %.4f   p = %.4f\n", $chi['chi2'], $chi['p']);

$or = $stats->oddsRatio(45, 18, 30, 22);
printf("  OR         = %.3f   95%% CI [%.3f, %.3f]\n",
    $or['or'], $or['ci_low'], $or['ci_high']
);

echo "\n─── Welch's t-test ───────────────────────────────────────\n";
$g1 = [10, 12, 11, 13, 14, 12, 11];
$g2 = [15, 17, 16, 18, 16, 17, 15];
$t = $stats->tTest($g1, $g2);
printf("  t          = %.4f   df = %.2f   p = %.4f\n",
    $t['t'], $t['df'], $t['p']
);

echo "\n─── Pearson correlation ──────────────────────────────────\n";
$xx = [10, 12, 14, 16, 18, 20, 22, 24, 26, 28];
$yy = [15, 18, 17, 22, 24, 26, 25, 29, 32, 30];
$pr = $stats->pearson($xx, $yy);
printf("  r          = %.4f   p = %.4f\n", $pr['r'], $pr['p']);

echo "\nDone.\n";
