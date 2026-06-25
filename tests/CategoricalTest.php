<?php

declare(strict_types=1);

namespace TouilElhadj\BiostatPhp\Tests;

/**
 * Tests for 2 × 2 contingency-table methods.
 *
 * Reference values were obtained in R 4.3.0 / SciPy 1.11:
 *
 *   m <- matrix(c(45, 18, 30, 22), nrow = 2, byrow = TRUE)
 *   chisq.test(m, correct = FALSE)        # X² = 2.3695, p = 0.1237
 *   chisq.test(m, correct = TRUE)         # Yates: X² = 1.8027, p = 0.1794
 *   epitools::oddsratio.wald(m)$measure   # OR = 1.833 [0.844, 3.984]
 *   fisher.test(m)                         # p = 0.1789
 *
 * Reproducer: see tests/fixtures/reference-values.R
 */
final class CategoricalTest extends BiostatTestCase
{
    public function testChi2WithoutYates(): void
    {
        $r = $this->stats->chi2Test2x2(45, 18, 30, 22, false);
        $this->assertNear(2.3695, (float) $r['chi2'], 0.001, 'chi2 (no Yates)');
        $this->assertNear(0.1237, (float) $r['p'],    0.001, 'p-value (no Yates)');
    }

    public function testChi2WithYates(): void
    {
        $r = $this->stats->chi2Test2x2(45, 18, 30, 22, true);
        $this->assertNear(1.8027, (float) $r['chi2'], 0.001, 'chi2 Yates');
        $this->assertNear(0.1794, (float) $r['p'],    0.001, 'p-value Yates');
    }

    public function testOddsRatioBasic(): void
    {
        $r = $this->stats->oddsRatio(45, 18, 30, 22);
        $this->assertNear(1.833, (float) $r['or'],      0.01, 'OR');
        $this->assertNear(0.844, (float) $r['ci_low'],  0.01, 'CI low');
        $this->assertNear(3.984, (float) $r['ci_high'], 0.01, 'CI high');
    }

    public function testOddsRatioHaldaneAnscombeWithZeroCell(): void
    {
        // Cell b = 0 → continuity correction (+0.5 to every cell)
        // After correction: (10.5, 0.5, 5.5, 10.5) → OR = (10.5 × 10.5)/(0.5 × 5.5) = 40.09
        $r = $this->stats->oddsRatio(10, 0, 5, 10);
        $this->assertIsArray($r);
        $this->assertGreaterThan(10.0, (float) $r['or'], 'OR should be large after correction');
        $this->assertNotNull($r['ci_low']);
        $this->assertNotNull($r['ci_high']);
    }
}
