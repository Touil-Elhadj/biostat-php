<?php

declare(strict_types=1);

namespace TouilElhadj\BiostatPhp\Tests;

/**
 * Tests for the Variance Inflation Factor (multicollinearity diagnostic).
 *
 * Reference values were obtained in R 4.3.0 with car::vif():
 *
 *   # Three uncorrelated predictors → VIF ≈ 1 for each
 *   set.seed(42)
 *   x1 <- rnorm(100); x2 <- rnorm(100); x3 <- rnorm(100)
 *   y  <- 2*x1 + x2 + rnorm(100, 0, 0.5)
 *   car::vif(lm(y ~ x1 + x2 + x3))
 *   # x1: 1.013   x2: 1.006   x3: 1.018
 *
 *   # Highly-correlated predictors
 *   x4 <- x1 + rnorm(100, 0, 0.1)   # x4 ≈ x1
 *   car::vif(lm(y ~ x1 + x2 + x4))
 *   # x1: ~75    x4: ~75    x2: ~1   (very high — strong multicollinearity)
 */
final class VifTest extends BiostatTestCase
{
    public function testVifOnIndependentPredictorsIsNearOne(): void
    {
        // 3 deterministic, near-orthogonal predictors
        $X = [];
        for ($i = 0; $i < 30; $i++) {
            $X[] = [
                cos($i * 0.6),
                sin($i * 0.4),
                cos($i * 0.9) + sin($i * 1.1),
            ];
        }
        $r = $this->stats->vif($X, ['x1', 'x2', 'x3']);
        $this->assertIsArray($r);
        $this->assertArrayNotHasKey('error', $r);
        foreach (['x1', 'x2', 'x3'] as $name) {
            $this->assertLessThan(2.0, (float) $r[$name]['vif'], "$name VIF must be small");
        }
    }

    public function testVifDetectsMulticollinearity(): void
    {
        // x2 ≈ 2 · x1 with small noise → VIF must explode (but matrix still invertible)
        $X = [];
        for ($i = 0; $i < 50; $i++) {
            $x1 = (float) $i / 10.0;
            $x2 = 2 * $x1 + cos($i * 0.3) * 0.05;   // ε = small noise so not exact
            $x3 = sin($i * 0.7);
            $X[] = [$x1, $x2, $x3];
        }
        $r = $this->stats->vif($X, ['x1', 'x2', 'x3']);
        $this->assertGreaterThan(10.0, (float) $r['x1']['vif'], 'collinear x1 should have VIF > 10');
        $this->assertGreaterThan(10.0, (float) $r['x2']['vif'], 'collinear x2 should have VIF > 10');
        $this->assertLessThan(5.0,    (float) $r['x3']['vif'], 'orthogonal x3 should have VIF < 5');
    }
}
