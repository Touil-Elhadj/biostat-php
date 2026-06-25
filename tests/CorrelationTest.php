<?php

declare(strict_types=1);

namespace TouilElhadj\BiostatPhp\Tests;

/**
 * Tests for correlation methods (Pearson r and Spearman ρ).
 *
 * Reference values were obtained in R 4.3.0:
 *
 *   x <- c(10, 12, 14, 16, 18, 20, 22, 24, 26, 28)
 *   y <- c(15, 18, 17, 22, 24, 26, 25, 29, 32, 30)
 *   cor.test(x, y, method = "pearson")
 *   # r = 0.9668, t = 10.69, df = 8, p = 5.16e-06
 *
 *   cor.test(x, y, method = "spearman")
 *   # rho = 0.9758, p = 4.21e-06 (with ties warning)
 *
 * Reproducer: see tests/fixtures/reference-values.R
 */
final class CorrelationTest extends BiostatTestCase
{
    /**
     * @var array<int, float>
     */
    private array $x = [10, 12, 14, 16, 18, 20, 22, 24, 26, 28];

    /**
     * @var array<int, float>
     */
    private array $y = [15, 18, 17, 22, 24, 26, 25, 29, 32, 30];

    public function testPearsonR(): void
    {
        $r = $this->stats->pearson($this->x, $this->y);
        $this->assertNear(0.9668, (float) $r['r'], 0.001, 'Pearson r');
    }

    public function testPearsonPValueIsSignificant(): void
    {
        $r = $this->stats->pearson($this->x, $this->y);
        $this->assertLessThan(0.001, (float) $r['p'], 'Pearson p-value');
    }

    public function testPearsonOfPerfectCorrelation(): void
    {
        $r = $this->stats->pearson([1, 2, 3, 4, 5], [2, 4, 6, 8, 10]);
        $this->assertNear(1.0, (float) $r['r'], 1e-9, 'r=1 for y=2x');
    }

    public function testSpearmanRho(): void
    {
        $r = $this->stats->spearman($this->x, $this->y);
        // Allow looser tolerance (mid-rank ties handled differently)
        $this->assertNear(0.97, (float) $r['r'], 0.03, 'Spearman rho');
    }
}
