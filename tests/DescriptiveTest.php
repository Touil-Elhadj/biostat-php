<?php

declare(strict_types=1);

namespace TouilElhadj\BiostatPhp\Tests;

/**
 * Tests for descriptive statistics (mean, std, median, quantile).
 *
 * All reference values were obtained in R 4.3.0:
 *   x <- c(2.3, 4.1, 5.7, 6.2, 7.8, 8.4, 9.1, 10.5, 11.2, 12.9)
 *   mean(x)                # 7.82
 *   sd(x)                  # 3.296732
 *   median(x)              # 8.10
 *   quantile(x, 0.25)      # 5.825   (R type = 7, default)
 *   quantile(x, 0.75)      # 10.150
 *
 * Reproducer: see tests/fixtures/reference-values.R
 */
final class DescriptiveTest extends BiostatTestCase
{
    /**
     * @var array<int, float>
     */
    private array $x = [2.3, 4.1, 5.7, 6.2, 7.8, 8.4, 9.1, 10.5, 11.2, 12.9];

    public function testMean(): void
    {
        $this->assertNear(7.82, (float) $this->stats->mean($this->x), 1e-4, 'mean');
    }

    public function testStd(): void
    {
        // R: sd(x) — verified value is 3.296732
        $this->assertNear(3.296732, (float) $this->stats->std($this->x), 1e-4, 'std (n-1)');
    }

    public function testMedianOdd(): void
    {
        $this->assertSame(3, $this->stats->median([1, 2, 3, 4, 5]));
    }

    public function testMedianEven(): void
    {
        $this->assertSame(8.10, round((float) $this->stats->median($this->x), 2));
    }

    public function testQuantileQ25(): void
    {
        $this->assertNear(5.825, (float) $this->stats->quantile($this->x, 0.25), 1e-3, 'Q25');
    }

    public function testQuantileQ75(): void
    {
        // R type=7: 10.15  (NOT 10.325 as previously written by mistake)
        $this->assertNear(10.150, (float) $this->stats->quantile($this->x, 0.75), 1e-3, 'Q75');
    }

    public function testMeanOfEmptyArrayIsZero(): void
    {
        // Defensive: empty input must not divide-by-zero
        $this->assertSame(0, $this->stats->mean([]));
    }
}
