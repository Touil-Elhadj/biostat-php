<?php

declare(strict_types=1);

namespace TouilElhadj\BiostatPhp\Tests;

/**
 * Tests for comparison-of-means methods.
 *
 * Reference values were obtained in R 4.3.0:
 *
 *   a <- c(10, 12, 11, 13, 14, 12, 11)
 *   b <- c(15, 17, 16, 18, 16, 17, 15)
 *   t.test(a, b)              # t = -6.7116, df = 11.6, p = 2.71e-05
 *
 *   g1 <- c(10, 12, 11)
 *   g2 <- c(15, 17, 16)
 *   g3 <- c(20, 22, 21)
 *   summary(aov(values ~ groups,
 *               data = stack(list(g1 = g1, g2 = g2, g3 = g3))))
 *   # F(2, 6) = 75, p = 6.34e-05
 */
final class MeansTest extends BiostatTestCase
{
    public function testWelchTStatistic(): void
    {
        $r = $this->stats->tTest(
            [10, 12, 11, 13, 14, 12, 11],
            [15, 17, 16, 18, 16, 17, 15]
        );
        $this->assertNear(-6.7116, (float) $r['t'], 0.01, 't-statistic');
    }

    public function testWelchSatterthwaiteDf(): void
    {
        $r = $this->stats->tTest(
            [10, 12, 11, 13, 14, 12, 11],
            [15, 17, 16, 18, 16, 17, 15]
        );
        $this->assertNear(11.6, (float) $r['df'], 0.1, 'Welch df');
    }

    public function testWelchPValueIsTiny(): void
    {
        $r = $this->stats->tTest(
            [10, 12, 11, 13, 14, 12, 11],
            [15, 17, 16, 18, 16, 17, 15]
        );
        // p ≈ 2.71e-5; the CDF approximation prints it as 0.0001 after rounding —
        // we just check it is < 0.001.
        $this->assertLessThan(0.001, (float) $r['p']);
    }

    public function testAnovaFStatistic(): void
    {
        $r = $this->stats->anova([
            [10, 12, 11],
            [15, 17, 16],
            [20, 22, 21],
        ]);
        $this->assertNear(75.0, (float) $r['F'], 0.01, 'F-statistic');
        $this->assertSame(2, $r['dfB']);
        $this->assertSame(6, $r['dfW']);
        $this->assertLessThan(0.001, (float) $r['p']);
    }
}
