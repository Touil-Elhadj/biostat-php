<?php

declare(strict_types=1);

namespace TouilElhadj\BiostatPhp\Tests;

use PHPUnit\Framework\TestCase;
use TouilElhadj\BiostatPhp\BiostatAnalysis;

/**
 * Base class for every BiostatPhp test case.
 *
 * Provides:
 *  • {@see assertNear()} — assertion with explicit numerical tolerance,
 *    used throughout the suite to compare PHP results to the closed-form
 *    reference values pre-computed in R 4.x or IBM SPSS 25.
 *  • a shared `$stats` instance so individual tests can stay focused on
 *    the assertions.
 */
abstract class BiostatTestCase extends TestCase
{
    protected BiostatAnalysis $stats;

    protected function setUp(): void
    {
        $this->stats = new BiostatAnalysis();
    }

    /**
     * Assert that |actual − expected| ≤ tolerance, with a descriptive
     * failure message that includes both values and the gap.
     *
     * @param float       $expected  the reference value (typically from R)
     * @param float       $actual    the value produced by BiostatPhp
     * @param float       $tolerance maximum acceptable absolute difference
     * @param string|null $label     optional label printed on failure
     */
    protected function assertNear(
        float $expected,
        float $actual,
        float $tolerance,
        ?string $label = null
    ): void {
        $diff = abs($expected - $actual);
        $msg  = sprintf(
            '%sexpected %.6g (±%.6g), got %.6g — |Δ| = %.6g',
            $label ? "[$label] " : '',
            $expected,
            $tolerance,
            $actual,
            $diff
        );
        $this->assertLessThanOrEqual($tolerance, $diff, $msg);
    }
}
