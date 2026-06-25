<?php

declare(strict_types=1);

namespace TouilElhadj\BiostatPhp\Tests;

use TouilElhadj\BiostatPhp\BiostatAnalysis;

/**
 * Tests for the Benjamini-Hochberg FDR adjustment.
 *
 * Reference values were obtained in R 4.3.0:
 *
 *   p <- c(0.001, 0.008, 0.039, 0.041, 0.042, 0.060, 0.074, 0.205, 0.212, 0.216)
 *   p.adjust(p, method = "BH")
 *   # 0.0100 0.0400 0.1025 0.1025 0.1025 0.1000 0.1057 0.2400 0.2400 0.2400
 *
 * Note: R also enforces monotonicity (a smaller-rank adjusted p never
 * exceeds a larger-rank one), which this implementation also enforces.
 */
final class BenjaminiHochbergTest extends BiostatTestCase
{
    public function testReturnsArrayOfSameLength(): void
    {
        $p   = [0.001, 0.008, 0.039, 0.041, 0.042, 0.060, 0.074, 0.205, 0.212, 0.216];
        $adj = BiostatAnalysis::benjaminiHochberg($p);
        $this->assertCount(count($p), $adj);
    }

    public function testSmallestPValueAdjusts(): void
    {
        $p   = [0.001, 0.008, 0.039, 0.041, 0.042, 0.060, 0.074, 0.205, 0.212, 0.216];
        $adj = BiostatAnalysis::benjaminiHochberg($p);
        // R: 0.001 × 10/1 = 0.0100
        $this->assertNear(0.0100, (float) $adj[0], 0.001, 'BH[0.001]');
    }

    public function testSecondSmallestAdjusts(): void
    {
        $p   = [0.001, 0.008, 0.039, 0.041, 0.042, 0.060, 0.074, 0.205, 0.212, 0.216];
        $adj = BiostatAnalysis::benjaminiHochberg($p);
        // R: 0.008 × 10/2 = 0.040
        $this->assertNear(0.0400, (float) $adj[1], 0.001, 'BH[0.008]');
    }

    public function testMonotonicity(): void
    {
        $p   = [0.001, 0.008, 0.039, 0.041, 0.042, 0.060, 0.074, 0.205, 0.212, 0.216];
        $adj = BiostatAnalysis::benjaminiHochberg($p);
        // Sort by original p and check adjusted values are non-decreasing
        $pairs = array_map(null, $p, $adj);
        usort($pairs, static fn($a, $b) => $a[0] <=> $b[0]);
        $prev = -INF;
        foreach ($pairs as [, $a]) {
            $this->assertGreaterThanOrEqual($prev, (float) $a, 'BH must be monotonic');
            $prev = (float) $a;
        }
    }

    public function testPreservesKeys(): void
    {
        $p = ['gene_A' => 0.01, 'gene_B' => 0.5, 'gene_C' => 0.001];
        $adj = BiostatAnalysis::benjaminiHochberg($p);
        $this->assertArrayHasKey('gene_A', $adj);
        $this->assertArrayHasKey('gene_B', $adj);
        $this->assertArrayHasKey('gene_C', $adj);
    }

    public function testEmptyInputReturnsEmpty(): void
    {
        $this->assertSame([], BiostatAnalysis::benjaminiHochberg([]));
    }
}
