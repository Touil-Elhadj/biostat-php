<?php

/**
 * ════════════════════════════════════════════════════════════════════
 *  run-tests-manually.php — Lightweight test harness
 * ────────────────────────────────────────────────────────────────────
 *  This script is a stop-gap for environments where Composer / PHPUnit
 *  cannot be installed (offline machines, restricted shared hosting).
 *
 *  In normal development you should use:
 *      composer install
 *      composer test
 *
 *  This file emulates a tiny subset of PHPUnit: it discovers Test*::test*
 *  methods and runs them with a shared `BiostatAnalysis` instance plus
 *  an assertNear() helper. It does NOT replace PHPUnit's coverage,
 *  groups, data providers or fixtures.
 * ════════════════════════════════════════════════════════════════════
 */

declare(strict_types=1);

require_once __DIR__ . '/autoload.php';
require_once __DIR__ . '/tests/_phpunit_shim.php';

use TouilElhadj\BiostatPhp\BiostatAnalysis;

class TinyTestRunner
{
    private int $passed = 0;
    private int $failed = 0;
    /** @var array<int, string> */
    private array $failures = [];

    public function run(): int
    {
        $testDir = __DIR__ . '/tests';
        // Load the abstract base class first
        require_once $testDir . '/BiostatTestCase.php';
        foreach (glob($testDir . '/*Test.php') as $file) {
            require_once $file;
        }

        echo "\nbiostat-php — manual test harness\n";
        echo str_repeat('━', 70) . "\n\n";

        $declared = get_declared_classes();
        $tests    = array_filter(
            $declared,
            fn($c) => str_starts_with($c, 'TouilElhadj\\BiostatPhp\\Tests\\')
                  && !str_ends_with($c, 'TestCase')
        );

        foreach ($tests as $cls) {
            $short = substr($cls, strrpos($cls, '\\') + 1);
            echo "── $short ";
            echo str_repeat('─', max(0, 60 - strlen($short))) . "\n";
            $this->runClass($cls);
        }

        echo "\n" . str_repeat('━', 70) . "\n";
        printf("Tests: %d passed, %d failed\n", $this->passed, $this->failed);
        if ($this->failed > 0) {
            echo "\nFailures:\n";
            foreach ($this->failures as $f) {
                echo "  • $f\n";
            }
            return 1;
        }
        echo "✓ All tests passed.\n";
        return 0;
    }

    private function runClass(string $cls): void
    {
        $rc = new ReflectionClass($cls);
        if ($rc->isAbstract()) return;

        $instance = $rc->newInstance();
        $stats    = new BiostatAnalysis();

        // Inject $stats into the protected property
        $statsProp = (new ReflectionClass($cls))->getParentClass()->getProperty('stats');
        $statsProp->setAccessible(true);
        $statsProp->setValue($instance, $stats);

        foreach ($rc->getMethods(ReflectionMethod::IS_PUBLIC) as $m) {
            if (!str_starts_with($m->name, 'test')) continue;

            try {
                $m->invoke($instance);
                $this->passed++;
                printf("  ✓ %s\n", $m->name);
            } catch (Throwable $e) {
                $this->failed++;
                $msg = sprintf('%s::%s — %s', $rc->getShortName(), $m->name, $e->getMessage());
                $this->failures[] = $msg;
                printf("  ✗ %s\n      %s\n", $m->name, $e->getMessage());
            }
        }
    }
}

// ── shim loaded at top of file ──────────────────────────────────────

exit((new TinyTestRunner())->run());
