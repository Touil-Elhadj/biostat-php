<?php

/**
 * Minimal PHPUnit shim used only by run-tests-manually.php in
 * environments where Composer cannot install the real PHPUnit.
 *
 * In a normal environment this file is NEVER loaded — the autoloader
 * picks up phpunit/phpunit from vendor/ instead.
 */

declare(strict_types=1);

namespace PHPUnit\Framework {

    use Throwable;

    class AssertionFailedError extends \RuntimeException {}

    abstract class TestCase
    {
        protected function setUp(): void {}
        protected function tearDown(): void {}

        public function __construct() {
            $this->setUp();
        }

        protected function assertTrue($cond, string $msg = ''): void
        {
            if ($cond !== true) {
                throw new AssertionFailedError($msg ?: 'expected true, got ' . var_export($cond, true));
            }
        }

        protected function assertFalse($cond, string $msg = ''): void
        {
            if ($cond !== false) {
                throw new AssertionFailedError($msg ?: 'expected false');
            }
        }

        protected function assertSame($expected, $actual, string $msg = ''): void
        {
            if ($expected !== $actual) {
                throw new AssertionFailedError($msg ?: "expected " . var_export($expected, true) . ", got " . var_export($actual, true));
            }
        }

        protected function assertEquals($expected, $actual, string $msg = '', float $delta = 0.0): void
        {
            if (is_numeric($expected) && is_numeric($actual) && $delta > 0) {
                if (abs((float)$expected - (float)$actual) > $delta) {
                    throw new AssertionFailedError($msg ?: "expected $expected ±$delta, got $actual");
                }
            } elseif ($expected != $actual) {
                throw new AssertionFailedError($msg ?: "expected " . var_export($expected, true) . ", got " . var_export($actual, true));
            }
        }

        protected function assertLessThan($limit, $actual, string $msg = ''): void
        {
            if (!($actual < $limit)) {
                throw new AssertionFailedError($msg ?: "expected $actual < $limit");
            }
        }

        protected function assertLessThanOrEqual($limit, $actual, string $msg = ''): void
        {
            if (!($actual <= $limit)) {
                throw new AssertionFailedError($msg ?: "expected $actual <= $limit");
            }
        }

        protected function assertGreaterThan($limit, $actual, string $msg = ''): void
        {
            if (!($actual > $limit)) {
                throw new AssertionFailedError($msg ?: "expected $actual > $limit");
            }
        }

        protected function assertGreaterThanOrEqual($limit, $actual, string $msg = ''): void
        {
            if (!($actual >= $limit)) {
                throw new AssertionFailedError($msg ?: "expected $actual >= $limit");
            }
        }

        protected function assertIsArray($v, string $msg = ''): void
        {
            if (!is_array($v)) {
                throw new AssertionFailedError($msg ?: 'expected array');
            }
        }

        protected function assertArrayHasKey($k, $arr, string $msg = ''): void
        {
            if (!is_array($arr) || !array_key_exists($k, $arr)) {
                throw new AssertionFailedError($msg ?: "missing key $k");
            }
        }

        protected function assertArrayNotHasKey($k, $arr, string $msg = ''): void
        {
            if (is_array($arr) && array_key_exists($k, $arr)) {
                throw new AssertionFailedError($msg ?: "unexpected key $k");
            }
        }

        protected function assertCount(int $n, $arr, string $msg = ''): void
        {
            if (count($arr) !== $n) {
                throw new AssertionFailedError($msg ?: "expected count $n, got " . count($arr));
            }
        }

        protected function assertNotEmpty($v, string $msg = ''): void
        {
            if (empty($v)) {
                throw new AssertionFailedError($msg ?: 'expected non-empty');
            }
        }

        protected function assertNotNull($v, string $msg = ''): void
        {
            if ($v === null) {
                throw new AssertionFailedError($msg ?: 'expected not null');
            }
        }
    }
}
