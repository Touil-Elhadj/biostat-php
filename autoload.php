<?php

/**
 * Minimal PSR-4 autoloader for development before `composer install`.
 * In production this file is replaced by Composer's autoloader.
 */

spl_autoload_register(function (string $class): void {
    $prefix = 'TouilElhadj\\BiostatPhp\\';
    if (strpos($class, $prefix) !== 0) {
        return;
    }
    $relative = substr($class, strlen($prefix));
    $file = __DIR__ . '/src/' . str_replace('\\', '/', $relative) . '.php';
    if (is_file($file)) {
        require $file;
    }
});
