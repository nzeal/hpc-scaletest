#!/usr/bin/env python3
"""
HPC-ScaleTest: Modular HPC Scaling Test Framework
Main CLI entry point.
"""

import argparse
import sys
import importlib.util
from pathlib import Path
from core.test_definition import Test
from engine.runner import TestRunner
from utils.logging_config import setup_logging
from core.factory import BackendFactory
from core.types import SchedulerBackend, LauncherBackend, ModuleBackend, BuildBackend

def load_test_file(test_file: str):
    spec = importlib.util.spec_from_file_location("test_module", test_file)
    module = importlib.util.module_from_spec(spec)
    if spec.loader is None:
        raise ImportError(f"Cannot load spec from {test_file}")
    spec.loader.exec_module(module)
    # Assume it defines create_test()
    if hasattr(module, 'create_test'):
        return module.create_test()
    raise ValueError("Test file must define create_test()")

def list_backends():
    print("Available Schedulers:", [e.name for e in SchedulerBackend])
    print("Available Launchers:", [e.name for e in LauncherBackend])
    print("Available Module Systems:", [e.name for e in ModuleBackend])
    print("Available Build Systems:", [e.name for e in BuildBackend])

def validate_test(test: Test):
    config = test.get_config()
    if getattr(config, 'scaling', None) is None:
        raise ValueError("Scaling config required")
    print("Test validation passed.")

def main():
    parser = argparse.ArgumentParser(description="HPC-ScaleTest CLI")
    subparsers = parser.add_subparsers(dest='command')

    run_parser = subparsers.add_parser('run')
    run_parser.add_argument('--test', required=True, help='Test definition file')
    run_parser.add_argument('--scaling', default='strong', choices=['strong', 'weak'])
    run_parser.add_argument('--max-nodes', type=int, default=4)
    run_parser.add_argument('--backend', default='local')
    run_parser.add_argument('--output')
    run_parser.add_argument('-v', '--verbose', action='store_true')

    validate_parser = subparsers.add_parser('validate')
    validate_parser.add_argument('--test', required=True)

    list_parser = subparsers.add_parser('list-backends')

    args = parser.parse_args()

    # Setup logging if verbose flag exists and is True
    verbose = getattr(args, 'verbose', False)
    setup_logging(verbose)

    if args.command == 'run':
        test = load_test_file(args.test)
        # Override scaling if provided
        if hasattr(test, "set_scaling"):
            test.set_scaling(args.scaling, args.max_nodes)
        if hasattr(test, "set_backend"):
            test.set_backend(args.backend)
        runner = TestRunner(test)
        output_path = Path(args.output) if args.output else None
        runner.run(output_path)

    elif args.command == 'validate':
        test = load_test_file(args.test)
        validate_test(test)

    elif args.command == 'list-backends':
        list_backends()

    else:
        parser.print_help()

if __name__ == '__main__':
    main()
