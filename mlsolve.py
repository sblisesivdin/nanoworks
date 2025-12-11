#!/usr/bin/env python3
"""
mlsolve.py - ML force-field solver

This is an intentionally tiny starter file for the project history.
It accepts a geometry file path and prints a message. Full features
will be added in the future.
"""

import sys
from pathlib import Path

def main():
    if len(sys.argv) < 2:
        print("Usage: python mlsolve.py <geometry-file>")
        sys.exit(1)

    geom = Path(sys.argv[1])
    if not geom.exists():
        print(f"Error: geometry file '{geom}' not found.")
        sys.exit(2)

    print("mlsolve")
    print(f"Geometry file: {geom}")
    print("No calculators available in this time.")

if __name__ == "__main__":
    main()
