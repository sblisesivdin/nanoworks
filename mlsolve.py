#!/usr/bin/env python3
"""
mlsolve.py - ML force-field solver

The script accepts a geometry file path and parses the args, and prints a message. Full features
will be added in the future.
"""

import sys
import argparse
import ast
from pathlib import Path

try:
    from ase.io import read, write
    from ase.optimize import BFGS
except ImportError:
    sys.exit("Error: ASE is required. Install with: pip install ase")


# ---------------------------------------------------------
# Calculator Factory
# ---------------------------------------------------------

def get_ml_calculator(model_type, device="cpu", **kwargs):
    model_type = model_type.lower()

    print(f"\n[Calculator Factory] Model: {model_type}, device: {device}")

    if model_type == "mace":
        try:
            from mace.calculators import mace_mp
        except ImportError:
            sys.exit("Error: MACE not installed. pip install mace-torch")
        return mace_mp(model="medium", device=device)

    elif model_type == "chgnet":
        try:
            from chgnet.model import CHGNet
            from chgnet.model.model import CHGNetCalculator
        except ImportError:
            sys.exit("Error: CHGNet not installed. pip install chgnet")
        model = CHGNet.load()
        return CHGNetCalculator(model=model, use_device=device)

    elif model_type == "sevennet":
        try:
            from sevenn.sevennet_calculator import SevenNetCalculator
        except ImportError:
            sys.exit("Error: SevenNet not installed. pip install sevenn")
        return SevenNetCalculator(model="7net-0", device=device)

    else:
        sys.exit(f"Error: unknown model '{model_type}'")


# ---------------------------------------------------------
# Task Functions
# ---------------------------------------------------------

def run_static(atoms, config):
    print("\n--- Static calculation ---")
    atoms.calc = get_ml_calculator(config["model"], config["device"], **config)

    try:
        energy = atoms.get_potential_energy()
        forces = atoms.get_forces()
        print(f"Potential energy: {energy:.6f} eV")
        print(f"Max force: {forces.max():.6f} eV/Å")
    except Exception as e:
        print(f"Static calculation failed: {e}")


def run_optimize(atoms, config):
    print("\n--- Geometry optimization ---")
    atoms.calc = get_ml_calculator(config["model"], config["device"], **config)

    fmax = config.get("fmax", 0.05)
    steps = config.get("steps", 200)

    dyn = BFGS(atoms)
    try:
        dyn.run(fmax=fmax, steps=steps)
        write("optimized.cif", atoms)
        print("Optimization finished.")
        print(f"Final energy: {atoms.get_potential_energy():.6f} eV")
        print("Output written to optimized.cif")
    except Exception as e:
        print(f"Optimization failed: {e}")


# ---------------------------------------------------------
# Argument Parsing
# ---------------------------------------------------------

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="MLSolve: unified ML force-field solver."
    )
    parser.add_argument("-g", "--geometry", required=True, type=str)
    parser.add_argument("-i", "--input", required=False, type=str)
    return parser.parse_args()


# ---------------------------------------------------------
# Main
# ---------------------------------------------------------

def main():
    args = parse_arguments()
    geom_path = Path(args.geometry)

    if not geom_path.exists():
        sys.exit(f"Error: geometry file '{geom_path}' not found.")

    try:
        atoms = read(geom_path)
    except Exception as e:
        sys.exit(f"Error reading geometry: {e}")

    user_config = {}
    if args.input:
        try:
            user_config = ast.literal_eval(args.input)
            if not isinstance(user_config, dict):
                raise ValueError
        except Exception:
            sys.exit("Error: -i must be a valid dictionary")

    config = {
        "model": "mace",
        "task": "static",
        "device": "cpu",
        "fmax": 0.05,
        "steps": 200,
    }
    config.update(user_config)

    print("MLSolve")
    print(f"Structure: {atoms.get_chemical_formula()}")
    print(f"Atoms: {len(atoms)}")

    task = config["task"].lower()
    if task == "static":
        run_static(atoms, config)
    elif task == "optimize":
        run_optimize(atoms, config)
    else:
        sys.exit(f"Error: unknown task '{task}'")


if __name__ == "__main__":
    main()

