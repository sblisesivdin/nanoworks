"""GPAW static-energy backend for convergence workflows."""

import math
from pathlib import Path

from nanoworks.convergence import StaticEnergyResult
from nanoworks.engine import load_engine_module


class GPAWStaticEnergyBackend:
    """Adapt Nanoworks GPAW PW calculators to the convergence interface."""

    name = 'GPAW'

    def __init__(self, engine_module=None):
        self._engine_module = engine_module

    @property
    def engine_module(self):
        """Import GPAW only when a GPAW calculation is requested."""
        if self._engine_module is None:
            self._engine_module = load_engine_module('GPAW')

        return self._engine_module

    def calculate_static_energy(
        self,
        atoms,
        *,
        cutoff_ev,
        kpoint_settings,
        workdir,
        settings,
        parallel_cores,
    ):
        """Run one GPAW PW calculation and return a common energy result."""
        workdir = Path(workdir)
        workdir.mkdir(parents=True, exist_ok=True)
        kpoint_settings = dict(kpoint_settings)
        settings = dict(settings)
        engine = self.engine_module

        calculation_atoms = atoms.copy()
        spinpol = bool(settings.get('spinpol', False))
        magnetic_moments = settings.get('magnetic_moments')

        if spinpol and magnetic_moments is not None:
            calculation_atoms.set_initial_magnetic_moments(
                magnetic_moments
            )

        xc_calc = settings.get('xc_calc', 'PBE')
        actual_xc, setups, _ = engine.resolve_xc_and_setups(
            xc_calc,
            settings.get('setup_params'),
        )

        mixer = settings.get('mixer')
        if mixer is None:
            mixer = engine.create_default_mixer()

        common = {
            'cutoff': float(cutoff_ev),
            'mixer': mixer,
            'charge': settings.get('total_charge', 0.0),
            'spinpol': spinpol,
            'txt': str(workdir / 'gpaw-scf.txt'),
            'convergence': settings.get('convergence', {}),
            'occupations': settings.get('occupation'),
            'kpoint_density': kpoint_settings.get('density'),
            'kpoint_size': kpoint_settings.get('size', (5, 5, 5)),
            'gamma': bool(kpoint_settings.get('gamma', False)),
            'nbands': settings.get('nbands'),
        }

        if engine.is_hybrid(actual_xc):
            calculator = engine.create_hybrid_pw_ground_calc(
                xc_calc=xc_calc,
                exx_fraction=settings.get('exx_fraction'),
                omega=settings.get('omega'),
                backend=settings.get('xc_backend', 'pw'),
                **common,
            )
        else:
            calculator = engine.create_regular_pw_ground_calc(
                xc=actual_xc,
                setups=setups,
                parallel={'domain': parallel_cores},
                **common,
            )

        calculation_atoms.calc = calculator
        total_energy_ev = float(
            calculation_atoms.get_potential_energy()
        )

        if not math.isfinite(total_energy_ev):
            raise RuntimeError(
                'GPAW returned a non-finite total energy.'
            )

        metadata = {
            'log_file': str(workdir / 'gpaw-scf.txt'),
            'kpoint_density': kpoint_settings.get('density'),
            'kpoint_size': kpoint_settings.get('size', (5, 5, 5)),
            'hybrid': bool(engine.is_hybrid(actual_xc)),
        }

        return StaticEnergyResult(
            engine=self.name,
            total_energy_ev=total_energy_ev,
            metadata=metadata,
        )
