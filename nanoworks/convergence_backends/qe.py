"""Quantum ESPRESSO static-energy backend for convergence workflows."""

import math
from pathlib import Path

from nanoworks.convergence import StaticEnergyResult
from nanoworks.engine import load_engine_module


class QEStaticEnergyBackend:
    """Adapt the existing QE SCF runner to the convergence interface."""

    name = 'QE'

    def __init__(
        self,
        pseudopotentials,
        pseudo_dir,
        executable='pw.x',
        engine_module=None,
    ):
        self.pseudopotentials = dict(pseudopotentials)
        self.pseudo_dir = Path(pseudo_dir)
        self.executable = str(executable)
        self._engine_module = engine_module

    @property
    def engine_module(self):
        """Load the QE implementation only when a calculation is requested."""
        if self._engine_module is None:
            self._engine_module = load_engine_module('QE')

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
        """Run one QE SCF calculation and return a common energy result."""
        workdir = Path(workdir)
        workdir.mkdir(parents=True, exist_ok=True)
        kpoint_settings = dict(kpoint_settings)
        settings = dict(settings)
        kpoint_density = kpoint_settings.get('density')
        kpoint_size = kpoint_settings.get('size')
        if kpoint_density is None and kpoint_size is None:
            kpoint_size = (5, 5, 5)

        calculation = self.engine_module.run_scf(
            atoms=atoms,
            input_file=workdir / 'qe-scf.in',
            output_file=workdir / 'qe-scf.out',
            state_dir=workdir / 'state',
            pseudopotentials=self.pseudopotentials,
            pseudo_dir=self.pseudo_dir,
            cutoff_ev=float(cutoff_ev),
            kpoint_density=kpoint_density,
            kpoint_size=kpoint_size,
            gamma=bool(kpoint_settings.get('gamma', False)),
            total_charge=settings.get('total_charge', 0.0),
            nbands=settings.get('nbands'),
            spinpol=bool(settings.get('spinpol', False)),
            magnetic_moments=settings.get('magnetic_moments'),
            setup_params=settings.get('setup_params'),
            xc_calc=settings.get('xc_calc', 'PBE'),
            pseudo_xc=settings.get('pseudo_xc', 'pbe'),
            exx_fraction=settings.get('exx_fraction'),
            omega=settings.get('omega'),
            occupation=settings.get('occupation'),
            parallel_cores=parallel_cores,
            executable=self.executable,
            prefix=settings.get('prefix', 'nanoworks-converge'),
            exx_additional_kpoints=settings.get(
                'exx_additional_kpoints'
            ),
        )

        qe_result = calculation.get('result', {})
        total_energy_ev = qe_result.get('total_energy_ev')

        if total_energy_ev is None:
            raise RuntimeError(
                'QE SCF output did not contain a total energy.'
            )

        total_energy_ev = float(total_energy_ev)
        if not math.isfinite(total_energy_ev):
            raise RuntimeError(
                'QE SCF output contained a non-finite total energy.'
            )

        return StaticEnergyResult(
            engine=self.name,
            total_energy_ev=total_energy_ev,
            metadata={
                'input_file': str(calculation['input_file']),
                'output_file': str(calculation['output_file']),
                'state_dir': str(calculation['state_dir']),
                'kpoint_size': tuple(calculation['kpoint_size']),
                'qe_version': qe_result.get('qe_version'),
            },
        )
