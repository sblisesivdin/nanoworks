import unittest
import tempfile
from pathlib import Path
from unittest.mock import patch

from nanoworks.engine import load_engine_module
from ase.units import Bohr
from ase import Atoms
from ase.build import bulk
from nanoworks.engine.qe import (
    QE_REFERENCE_VERSION,
    THZ_PER_CM_MINUS_ONE,
    ev_to_rydberg,
    build_control_settings,
    build_system_settings,
    build_cell_parameters,
    build_atomic_positions,
    build_atomic_species,
    build_kpoint_settings,
    build_band_path,
    render_band_kpoints,
    build_qe_exx_additional_kpoints,
    render_qe_exx_additional_kpoints,
    build_occupation_settings,
    build_electrons_settings,
    format_qe_value,
    render_namelist,
    render_pw_input,
    render_scf_input,
    render_nscf_input,
    render_bands_input,
    render_bands_postprocess_input,
    render_relax_input,
    rydberg_to_ev,
    build_qe_launcher,
    resolve_qe_executable,
    build_qe_command,
    run_qe_program,
    parse_qe_auxiliary_output,
    parse_epsilon_data_file,
    parse_matdyn_frequency_file,
    parse_matdyn_dos_file,
    write_matdyn_band_data,
    write_matdyn_dos_data,
    calculate_phonon_thermal_properties,
    write_phonon_thermal_properties,
    parse_pw_output,
    resolve_qe_band_reference,
    parse_pw_bands_output,
    parse_bands_x_output,
    resolve_qe_kpoint_size,
    resolve_qe_occupation,
    validate_qe_version,
    resolve_qe_xc_settings,
    validate_qe_xc,
    run_scf,
    run_nscf,
    run_bands,
    run_bands_postprocess,
    run_hybrid_bands,
    run_hybrid_dos,
    has_qe_state,
    render_dos_input,
    render_epsilon_input,
    run_epsilon,
    run_dos,
    parse_dos_output,
    render_projwfc_input,
    run_projwfc,
    parse_projwfc_pdos_file,
    aggregate_projwfc_pdos,
    prepare_qe_band_data,
    resolve_qe_cell_dofree,
    resolve_qe_relaxation_settings,
    parse_pw_relaxed_structure,
    run_relax,
    build_qe_magnetic_species,
    run_band_projections,
    parse_projwfc_band_file,
    prepare_qe_band_projection_data,
    render_pp_input,
    run_pp_density,
    render_qe_hubbard_card,
    resolve_qe_hubbard,
    resolve_qe_phonon_qpoint_grid,
    build_ph_settings,
    render_ph_input,
    build_q2r_settings,
    render_q2r_input,
    build_matdyn_band_settings,
    render_matdyn_qpoints,
    render_matdyn_band_input,
    build_matdyn_dos_settings,
    render_matdyn_dos_input,
    run_ph,
    run_q2r,
    run_matdyn_band,
    run_matdyn_dos,
)


class TestQEEngine(unittest.TestCase):

    def test_reference_version_is_qe_72(self):
        self.assertEqual(QE_REFERENCE_VERSION, (7, 2))

    def test_build_and_render_qe_exx_additional_kpoints(self):
        settings = build_qe_exx_additional_kpoints(
            band_path={
                'kpoints': [
                    (0.1, 0.0, 0.0),
                    (0.2, 0.0, 0.0),
                ],
                'npoints': 2,
            },
            qpoint_grid=(2, 1, 1),
        )

        self.assertEqual(
            settings['band_indices'],
            [0, 1],
        )
        self.assertEqual(
            settings['helper_indices'],
            [2, 3],
        )
        self.assertEqual(
            settings['kpoints'][2],
            (0.6, 0.0, 0.0),
        )
        self.assertEqual(
            settings['kpoints'][3],
            (0.7, 0.0, 0.0),
        )

        text = render_qe_exx_additional_kpoints(
            settings
        )
        lines = text.splitlines()

        self.assertEqual(
            lines[0],
            'ADDITIONAL_K_POINTS crystal',
        )
        self.assertEqual(
            lines[1],
            '4',
        )
        self.assertTrue(
            all(
                line.endswith(' 0.0')
                for line in lines[2:]
            )
        )

    def test_resolve_qe_xc_settings_for_pbe(self):
        settings = resolve_qe_xc_settings(
            'PBE'
        )

        self.assertEqual(
            settings,
            {
                'name': 'PBE',
                'input_dft': 'PBE',
                'hybrid': False,
                'exx_fraction': None,
                'screening_parameter': None,
                'pseudo_xc': 'pbe',
            },
        )

    def test_resolve_qe_xc_settings_for_hse06(self):
        settings = resolve_qe_xc_settings(
            'HSE06'
        )

        self.assertEqual(
            settings['input_dft'],
            'HSE',
        )
        self.assertTrue(
            settings['hybrid']
        )
        self.assertEqual(
            settings['exx_fraction'],
            0.25,
        )
        self.assertEqual(
            settings['screening_parameter'],
            0.106,
        )

    def test_resolve_qe_xc_settings_for_hse03_overrides(self):
        settings = resolve_qe_xc_settings(
            'HSE-03',
            exx_fraction=0.30,
            omega=0.16,
        )

        self.assertEqual(
            settings['name'],
            'HSE03',
        )
        self.assertEqual(
            settings['exx_fraction'],
            0.30,
        )
        self.assertEqual(
            settings['screening_parameter'],
            0.16,
        )

    def test_resolve_qe_xc_settings_for_pbe0(self):
        settings = resolve_qe_xc_settings(
            'PBE-0'
        )

        self.assertEqual(
            settings['input_dft'],
            'PBE0',
        )
        self.assertEqual(
            settings['exx_fraction'],
            0.25,
        )
        self.assertIsNone(
            settings['screening_parameter']
        )

    def test_validate_qe_xc_requires_explicit_hybrid_enablement(self):
        with self.assertRaisesRegex(
            ValueError,
            'is not enabled for this calculation stage yet',
        ):
            validate_qe_xc(
                'HSE06'
            )

        self.assertEqual(
            validate_qe_xc(
                'HSE06',
                allow_hybrid=True,
            ),
            'hse06',
        )

    def test_build_system_settings_adds_hybrid_controls(self):
        settings = build_system_settings(
            cutoff_ev=500.0,
            nat=2,
            ntyp=1,
            xc_calc='HSE03',
            exx_fraction=0.30,
            omega=0.16,
        )

        self.assertEqual(
            settings['input_dft'],
            'HSE',
        )
        self.assertEqual(
            settings['exx_fraction'],
            0.30,
        )
        self.assertEqual(
            settings['screening_parameter'],
            0.16,
        )

    def test_render_scf_input_adds_pbe0_controls(self):
        atoms = bulk(
            'Si',
            'diamond',
            a=5.43,
        )

        text = render_scf_input(
            atoms=atoms,
            pseudopotentials={
                'Si': 'Si.upf',
            },
            cutoff_ev=500.0,
            kpoint_size=(2, 2, 2),
            xc_calc='PBE0',
        )

        self.assertIn(
            "  input_dft = 'PBE0',",
            text,
        )
        self.assertIn(
            '  exx_fraction = 0.25,',
            text,
        )
        self.assertNotIn(
            'screening_parameter',
            text,
        )

    def test_render_nscf_input_rejects_hybrid_functional(self):
        atoms = bulk(
            'Si',
            'diamond',
            a=5.43,
        )

        with self.assertRaisesRegex(
            NotImplementedError,
            'does not support separate NSCF or bands',
        ):
            render_nscf_input(
                atoms=atoms,
                pseudopotentials={
                    'Si': 'Si.upf',
                },
                cutoff_ev=500.0,
                kpoint_size=(4, 4, 4),
                xc_calc='HSE06',
            )

    def test_render_bands_input_rejects_hybrid_functional(self):
        atoms = bulk(
            'Si',
            'diamond',
            a=5.43,
        )
        band_path = build_band_path(
            atoms=atoms,
            path='GX',
            npoints=5,
        )

        with self.assertRaisesRegex(
            NotImplementedError,
            'does not support separate NSCF or bands',
        ):
            render_bands_input(
                atoms=atoms,
                pseudopotentials={
                    'Si': 'Si.upf',
                },
                cutoff_ev=500.0,
                band_path=band_path,
                xc_calc='HSE03',
            )

    def test_resolve_qe_xc_settings_rejects_invalid_controls(self):
        invalid_requests = (
            (
                {
                    'xc_calc': 'PBE',
                    'exx_fraction': 0.25,
                },
                'only be used with',
            ),
            (
                {
                    'xc_calc': 'HSE06',
                    'exx_fraction': 0.0,
                },
                '0 < XC_exx_fraction',
            ),
            (
                {
                    'xc_calc': 'PBE0',
                    'omega': 0.11,
                },
                'only valid for screened HSE',
            ),
            (
                {
                    'xc_calc': 'B3LYP',
                },
                'supports PBE, HSE06, HSE03, and PBE0',
            ),
        )

        for request, message in invalid_requests:
            with self.subTest(request=request):
                with self.assertRaisesRegex(
                    ValueError,
                    message,
                ):
                    resolve_qe_xc_settings(
                        **request
                    )

    def test_resolve_qe_phonon_grid_from_diagonal_matrix(self):
        grid = resolve_qe_phonon_qpoint_grid(
            [
                [2, 0, 0],
                [0, 3, 0],
                [0, 0, 4],
            ]
        )

        self.assertEqual(
            grid,
            (2, 3, 4),
        )

    def test_resolve_qe_phonon_grid_from_diagonal_values(self):
        grid = resolve_qe_phonon_qpoint_grid(
            (2, 2, 1)
        )

        self.assertEqual(
            grid,
            (2, 2, 1),
        )

    def test_resolve_qe_phonon_grid_rejects_nondiagonal_matrix(self):
        with self.assertRaisesRegex(
            ValueError,
            'non-diagonal matrices',
        ):
            resolve_qe_phonon_qpoint_grid(
                [
                    [2, 1, 0],
                    [0, 2, 0],
                    [0, 0, 2],
                ]
            )

    def test_resolve_qe_phonon_grid_rejects_invalid_values(self):
        invalid_supercells = (
            (2, 2),
            (2, 0, 2),
            (2, 2.5, 2),
            (True, 2, 2),
            [
                [2, 0],
                [0, 2],
                [0, 0],
            ],
            [
                [2, 0, 0],
                2,
                [0, 0, 2],
            ],
        )

        for supercell in invalid_supercells:
            with self.subTest(supercell=supercell):
                with self.assertRaises(
                    (TypeError, ValueError)
                ):
                    resolve_qe_phonon_qpoint_grid(
                        supercell
                    )

    def test_ev_to_rydberg(self):
        self.assertAlmostEqual(
            ev_to_rydberg(13.605693122994),
            1.0,
        )

    def test_control_settings_defaults(self):
        settings = build_control_settings()

        self.assertEqual(settings['calculation'], 'scf')
        self.assertEqual(settings['prefix'], 'nanoworks')
        self.assertNotIn('pseudo_dir', settings)
        self.assertNotIn('outdir', settings)

    def test_control_settings_optional_paths(self):
        settings = build_control_settings(
            calculation='SCF',
            prefix='GaAs',
            pseudo_dir='/tmp/pseudos',
            outdir='/tmp/qe',
        )

        self.assertEqual(settings['calculation'], 'scf')
        self.assertEqual(settings['prefix'], 'GaAs')
        self.assertEqual(settings['pseudo_dir'], '/tmp/pseudos')
        self.assertEqual(settings['outdir'], '/tmp/qe')

    def test_system_settings_basic(self):
        settings = build_system_settings(
            cutoff_ev=340,
            nat=2,
            ntyp=2,
        )

        self.assertEqual(settings['ibrav'], 0)
        self.assertEqual(settings['nat'], 2)
        self.assertEqual(settings['ntyp'], 2)
        self.assertAlmostEqual(
            settings['ecutwfc'],
            340 / 13.605693122994,
        )

        self.assertNotIn('tot_charge', settings)
        self.assertNotIn('nbnd', settings)
        self.assertNotIn('nspin', settings)
        self.assertNotIn('ecutrho', settings)

    def test_system_settings_optional_values(self):
        settings = build_system_settings(
            cutoff_ev=400,
            nat=4,
            ntyp=2,
            total_charge=1.0,
            nbands=20,
            spinpol=True,
        )

        self.assertEqual(settings['tot_charge'], 1.0)
        self.assertEqual(settings['nbnd'], 20)
        self.assertEqual(settings['nspin'], 2)

    def test_qe_engine_can_be_loaded(self):
        engine = load_engine_module('qe')

        self.assertEqual(
            engine.QE_REFERENCE_VERSION,
            (7, 2),
        )
    
    def test_cell_parameters_are_built_in_angstrom(self):
        atoms = Atoms(
            'GaAs',
            positions=[
                [0.0, 0.0, 0.0],
                [1.0, 1.0, 1.0],
            ],
            cell=[
                [5.65, 0.0, 0.0],
                [0.0, 5.65, 0.0],
                [0.0, 0.0, 5.65],
            ],
            pbc=True,
        )

        card = build_cell_parameters(atoms)

        self.assertEqual(card['option'], 'angstrom')
        self.assertEqual(
            card['vectors'],
            [
                (5.65, 0.0, 0.0),
                (0.0, 5.65, 0.0),
                (0.0, 0.0, 5.65),
            ],
        )
    
    def test_atomic_positions_are_built_in_angstrom(self):
        atoms = Atoms(
            'GaAs',
            positions=[
                [0.0, 0.0, 0.0],
                [1.4125, 1.4125, 1.4125],
            ],
            cell=[5.65, 5.65, 5.65],
            pbc=True,
        )

        card = build_atomic_positions(atoms)

        self.assertEqual(card['option'], 'angstrom')
        self.assertEqual(
            card['positions'],
            [
                ('Ga', 0.0, 0.0, 0.0),
                ('As', 1.4125, 1.4125, 1.4125),
            ],
        )
    
    def test_atomic_species_use_pseudopotential_mapping(self):
        atoms = Atoms(
            'GaAsGa',
            positions=[
                [0.0, 0.0, 0.0],
                [1.0, 1.0, 1.0],
                [2.0, 2.0, 2.0],
            ],
        )

        species = build_atomic_species(
            atoms,
            {
                'Ga': 'Ga.upf',
                'As': 'As.upf',
            },
        )

        self.assertEqual(len(species), 2)

        self.assertEqual(species[0][0], 'Ga')
        self.assertEqual(species[0][2], 'Ga.upf')

        self.assertEqual(species[1][0], 'As')
        self.assertEqual(species[1][2], 'As.upf')

        self.assertGreater(species[0][1], 0.0)
        self.assertGreater(species[1][1], 0.0)
    
    def test_atomic_species_reject_missing_pseudopotential(self):
        atoms = Atoms(
            'GaAs',
            positions=[
                [0.0, 0.0, 0.0],
                [1.0, 1.0, 1.0],
            ],
        )

        with self.assertRaisesRegex(
            ValueError,
            'Missing pseudopotential mapping for: As',
        ):
            build_atomic_species(
                atoms,
                {
                    'Ga': 'Ga.upf',
                },
            )
    
    def test_gamma_centered_kpoint_mesh(self):
        settings = build_kpoint_settings(
            (4, 4, 4),
            gamma=True,
        )

        self.assertEqual(settings['option'], 'automatic')
        self.assertEqual(settings['size'], (4, 4, 4))
        self.assertEqual(settings['shift'], (0, 0, 0))
    
    def test_monkhorst_pack_even_mesh_is_shifted(self):
        settings = build_kpoint_settings(
            (4, 6, 8),
            gamma=False,
        )

        self.assertEqual(settings['size'], (4, 6, 8))
        self.assertEqual(settings['shift'], (1, 1, 1))
    
    def test_monkhorst_pack_shift_depends_on_mesh_parity(self):
        settings = build_kpoint_settings(
            (4, 5, 6),
            gamma=False,
        )

        self.assertEqual(settings['shift'], (1, 0, 1))
    
    def test_kpoint_mesh_rejects_nonpositive_values(self):
        with self.assertRaises(ValueError):
            build_kpoint_settings((4, 0, 4))

    def test_band_path_builds_explicit_crystal_kpoints(self):
        atoms = Atoms(
            'Si',
            cell=[4.0, 4.0, 4.0],
            pbc=True,
        )

        settings = build_band_path(
            atoms=atoms,
            path='GX',
            npoints=5,
        )

        self.assertEqual(settings['option'], 'crystal')
        self.assertEqual(settings['path'], 'GX')
        self.assertEqual(settings['npoints'], 5)
        self.assertEqual(len(settings['kpoints']), 5)
        self.assertEqual(len(settings['distances']), 5)
        self.assertEqual(settings['labels'], ['G', 'X'])

        self.assertEqual(
            settings['kpoints'][0],
            (0.0, 0.0, 0.0),
        )

        for actual, expected in zip(
            settings['kpoints'][-1],
            (0.0, 0.5, 0.0),
        ):
            self.assertAlmostEqual(actual, expected)

        self.assertAlmostEqual(
            settings['distances'][0],
            settings['special_distances'][0],
        )
        self.assertAlmostEqual(
            settings['distances'][-1],
            settings['special_distances'][-1],
        )

    def test_band_path_rejects_too_few_points(self):
        atoms = Atoms(
            'Si',
            cell=[4.0, 4.0, 4.0],
            pbc=True,
        )

        with self.assertRaisesRegex(
            ValueError,
            'at least 2 k-points',
        ):
            build_band_path(
                atoms=atoms,
                path='GX',
                npoints=1,
            )

    def test_render_explicit_band_kpoints(self):
        rendered = render_band_kpoints(
            {
                'option': 'crystal',
                'kpoints': [
                    (0.0, 0.0, 0.0),
                    (0.0, 0.25, 0.0),
                    (0.0, 0.5, 0.0),
                ],
                'npoints': 3,
            }
        )

        self.assertEqual(
            rendered,
            (
                "K_POINTS crystal\n"
                "3\n"
                "0.000000000000 0.000000000000 "
                "0.000000000000 1.0\n"
                "0.000000000000 0.250000000000 "
                "0.000000000000 1.0\n"
                "0.000000000000 0.500000000000 "
                "0.000000000000 1.0"
            ),
        )

    def test_render_band_kpoints_rejects_count_mismatch(self):
        with self.assertRaisesRegex(
            ValueError,
            'point count does not match',
        ):
            render_band_kpoints(
                {
                    'option': 'crystal',
                    'kpoints': [
                        (0.0, 0.0, 0.0),
                        (0.0, 0.5, 0.0),
                    ],
                    'npoints': 3,
                }
            )

    def test_fixed_occupation_settings(self):
        settings = build_occupation_settings('fixed')

        self.assertEqual(
            settings,
            {'occupations': 'fixed'},
        )
    
    def test_smearing_settings_convert_width_to_rydberg(self):
        settings = build_occupation_settings(
            occupations='smearing',
            smearing='fermi-dirac',
            width_ev=0.1,
        )

        self.assertEqual(
            settings['occupations'],
            'smearing',
        )
        self.assertEqual(
            settings['smearing'],
            'fermi-dirac',
        )
        self.assertAlmostEqual(
            settings['degauss'],
            0.1 / 13.605693122994,
        )
    
    def test_cold_smearing_alias(self):
        settings = build_occupation_settings(
            occupations='smearing',
            smearing='cold',
            width_ev=0.05,
        )

        self.assertEqual(
            settings['smearing'],
            'marzari-vanderbilt',
        )
    
    def test_tetrahedra_has_no_smearing_parameters(self):
        settings = build_occupation_settings(
            'tetrahedra'
        )

        self.assertEqual(
            settings,
            {'occupations': 'tetrahedra'},
        )
    
    def test_smearing_requires_width(self):
        with self.assertRaises(ValueError):
            build_occupation_settings(
                occupations='smearing',
                smearing='gaussian',
            )
    
    def test_electrons_settings_can_use_qe_defaults(self):
        settings = build_electrons_settings()

        self.assertEqual(settings, {})

    def test_electrons_settings_explicit_values(self):
        settings = build_electrons_settings(
            conv_thr=1.0e-8,
            mixing_beta=0.3,
            electron_maxstep=200,
            diagonalization='david',
        )

        self.assertEqual(settings['conv_thr'], 1.0e-8)
        self.assertEqual(settings['mixing_beta'], 0.3)
        self.assertEqual(settings['electron_maxstep'], 200)
        self.assertEqual(settings['diagonalization'], 'david')

    def test_resolve_qe_cell_dofree(self):
        cases = [
            (
                [
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                ],
                None,
            ),
            (
                [
                    True,
                    True,
                    False,
                    False,
                    False,
                    False,
                ],
                'xy',
            ),
            (
                [
                    True,
                    True,
                    True,
                    False,
                    False,
                    False,
                ],
                'xyz',
            ),
            (
                [
                    True,
                    True,
                    False,
                    False,
                    False,
                    True,
                ],
                '2Dxy',
            ),
            (
                [
                    True,
                    True,
                    True,
                    True,
                    True,
                    True,
                ],
                'all',
            ),
        ]

        for mask, expected in cases:
            with self.subTest(
                mask=mask
            ):
                self.assertEqual(
                    resolve_qe_cell_dofree(
                        mask
                    ),
                    expected,
                )

    def test_resolve_qe_atomic_relaxation_settings(self):
        settings = resolve_qe_relaxation_settings(
            optimizer='LBFGS',
            max_force=0.05,
            max_step=0.2,
            relax_cell=[
                False,
                False,
                False,
                False,
                False,
                False,
            ],
            fix_symmetry=False,
        )

        self.assertEqual(
            settings['calculation'],
            'relax',
        )

        self.assertIsNone(
            settings['cell']
        )

        self.assertEqual(
            settings['ions']['ion_dynamics'],
            'bfgs',
        )

        self.assertTrue(
            settings['system']['nosym']
        )

        self.assertAlmostEqual(
            settings['control']['forc_conv_thr'],
            ev_to_rydberg(
                0.05
                * Bohr
            ),
        )

        self.assertAlmostEqual(
            settings['ions']['trust_radius_max'],
            0.2 / Bohr,
        )

    def test_resolve_qe_variable_cell_settings(self):
        settings = resolve_qe_relaxation_settings(
            optimizer='QuasiNewton',
            max_force=0.03,
            max_step=0.1,
            relax_cell=[
                True,
                True,
                True,
                False,
                False,
                False,
            ],
            hydrostatic_pressure=2.5,
            fix_symmetry=True,
        )

        self.assertEqual(
            settings['calculation'],
            'vc-relax',
        )

        self.assertFalse(
            settings['system']['nosym']
        )

        self.assertEqual(
            settings['cell']['cell_dynamics'],
            'bfgs',
        )

        self.assertEqual(
            settings['cell']['cell_dofree'],
            'xyz',
        )

        self.assertAlmostEqual(
            settings['cell']['press'],
            25.0,
        )

    def test_qe_relaxation_rejects_unsupported_optimizer(self):
        with self.assertRaisesRegex(
            NotImplementedError,
            'QuasiNewton and LBFGS',
        ):
            resolve_qe_relaxation_settings(
                optimizer='FIRE',
                max_force=0.05,
                max_step=0.2,
                relax_cell=[
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                ],
            )

    def test_qe_value_formatter(self):
        self.assertEqual(
            format_qe_value('scf'),
            "'scf'",
        )

        self.assertEqual(
            format_qe_value(True),
            '.true.',
        )

        self.assertEqual(
            format_qe_value(False),
            '.false.',
        )

        self.assertEqual(
            format_qe_value(4),
            '4',
        )

    def test_render_namelist(self):
        text = render_namelist(
            'control',
            {
                'calculation': 'scf',
                'prefix': 'GaAs',
            },
        )

        self.assertEqual(
            text,
            "\n".join([
                "&CONTROL",
                "  calculation = 'scf',",
                "  prefix = 'GaAs',",
                "/",
            ]),
        )

    def test_render_complete_scf_input(self):
        atoms = Atoms(
            'GaAs',
            positions=[
                [0.0, 0.0, 0.0],
                [1.4125, 1.4125, 1.4125],
            ],
            cell=[
                [2.825, 2.825, 0.0],
                [2.825, 0.0, 2.825],
                [0.0, 2.825, 2.825],
            ],
            pbc=True,
        )

        text = render_scf_input(
            atoms=atoms,
            pseudopotentials={
                'Ga': 'Ga.upf',
                'As': 'As.upf',
            },
            cutoff_ev=340,
            kpoint_size=(4, 4, 4),
            gamma=False,
            occupations='fixed',
            prefix='GaAs',
        )

        self.assertIn("&CONTROL", text)
        self.assertIn("&SYSTEM", text)
        self.assertIn("&ELECTRONS", text)

        self.assertIn(
            "calculation = 'scf'",
            text,
        )

        self.assertIn(
            "prefix = 'GaAs'",
            text,
        )

        self.assertIn(
            "ibrav = 0",
            text,
        )

        self.assertIn(
            "nat = 2",
            text,
        )

        self.assertIn(
            "ntyp = 2",
            text,
        )

        self.assertIn(
            "occupations = 'fixed'",
            text,
        )

        self.assertIn(
            "ATOMIC_SPECIES",
            text,
        )

        self.assertIn(
            "Ga.upf",
            text,
        )

        self.assertIn(
            "As.upf",
            text,
        )

        self.assertIn(
            "ATOMIC_POSITIONS angstrom",
            text,
        )

        self.assertIn(
            "K_POINTS automatic",
            text,
        )

        self.assertIn(
            "4 4 4 1 1 1",
            text,
        )

        self.assertIn(
            "CELL_PARAMETERS angstrom",
            text,
        )

    def test_render_scf_input_with_smearing(self):
        atoms = Atoms(
            'Al',
            positions=[[0.0, 0.0, 0.0]],
            cell=[4.05, 4.05, 4.05],
            pbc=True,
        )

        text = render_scf_input(
            atoms=atoms,
            pseudopotentials={
                'Al': 'Al.upf',
            },
            cutoff_ev=400,
            kpoint_size=(6, 6, 6),
            occupations='smearing',
            smearing='cold',
            width_ev=0.1,
        )

        self.assertIn(
            "occupations = 'smearing'",
            text,
        )

        self.assertIn(
            "smearing = 'marzari-vanderbilt'",
            text,
        )

        self.assertIn(
            "degauss = ",
            text,
        )

    def test_render_hybrid_scf_with_additional_kpoints(self):
        atoms = bulk(
            'Si',
            'diamond',
            a=5.43,
        )

        band_path = build_band_path(
            atoms=atoms,
            path='GX',
            npoints=2,
        )
        additional_kpoints = (
            build_qe_exx_additional_kpoints(
                band_path=band_path,
                qpoint_grid=(2, 2, 2),
            )
        )

        text = render_scf_input(
            atoms=atoms,
            pseudopotentials={
                'Si': 'Si.upf',
            },
            cutoff_ev=400.0,
            kpoint_size=(2, 2, 2),
            xc_calc='HSE06',
            occupations='fixed',
            exx_additional_kpoints=(
                additional_kpoints
            ),
        )

        self.assertIn(
            "input_dft = 'HSE'",
            text,
        )
        self.assertIn(
            'nqx1 = 2',
            text,
        )
        self.assertIn(
            'nqx2 = 2',
            text,
        )
        self.assertIn(
            'nqx3 = 2',
            text,
        )
        self.assertIn(
            'K_POINTS automatic',
            text,
        )
        self.assertIn(
            'ADDITIONAL_K_POINTS crystal',
            text,
        )

    def test_render_additional_kpoints_requires_hybrid(self):
        atoms = bulk(
            'Si',
            'diamond',
            a=5.43,
        )
        additional_kpoints = (
            build_qe_exx_additional_kpoints(
                band_path={
                    'kpoints': [
                        (0.0, 0.0, 0.0),
                    ],
                },
                qpoint_grid=(1, 1, 1),
            )
        )

        with self.assertRaisesRegex(
            ValueError,
            'require a hybrid functional',
        ):
            render_scf_input(
                atoms=atoms,
                pseudopotentials={
                    'Si': 'Si.upf',
                },
                cutoff_ev=400.0,
                kpoint_size=(2, 2, 2),
                xc_calc='PBE',
                exx_additional_kpoints=(
                    additional_kpoints
                ),
            )

    def test_render_dos_input_for_tetrahedra(self):
        text = render_dos_input(
            prefix='nanoworks',
            outdir='/tmp/qe-state',
            fildos='/tmp/gaas.dos',
            bz_sum='tetrahedra',
            emin=-10.0,
            emax=10.0,
            delta_e=0.02,
        )

        self.assertIn(
            '&DOS',
            text,
        )

        self.assertIn(
            "prefix = 'nanoworks'",
            text,
        )

        self.assertIn(
            "outdir = '/tmp/qe-state'",
            text,
        )

        self.assertIn(
            "bz_sum = 'tetrahedra'",
            text,
        )

        self.assertIn(
            'Emin = -10',
            text,
        )

        self.assertIn(
            'Emax = 10',
            text,
        )

        self.assertIn(
            'DeltaE = 0.02',
            text,
        )

        self.assertIn(
            "fildos = '/tmp/gaas.dos'",
            text,
        )

        self.assertNotIn(
            'degauss',
            text,
        )

    def test_render_dos_input_rejects_invalid_energy_range(self):
        with self.assertRaisesRegex(
            ValueError,
            'Emax must be greater than Emin',
        ):
            render_dos_input(
                emin=5.0,
                emax=-5.0,
            )

    def test_render_epsilon_input(self):
        text = render_epsilon_input(
            prefix='nanoworks',
            outdir='/tmp/qe-optical-state',
            calculation='eps',
            smeartype='gauss',
            intersmear=0.1,
            intrasmear=0.0,
            wmin=0.0,
            wmax=20.0,
            nw=401,
            nbndmin=1,
            nbndmax=16,
            shift=0.25,
        )

        self.assertIn(
            '&INPUTPP',
            text,
        )
        self.assertIn(
            "prefix = 'nanoworks'",
            text,
        )
        self.assertIn(
            "outdir = '/tmp/qe-optical-state'",
            text,
        )
        self.assertIn(
            "calculation = 'eps'",
            text,
        )
        self.assertIn(
            '&ENERGY_GRID',
            text,
        )
        self.assertIn(
            '/\n&ENERGY_GRID',
            text,
        )
        self.assertIn(
            "smeartype = 'gauss'",
            text,
        )
        self.assertIn(
            'intersmear = 0.1',
            text,
        )
        self.assertIn(
            'intrasmear = 0',
            text,
        )
        self.assertIn(
            'wmin = 0',
            text,
        )
        self.assertIn(
            'wmax = 20',
            text,
        )
        self.assertIn(
            'nw = 401',
            text,
        )
        self.assertIn(
            'nbndmin = 1',
            text,
        )
        self.assertIn(
            'nbndmax = 16',
            text,
        )
        self.assertIn(
            'shift = 0.25',
            text,
        )

    def test_render_epsilon_input_rejects_invalid_calculation(self):
        with self.assertRaisesRegex(
            ValueError,
            'Unsupported QE epsilon.x calculation',
        ):
            render_epsilon_input(
                calculation='occ',
            )

    def test_render_epsilon_input_rejects_invalid_frequency_grid(self):
        with self.assertRaisesRegex(
            ValueError,
            'wmax must be greater than wmin',
        ):
            render_epsilon_input(
                wmin=10.0,
                wmax=5.0,
            )

        with self.assertRaisesRegex(
            ValueError,
            'at least two points',
        ):
            render_epsilon_input(
                nw=1,
            )

    def test_render_epsilon_input_rejects_invalid_band_range(self):
        with self.assertRaisesRegex(
            ValueError,
            'nbndmax must not be smaller than nbndmin',
        ):
            render_epsilon_input(
                nbndmin=8,
                nbndmax=4,
            )

    def test_run_epsilon_requires_qe_state(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(
                tmpdir
            )

            with self.assertRaisesRegex(
                FileNotFoundError,
                'valid QE electronic state',
            ):
                run_epsilon(
                    input_file=tmpdir / 'epsilon.in',
                    output_file=tmpdir / 'epsilon.out',
                    state_dir=tmpdir / 'state',
                    result_dir=tmpdir / 'optical',
                )

    def test_run_epsilon_validates_dielectric_outputs(self):
        output_text = """
         Program epsilon v.7.6 starts
         JOB DONE.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(
                tmpdir
            )
            state_dir = tmpdir / 'state'
            save_dir = state_dir / 'nanoworks.save'
            save_dir.mkdir(
                parents=True
            )
            (
                save_dir
                / 'data-file-schema.xml'
            ).write_text(
                '<qes/>',
                encoding='utf-8',
            )
            result_dir = tmpdir / 'optical'

            def fake_run_qe_program(**kwargs):
                Path(
                    kwargs['output_file']
                ).write_text(
                    output_text,
                    encoding='utf-8',
                )

                for name in (
                    'epsr.dat',
                    'epsi.dat',
                    'eels.dat',
                    'ieps.dat',
                ):
                    (
                        Path(kwargs['cwd'])
                        / name
                    ).write_text(
                        '# epsilon data\n',
                        encoding='utf-8',
                    )

                return {
                    'returncode': 0,
                }

            with patch(
                'nanoworks.engine.qe.run_qe_program',
                side_effect=fake_run_qe_program,
            ) as run:
                workflow = run_epsilon(
                    input_file=tmpdir / 'epsilon.in',
                    output_file=tmpdir / 'epsilon.out',
                    state_dir=state_dir,
                    result_dir=result_dir,
                    wmax=20.0,
                    nw=401,
                )

            run.assert_called_once()
            self.assertEqual(
                run.call_args.kwargs['cwd'],
                result_dir.resolve(),
            )
            self.assertEqual(
                workflow['metadata']['program'],
                'EPSILON',
            )
            self.assertTrue(
                workflow['metadata']['job_done']
            )
            self.assertEqual(
                set(workflow['result_files']),
                {
                    'epsr.dat',
                    'epsi.dat',
                    'eels.dat',
                    'ieps.dat',
                },
            )
            self.assertIn(
                'wmax = 20',
                (
                    tmpdir
                    / 'epsilon.in'
                ).read_text(
                    encoding='utf-8',
                ),
            )

    def test_run_epsilon_rejects_missing_data_files(self):
        output_text = """
         Program epsilon v.7.6 starts
         JOB DONE.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(
                tmpdir
            )
            state_dir = tmpdir / 'state'
            save_dir = state_dir / 'nanoworks.save'
            save_dir.mkdir(
                parents=True
            )
            (
                save_dir
                / 'data-file-schema.xml'
            ).write_text(
                '<qes/>',
                encoding='utf-8',
            )

            def fake_run_qe_program(**kwargs):
                Path(
                    kwargs['output_file']
                ).write_text(
                    output_text,
                    encoding='utf-8',
                )

                return {
                    'returncode': 0,
                }

            with (
                patch(
                    'nanoworks.engine.qe.run_qe_program',
                    side_effect=fake_run_qe_program,
                ),
                self.assertRaisesRegex(
                    RuntimeError,
                    'expected data files were not created',
                ),
            ):
                run_epsilon(
                    input_file=tmpdir / 'epsilon.in',
                    output_file=tmpdir / 'epsilon.out',
                    state_dir=state_dir,
                    result_dir=tmpdir / 'optical',
                )

    def test_parse_epsilon_data_file(self):
        content = """# energy grid [eV] epsr_x epsr_y epsr_z
# plasmon frequencies [eV]
 0.000000000  1.0D+00  2.0D+00  3.0D+00
 0.500000000  1.1D+00  2.1D+00  3.1D+00
 1.000000000  1.2D+00  2.2D+00  3.2D+00
"""

        with tempfile.TemporaryDirectory() as tmpdir:
            data_file = Path(tmpdir) / 'epsr.dat'
            data_file.write_text(
                content,
                encoding='utf-8',
            )

            result = parse_epsilon_data_file(
                data_file,
                expected_components=3,
            )

        self.assertEqual(
            result['energies_ev'],
            [0.0, 0.5, 1.0],
        )
        self.assertEqual(
            result['components'],
            [
                [1.0, 1.1, 1.2],
                [2.0, 2.1, 2.2],
                [3.0, 3.1, 3.2],
            ],
        )
        self.assertEqual(
            result['npoints'],
            3,
        )
        self.assertEqual(
            result['ncomponents'],
            3,
        )

    def test_parse_epsilon_data_file_rejects_inconsistent_columns(self):
        content = """# epsilon data
 0.0  1.0  2.0  3.0
 1.0  1.1  2.1
"""

        with tempfile.TemporaryDirectory() as tmpdir:
            data_file = Path(tmpdir) / 'epsr.dat'
            data_file.write_text(
                content,
                encoding='utf-8',
            )

            with self.assertRaisesRegex(
                ValueError,
                'inconsistent column counts',
            ):
                parse_epsilon_data_file(
                    data_file
                )

    def test_parse_epsilon_data_file_validates_grid_and_components(self):
        content = """# epsilon data
 0.0  1.0  2.0  3.0
 0.0  1.1  2.1  3.1
"""

        with tempfile.TemporaryDirectory() as tmpdir:
            data_file = Path(tmpdir) / 'epsr.dat'
            data_file.write_text(
                content,
                encoding='utf-8',
            )

            with self.assertRaisesRegex(
                ValueError,
                'component count does not match',
            ):
                parse_epsilon_data_file(
                    data_file,
                    expected_components=2,
                )

            with self.assertRaisesRegex(
                ValueError,
                'strictly increasing',
            ):
                parse_epsilon_data_file(
                    data_file,
                    expected_components=3,
                )

    def test_render_dos_input_rejects_invalid_bz_sum(self):
        with self.assertRaisesRegex(
            ValueError,
            'Unsupported QE DOS BZ summation method',
        ):
            render_dos_input(
                bz_sum='invalid',
            )

    def test_run_dos_requires_qe_state(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(
                tmpdir
            )

            with self.assertRaisesRegex(
                FileNotFoundError,
                'valid QE electronic state',
            ):
                run_dos(
                    input_file=tmpdir / 'dos.in',
                    output_file=tmpdir / 'dos.out',
                    state_dir=tmpdir / 'state',
                    dos_file=tmpdir / 'result.dos',
                )

    def test_render_nscf_input(self):
        atoms = bulk(
            'Si',
            'diamond',
            a=5.43,
        )

        text = render_nscf_input(
            atoms=atoms,
            pseudopotentials={
                'Si': 'Si.upf',
            },
            cutoff_ev=400.0,
            kpoint_size=(8, 8, 8),
            occupations='tetrahedra',
            prefix='nanoworks',
            pseudo_dir='/tmp/pseudos',
            outdir='/tmp/qe-state',
        )

        self.assertIn(
            "calculation = 'nscf'",
            text,
        )

        self.assertIn(
            "occupations = 'tetrahedra'",
            text,
        )

        self.assertIn(
            "8 8 8 1 1 1",
            text,
        )

        self.assertIn(
            "prefix = 'nanoworks'",
            text,
        )

        self.assertIn(
            "outdir = '/tmp/qe-state'",
            text,
        )

    def test_render_complete_bands_input(self):
        atoms = bulk(
            'Si',
            'diamond',
            a=5.43,
        )

        band_path = build_band_path(
            atoms=atoms,
            path='GX',
            npoints=5,
        )

        text = render_pw_input(
            calculation='bands',
            atoms=atoms,
            pseudopotentials={
                'Si': 'Si.upf',
            },
            cutoff_ev=400.0,
            kpoint_size=None,
            occupations='fixed',
            prefix='nanoworks',
            outdir='/tmp/qe-state',
            band_path=band_path,
        )

        self.assertIn(
            "calculation = 'bands'",
            text,
        )

        self.assertIn(
            "verbosity = 'high'",
            text,
        )

        self.assertIn(
            "K_POINTS crystal",
            text,
        )

        self.assertIn(
            (
                "0.000000000000 0.000000000000 "
                "0.000000000000 1.0"
            ),
            text,
        )

        self.assertIn(
            (
                "0.500000000000 0.000000000000 "
                "0.500000000000 1.0"
            ),
            text,
        )

        self.assertNotIn(
            "K_POINTS automatic",
            text,
        )

    def test_render_bands_input_wrapper(self):
        atoms = bulk(
            'Si',
            'diamond',
            a=5.43,
        )

        band_path = build_band_path(
            atoms=atoms,
            path='GX',
            npoints=5,
        )

        text = render_bands_input(
            atoms=atoms,
            pseudopotentials={
                'Si': 'Si.upf',
            },
            cutoff_ev=400.0,
            band_path=band_path,
            nbands=12,
            occupations='fixed',
            prefix='nanoworks',
            pseudo_dir='/tmp/pseudos',
            outdir='/tmp/qe-state',
        )

        self.assertIn(
            "calculation = 'bands'",
            text,
        )

        self.assertIn(
            "nbnd = 12",
            text,
        )

        self.assertIn(
            "K_POINTS crystal",
            text,
        )

        self.assertIn(
            "\n5\n",
            text,
        )

        self.assertIn(
            "pseudo_dir = '/tmp/pseudos'",
            text,
        )

        self.assertIn(
            "outdir = '/tmp/qe-state'",
            text,
        )

        self.assertNotIn(
            "K_POINTS automatic",
            text,
        )

    def test_render_bands_postprocess_input(self):
        text = render_bands_postprocess_input(
            prefix='nanoworks',
            outdir='/tmp/qe-band-state',
            filband='/tmp/si-hse.bands',
        )

        self.assertIn(
            '&BANDS',
            text,
        )
        self.assertIn(
            "prefix = 'nanoworks'",
            text,
        )
        self.assertIn(
            "outdir = '/tmp/qe-band-state'",
            text,
        )
        self.assertIn(
            "filband = '/tmp/si-hse.bands'",
            text,
        )
        self.assertIn(
            'lsym = .false.',
            text,
        )

    def test_render_bands_input_requires_band_path(self):
        atoms = bulk(
            'Si',
            'diamond',
            a=5.43,
        )

        with self.assertRaisesRegex(
            ValueError,
            'requires an explicit band path',
        ):
            render_pw_input(
                calculation='bands',
                atoms=atoms,
                pseudopotentials={
                    'Si': 'Si.upf',
                },
                cutoff_ev=400.0,
                kpoint_size=None,
            )

    def test_render_scf_input_rejects_band_path(self):
        atoms = bulk(
            'Si',
            'diamond',
            a=5.43,
        )

        band_path = build_band_path(
            atoms=atoms,
            path='GX',
            npoints=5,
        )

        with self.assertRaisesRegex(
            ValueError,
            'can only be used',
        ):
            render_pw_input(
                calculation='scf',
                atoms=atoms,
                pseudopotentials={
                    'Si': 'Si.upf',
                },
                cutoff_ev=400.0,
                kpoint_size=(4, 4, 4),
                band_path=band_path,
            )

    def test_render_pw_input_rejects_unsupported_calculation(self):
        atoms = bulk(
            'Si',
            'diamond',
            a=5.43,
        )

        with self.assertRaisesRegex(
            ValueError,
            'Unsupported QE pw.x calculation type',
        ):
            render_pw_input(
                calculation='md',
                atoms=atoms,
                pseudopotentials={
                    'Si': 'Si.upf',
                },
                cutoff_ev=400.0,
                kpoint_size=(4, 4, 4),
            )

    def test_render_atomic_relax_input(self):
        atoms = bulk(
            'Si',
            'diamond',
            a=5.43,
        )

        text = render_relax_input(
            atoms=atoms,
            pseudopotentials={
                'Si': 'Si.upf',
            },
            cutoff_ev=400.0,
            kpoint_size=(4, 4, 4),
            optimizer='LBFGS',
            max_force=0.05,
            max_step=0.20,
            relax_cell=[
                False,
                False,
                False,
                False,
                False,
                False,
            ],
            hydrostatic_pressure=0.0,
            fix_symmetry=False,
            prefix='nanoworks',
            pseudo_dir='/tmp/pseudos',
            outdir='/tmp/qe-state',
        )

        self.assertIn(
            "calculation = 'relax'",
            text,
        )
        self.assertIn(
            'forc_conv_thr =',
            text,
        )
        self.assertIn(
            'nosym = .true.',
            text,
        )
        self.assertIn(
            '&IONS',
            text,
        )
        self.assertIn(
            "ion_dynamics = 'bfgs'",
            text,
        )
        self.assertIn(
            'trust_radius_max =',
            text,
        )
        self.assertIn(
            'trust_radius_ini =',
            text,
        )
        self.assertNotIn(
            '&CELL',
            text,
        )

    def test_render_variable_cell_relax_input(self):
        atoms = bulk(
            'Si',
            'diamond',
            a=5.43,
        )

        text = render_relax_input(
            atoms=atoms,
            pseudopotentials={
                'Si': 'Si.upf',
            },
            cutoff_ev=400.0,
            kpoint_size=(4, 4, 4),
            optimizer='QuasiNewton',
            max_force=0.05,
            max_step=0.20,
            relax_cell=[
                True,
                True,
                True,
                False,
                False,
                False,
            ],
            hydrostatic_pressure=2.5,
            fix_symmetry=True,
            prefix='nanoworks',
            pseudo_dir='/tmp/pseudos',
            outdir='/tmp/qe-state',
        )

        self.assertIn(
            "calculation = 'vc-relax'",
            text,
        )
        self.assertIn(
            'nosym = .false.',
            text,
        )
        self.assertIn(
            '&IONS',
            text,
        )
        self.assertIn(
            '&CELL',
            text,
        )
        self.assertIn(
            "cell_dynamics = 'bfgs'",
            text,
        )
        self.assertIn(
            "cell_dofree = 'all'",
            text,
        )
        
        self.assertIn(
            "! NOTICE: The normal-strain Relax_cell mask",
            text,
        )
        
        self.assertIn(
            'press = 25',
            text,
        )
        self.assertLess(
            text.index('&IONS'),
            text.index('&CELL'),
        )
        self.assertLess(
            text.index('&CELL'),
            text.index('ATOMIC_SPECIES'),
        )

    def test_render_orthogonal_cell_keeps_xyz_relaxation(self):
        atoms = Atoms(
            'Si',
            positions=[
                [
                    0.0,
                    0.0,
                    0.0,
                ],
            ],
            cell=[
                [
                    5.0,
                    0.0,
                    0.0,
                ],
                [
                    0.0,
                    5.0,
                    0.0,
                ],
                [
                    0.0,
                    0.0,
                    5.0,
                ],
            ],
            pbc=True,
        )

        text = render_relax_input(
            atoms=atoms,
            pseudopotentials={
                'Si': 'Si.upf',
            },
            cutoff_ev=400.0,
            kpoint_size=(4, 4, 4),
            optimizer='LBFGS',
            max_force=0.05,
            max_step=0.20,
            relax_cell=[
                True,
                True,
                True,
                False,
                False,
                False,
            ],
            fix_symmetry=True,
        )

        self.assertIn(
            "cell_dofree = 'xyz'",
            text,
        )
        self.assertNotIn(
            '! NOTICE:',
            text,
        )

    def test_render_pw_input_rejects_mismatched_relaxation_type(self):
        atoms = bulk(
            'Si',
            'diamond',
            a=5.43,
        )

        settings = resolve_qe_relaxation_settings(
            optimizer='LBFGS',
            max_force=0.05,
            max_step=0.20,
            relax_cell=[
                False,
                False,
                False,
                False,
                False,
                False,
            ],
        )

        with self.assertRaisesRegex(
            ValueError,
            'does not match',
        ):
            render_pw_input(
                calculation='vc-relax',
                atoms=atoms,
                pseudopotentials={
                    'Si': 'Si.upf',
                },
                cutoff_ev=400.0,
                kpoint_size=(4, 4, 4),
                relaxation_settings=settings,
            )
        
    def test_rydberg_to_ev(self):
        self.assertAlmostEqual(
            rydberg_to_ev(1.0),
            13.605693122994,
        )

    def test_has_qe_state_accepts_complete_saved_state(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            state_dir = Path(tmpdir)

            save_dir = (
                state_dir
                / 'nanoworks.save'
            )

            save_dir.mkdir()

            (
                save_dir
                / 'data-file-schema.xml'
            ).write_text(
                '<espresso/>',
                encoding='utf-8',
            )

            self.assertTrue(
                has_qe_state(
                    state_dir
                )
            )

    def test_has_qe_state_rejects_incomplete_saved_state(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            state_dir = Path(tmpdir)

            (
                state_dir
                / 'nanoworks.save'
            ).mkdir()

            self.assertFalse(
                has_qe_state(
                    state_dir
                )
            )
    
    def test_has_qe_state_rejects_missing_state_directory(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            state_dir = (
                Path(tmpdir)
                / 'missing-state'
            )

            self.assertFalse(
                has_qe_state(
                    state_dir
                )
            )
            
    @patch('nanoworks.engine.qe.shutil.which')
    def test_resolve_qe_executable_from_path(
        self,
        which,
    ):
        which.return_value = '/opt/qe/bin/pw.x'

        resolved = resolve_qe_executable(
            'pw.x'
        )

        self.assertEqual(
            resolved,
            '/opt/qe/bin/pw.x',
        )

        which.assert_called_once_with(
            'pw.x'
        )

    @patch('nanoworks.engine.qe.resolve_qe_executable')
    def test_build_qe_command_with_mpi(
        self,
        resolve,
    ):
        resolve.return_value = '/opt/qe/bin/pw.x'

        command = build_qe_command(
            input_file='si.in',
            executable='pw.x',
            launcher=[
                'mpiexec',
                '-np',
                '8',
            ],
        )

        self.assertEqual(
            command,
            [
                'mpiexec',
                '-np',
                '8',
                '/opt/qe/bin/pw.x',
                '-i',
                'si.in',
            ],
        )

    @patch('nanoworks.engine.qe.resolve_qe_executable')
    def test_qe_launcher_rejects_shell_string(
        self,
        resolve,
    ):
        resolve.return_value = '/opt/qe/bin/pw.x'

        with self.assertRaises(TypeError):
            build_qe_command(
                input_file='si.in',
                launcher='mpiexec -np 8',
            )

    def test_parse_pw_output(self):
        output_text = """
         Program PWSCF v.7.2 starts

    !    total energy              =     -15.12345678 Ry
    
         the Fermi energy is     5.4321 ev

         JOB DONE.
    """

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'si.out'

            path.write_text(
                output_text
            )

            result = parse_pw_output(
                path
            )

            self.assertTrue(
                result['job_done']
            )
            
            self.assertEqual(
                result['qe_version'],
                (7, 2),
            )

            self.assertAlmostEqual(
                result['fermi_energy_ev'],
                5.4321,
            )
        
            self.assertAlmostEqual(
                result['total_energy_ry'],
                -15.12345678,
            )

            self.assertAlmostEqual(
                result['total_energy_ev'],
                -15.12345678
                * 13.605693122994,
            )
    
    def test_parse_pw_output_with_patch_version(self):
        output_text = """
         Program PWSCF v.7.2.1 starts

         JOB DONE.
        """

        with tempfile.TemporaryDirectory() as tmp:
            path = (
                Path(tmp)
                / 'pw.out'
            )

            path.write_text(
                output_text,
                encoding='utf-8',
            )

            result = parse_pw_output(
                path
            )

            self.assertEqual(
                result['qe_version'],
                (7, 2, 1),
            )
        
    def test_parse_pw_output_without_version(self):
        output_text = """
         JOB DONE.
        """

        with tempfile.TemporaryDirectory() as tmp:
            path = (
                Path(tmp)
                / 'pw.out'
            )

            path.write_text(
                output_text,
                encoding='utf-8',
            )

            result = parse_pw_output(
                path
            )

            self.assertIsNone(
                result['qe_version']
            )
            
    def test_parse_pw_output_without_fermi_energy(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = (
                Path(tmpdir)
                / 'pw.out'
            )

            output_file.write_text(
                "\n".join([
                    "!    total energy = -10.00000000 Ry",
                    "JOB DONE.",
                ]),
                encoding='utf-8',
            )

            result = parse_pw_output(
                output_file
            )

            self.assertIsNone(
                result['fermi_energy_ev']
            )

    def test_parse_pw_output_resolves_band_edges(self):
        output_text = """
         Program PWSCF v.7.2 starts

         highest occupied, lowest unoccupied level (ev):
             5.0000  6.0000

         JOB DONE.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = (
                Path(tmpdir)
                / 'pw.out'
            )

            output_file.write_text(
                output_text,
                encoding='utf-8',
            )

            result = parse_pw_output(
                output_file
            )

        self.assertAlmostEqual(
            result['highest_occupied_ev'],
            5.0,
        )

        self.assertAlmostEqual(
            result['lowest_unoccupied_ev'],
            6.0,
        )

        reference = resolve_qe_band_reference(
            result
        )

        self.assertAlmostEqual(
            reference['energy_ev'],
            5.5,
        )

        self.assertEqual(
            reference['source'],
            'midgap',
        )

    def test_qe_band_reference_prefers_fermi_energy(self):
        reference = resolve_qe_band_reference(
            {
                'fermi_energy_ev': 5.25,
                'highest_occupied_ev': 5.0,
                'lowest_unoccupied_ev': 6.0,
            }
        )

        self.assertAlmostEqual(
            reference['energy_ev'],
            5.25,
        )

        self.assertEqual(
            reference['source'],
            'fermi',
        )

    def test_qe_band_reference_uses_highest_occupied_fallback(self):
        reference = resolve_qe_band_reference(
            {
                'fermi_energy_ev': None,
                'highest_occupied_ev': 5.0,
                'lowest_unoccupied_ev': None,
            }
        )

        self.assertAlmostEqual(
            reference['energy_ev'],
            5.0,
        )

        self.assertEqual(
            reference['source'],
            'highest_occupied',
        )

    def test_parse_pw_bands_output(self):
        output_text = """
         End of band structure calculation

              k = 0.0000 0.0000 0.0000 ( 123 PWs)   bands (ev):

            -5.0000  -1.0000   1.0000
             3.0000   5.0000

              k =-0.5000-0.2887 0.4083 ( 120 PWs)   bands (ev):

            -4.5000  -0.5000   1.5000
             3.5000   5.5000

         JOB DONE.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = (
                Path(tmpdir)
                / 'bands.out'
            )

            output_file.write_text(
                output_text,
                encoding='utf-8',
            )

            result = parse_pw_bands_output(
                output_file
            )

        self.assertFalse(
            result['spin_polarized']
        )

        self.assertEqual(
            result['nspins'],
            1,
        )

        self.assertEqual(
            result['nkpoints'],
            2,
        )

        self.assertEqual(
            result['nbands'],
            5,
        )

        self.assertEqual(
            result['kpoints'],
            [
                (0.0, 0.0, 0.0),
                (-0.5, -0.2887, 0.4083),
            ],
        )

        self.assertEqual(
            result['eigenvalues_ev'],
            [[
                [-5.0, -1.0, 1.0, 3.0, 5.0],
                [-4.5, -0.5, 1.5, 3.5, 5.5],
            ]],
        )

    def test_parse_bands_x_output_selects_band_path_points(self):
        band_text = """
 &plot nbnd= 3, nks= 3 /
  0.000000  0.000000  0.000000
 -5.000000 -1.000000  1.000000
  0.250000  0.000000  0.000000
 -4.000000 -0.500000  1.500000
  0.500000  0.000000  0.000000
 -3.000000  0.000000  2.000000
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            band_file = (
                Path(tmpdir)
                / 'si-hse.bands'
            )
            band_file.write_text(
                band_text,
                encoding='utf-8',
            )

            result = parse_bands_x_output(
                band_file,
                kpoint_indices=[0, 2],
            )

        self.assertFalse(
            result['spin_polarized']
        )
        self.assertEqual(
            result['nkpoints'],
            2,
        )
        self.assertEqual(
            result['nbands'],
            3,
        )
        self.assertEqual(
            result['kpoints'],
            [
                (0.0, 0.0, 0.0),
                (0.5, 0.0, 0.0),
            ],
        )
        self.assertEqual(
            result['eigenvalues_ev'],
            [[
                [-5.0, -1.0, 1.0],
                [-3.0, 0.0, 2.0],
            ]],
        )

    def test_prepare_qe_band_data_shifts_reference(self):
        bands = {
            'spin_polarized': False,
            'nspins': 1,
            'nkpoints': 2,
            'nbands': 3,
            'eigenvalues_ev': [[
                [0.0, 1.0, 2.0],
                [0.5, 1.5, 2.5],
            ]],
        }

        band_path = {
            'distances': [
                0.0,
                0.5,
            ],
            'special_distances': [
                0.0,
                0.5,
            ],
            'labels': [
                'G',
                'X',
            ],
        }

        result = prepare_qe_band_data(
            bands=bands,
            band_path=band_path,
            reference_energy=1.0,
        )

        self.assertEqual(
            result['nkpoints'],
            2,
        )

        self.assertEqual(
            result['nbands'],
            3,
        )

        self.assertEqual(
            result['eigenvalues_ev'],
            [
                [-1.0, 0.0, 1.0],
                [-0.5, 0.5, 1.5],
            ],
        )

        self.assertEqual(
            result['distances'],
            [
                0.0,
                0.5,
            ],
        )

        self.assertEqual(
            result['labels'],
            [
                'G',
                'X',
            ],
        )

    def test_parse_pw_bands_output_rejects_inconsistent_band_counts(self):
        output_text = """
              k = 0.0000 0.0000 0.0000 ( 123 PWs)   bands (ev):

            -5.0000  -1.0000   1.0000

              k = 0.5000 0.0000 0.5000 ( 120 PWs)   bands (ev):

            -4.5000  -0.5000

         JOB DONE.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = (
                Path(tmpdir)
                / 'bands.out'
            )

            output_file.write_text(
                output_text,
                encoding='utf-8',
            )

            with self.assertRaisesRegex(
                ValueError,
                'inconsistent numbers of bands',
            ):
                parse_pw_bands_output(
                    output_file
                )

    @patch('nanoworks.engine.qe.subprocess.run')
    @patch('nanoworks.engine.qe.resolve_qe_executable')
    def test_run_qe_program(
        self,
        resolve,
        run,
    ):
        resolve.return_value = '/opt/qe/bin/pw.x'

        run.return_value.returncode = 0

        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)

            input_file = (
                tmp / 'si.in'
            )

            output_file = (
                tmp / 'si.out'
            )

            input_file.write_text(
                '&CONTROL\n/\n'
            )

            result = run_qe_program(
                input_file=input_file,
                output_file=output_file,
            )

            self.assertEqual(
                result['returncode'],
                0,
            )

            self.assertTrue(
                output_file.exists()
            )

            command = (
                run.call_args.args[0]
            )

            self.assertEqual(
                command[-2:],
                [
                    '-i',
                    str(input_file),
                ],
            )

            environment = (
                run.call_args.kwargs['env']
            )

            self.assertEqual(
                environment['OMP_NUM_THREADS'],
                '1',
            )

            self.assertEqual(
                environment['OPENBLAS_NUM_THREADS'],
                '1',
            )

            self.assertEqual(
                environment['MKL_NUM_THREADS'],
                '1',
            )

            self.assertEqual(
                environment['VECLIB_MAXIMUM_THREADS'],
                '1',
            )

            self.assertEqual(
                environment['NUMEXPR_NUM_THREADS'],
                '1',
            )

            self.assertEqual(
                environment['OMP_DYNAMIC'],
                'FALSE',
            )
    
    def test_qe_launcher_is_none_for_single_core(self):
        launcher = build_qe_launcher(
            parallel_cores=1
        )

        self.assertIsNone(
            launcher
        )
    
    @patch('nanoworks.engine.qe.shutil.which')
    def test_qe_launcher_uses_mpiexec(
        self,
        which,
    ):
        which.side_effect = (
            lambda name: (
                '/usr/bin/mpiexec'
                if name == 'mpiexec'
                else None
            )
        )

        launcher = build_qe_launcher(
            parallel_cores=8
        )

        self.assertEqual(
            launcher,
            [
                '/usr/bin/mpiexec',
                '-np',
                '8',
            ],
        )

    @patch('nanoworks.engine.qe.shutil.which')
    def test_qe_launcher_uses_srun(
        self,
        which,
    ):
        available = {
            'mpiexec': None,
            'mpirun': None,
            'srun': '/usr/bin/srun',
        }

        which.side_effect = (
            lambda name: available[name]
        )

        launcher = build_qe_launcher(
            parallel_cores=16
        )

        self.assertEqual(
            launcher,
            [
                '/usr/bin/srun',
                '-n',
                '16',
            ],
        )

    def test_qe_launcher_rejects_nonpositive_core_count(self):
        with self.assertRaises(ValueError):
            build_qe_launcher(
                parallel_cores=0
            )

    @patch('nanoworks.engine.qe.shutil.which')
    def test_qe_launcher_requires_mpi_command(
        self,
        which,
    ):
        which.return_value = None

        with self.assertRaises(
            FileNotFoundError
        ):
            build_qe_launcher(
                parallel_cores=8
            )

    def test_qe_kpoint_density_resolves_mesh(self):
        atoms = Atoms(
            'Si2',
            positions=[
                [0.0, 0.0, 0.0],
                [1.35, 1.35, 1.35],
            ],
            cell=[
                [5.4, 0.0, 0.0],
                [0.0, 5.4, 0.0],
                [0.0, 0.0, 5.4],
            ],
            pbc=True,
        )

        mesh = resolve_qe_kpoint_size(
            atoms,
            density=2.5,
            size=(1, 1, 1),
        )

        self.assertEqual(
            len(mesh),
            3,
        )

        self.assertTrue(
            all(value > 0 for value in mesh)
        )

    def test_resolve_qe_occupation_supports_tetrahedra(self):
        result = resolve_qe_occupation(
            'tetrahedra'
        )

        self.assertEqual(
            result,
            {
                'occupations': 'tetrahedra',
                'smearing': None,
                'width_ev': None,
            },
        )

    def test_resolve_qe_occupation_supports_optimized_tetrahedra(self):
        result = resolve_qe_occupation(
            'tetrahedra-opt'
        )

        self.assertEqual(
            result['occupations'],
            'tetrahedra_opt',
        )

    def test_qe_occupation_translates_fermi_dirac(self):
        settings = resolve_qe_occupation(
            {
                'name': 'fermi-dirac',
                'width': 0.05,
            }
        )

        self.assertEqual(
            settings,
            {
                'occupations': 'smearing',
                'smearing': 'fermi-dirac',
                'width_ev': 0.05,
            },
        )

    def test_qe_xc_accepts_pbe(self):
        self.assertEqual(
            validate_qe_xc('PBE'),
            'pbe',
        )

    def test_qe_xc_accepts_pbe(self):
        self.assertEqual(
            validate_qe_xc('PBE'),
            'pbe',
        )
    
    def test_validate_qe_version_accepts_reference_version(self):
        result = validate_qe_version(
            (7, 2)
        )

        self.assertEqual(
            result,
            (7, 2),
        )
    
    def test_validate_qe_version_accepts_patch_version(self):
        result = validate_qe_version(
            (7, 2, 1)
        )

        self.assertEqual(
            result,
            (7, 2, 1),
        )

    def test_validate_qe_version_rejects_older_version(self):
        with self.assertRaisesRegex(
            ValueError,
            'requires Quantum ESPRESSO 7.2 or newer',
        ):
            validate_qe_version(
                (7, 1)
            )

    def test_validate_qe_version_rejects_missing_version(self):
        with self.assertRaisesRegex(
            ValueError,
            'version could not be detected',
        ):
            validate_qe_version(
                None
            )

    def test_run_nscf_requires_ground_state(self):
        atoms = bulk(
            'Si',
            'diamond',
            a=5.43,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            with self.assertRaisesRegex(
                FileNotFoundError,
                'valid QE ground-state result',
            ):
                run_nscf(
                    atoms=atoms,
                    input_file=tmpdir / 'nscf.in',
                    output_file=tmpdir / 'nscf.out',
                    state_dir=tmpdir / 'state',
                    pseudopotentials={
                        'Si': 'Si.upf',
                    },
                    pseudo_dir='/tmp/pseudos',
                    cutoff_ev=400.0,
                )

    def test_run_bands_requires_ground_state(self):
        atoms = bulk(
            'Si',
            'diamond',
            a=5.43,
        )

        band_path = build_band_path(
            atoms=atoms,
            path='GX',
            npoints=5,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            with self.assertRaisesRegex(
                FileNotFoundError,
                'valid QE ground-state result',
            ):
                run_bands(
                    atoms=atoms,
                    input_file=tmpdir / 'bands.in',
                    output_file=tmpdir / 'bands.out',
                    state_dir=tmpdir / 'state',
                    pseudopotentials={
                        'Si': 'Si.upf',
                    },
                    pseudo_dir='/tmp/pseudos',
                    cutoff_ev=400.0,
                    band_path=band_path,
                )

    def test_run_bands_postprocess_creates_band_file(self):
        output_text = """
         Program BANDS v.7.6 starts
         JOB DONE.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(
                tmpdir
            )
            state_dir = (
                tmpdir
                / 'state'
            )
            save_dir = (
                state_dir
                / 'nanoworks.save'
            )
            save_dir.mkdir(
                parents=True
            )
            (
                save_dir
                / 'data-file-schema.xml'
            ).write_text(
                '<qes/>',
                encoding='utf-8',
            )
            band_file = (
                tmpdir
                / 'si-hse.bands'
            )

            def fake_run_qe_program(**kwargs):
                Path(
                    kwargs['output_file']
                ).write_text(
                    output_text,
                    encoding='utf-8',
                )
                band_file.write_text(
                    '&plot nbnd=8, nks=2 /\n',
                    encoding='utf-8',
                )

                return {
                    'returncode': 0,
                }

            with patch(
                'nanoworks.engine.qe.run_qe_program',
                side_effect=fake_run_qe_program,
            ) as run:
                workflow = run_bands_postprocess(
                    input_file=(
                        tmpdir
                        / 'bands.in'
                    ),
                    output_file=(
                        tmpdir
                        / 'bands.out'
                    ),
                    state_dir=state_dir,
                    band_file=band_file,
                    parallel_cores=4,
                )

            self.assertEqual(
                run.call_count,
                1,
            )
            self.assertEqual(
                run.call_args.kwargs[
                    'executable'
                ],
                'bands.x',
            )
            self.assertTrue(
                workflow['band_file'].is_file()
            )
            self.assertEqual(
                workflow['result']['program'],
                'BANDS',
            )

            input_text = (
                workflow['input_file']
                .read_text(
                    encoding='utf-8'
                )
            )

        self.assertIn(
            'lsym = .false.',
            input_text,
        )

    def test_run_hybrid_bands_composes_scf_and_bands(self):
        band_path = {
            'kpoints': [
                (0.0, 0.0, 0.0),
                (0.5, 0.0, 0.0),
            ],
            'npoints': 2,
        }
        scf_workflow = {
            'result': {
                'job_done': True,
            },
        }
        bands_workflow = {
            'result': {
                'job_done': True,
            },
        }
        band_data = {
            'nkpoints': 2,
            'nbands': 3,
            'eigenvalues_ev': [[
                [-5.0, -1.0, 1.0],
                [-4.0, 0.0, 2.0],
            ]],
        }

        with patch(
            'nanoworks.engine.qe.run_scf',
            return_value=scf_workflow,
        ) as run_scf_mock, patch(
            'nanoworks.engine.qe.run_bands_postprocess',
            return_value=bands_workflow,
        ) as run_bands_mock, patch(
            'nanoworks.engine.qe.parse_bands_x_output',
            return_value=band_data,
        ) as parse_bands_mock:
            workflow = run_hybrid_bands(
                atoms=Atoms('Si'),
                scf_input_file='scf.in',
                scf_output_file='scf.out',
                bands_input_file='bands.in',
                bands_output_file='bands.out',
                state_dir='state',
                band_file='bands.dat',
                pseudopotentials={
                    'Si': 'Si.upf',
                },
                pseudo_dir='/tmp/pseudos',
                cutoff_ev=400.0,
                band_path=band_path,
                qpoint_grid=(2, 1, 1),
                xc_calc='HSE06',
            )

        run_scf_mock.assert_called_once()
        parse_bands_mock.assert_called_once_with(
            'bands.dat',
            kpoint_indices=[0, 1],
        )
        run_bands_mock.assert_called_once_with(
            input_file='bands.in',
            output_file='bands.out',
            state_dir='state',
            band_file='bands.dat',
            parallel_cores=1,
            executable='bands.x',
            prefix='nanoworks',
            lsym=False,
        )
        self.assertEqual(
            run_scf_mock.call_args.kwargs[
                'exx_additional_kpoints'
            ]['qpoint_grid'],
            (2, 1, 1),
        )
        self.assertIs(
            workflow['scf'],
            scf_workflow,
        )
        self.assertIs(
            workflow['bands'],
            bands_workflow,
        )
        self.assertIs(
            workflow['band_data'],
            band_data,
        )
        self.assertEqual(
            workflow['band_path'],
            band_path,
        )

    def test_run_hybrid_bands_prepares_projected_band_data(self):
        band_path = {
            'kpoints': [
                (0.0, 0.0, 0.0),
                (0.5, 0.0, 0.0),
            ],
            'npoints': 2,
        }
        scf_workflow = {
            'result': {
                'job_done': True,
            },
        }
        bands_workflow = {
            'result': {
                'job_done': True,
            },
        }
        band_data = {
            'nkpoints': 2,
            'nbands': 3,
            'eigenvalues_ev': [[
                [-5.0, -1.0, 1.0],
                [-4.0, 0.0, 2.0],
            ]],
        }
        projection_workflow = {
            'projection_up_file': 'projection.projwfc_up',
            'projection_down_file': None,
        }
        raw_projection = {
            'nkpoints': 2,
            'nbands': 3,
            'states': [],
        }
        prepared_projection = {
            'nkpoints': 2,
            'nbands': 3,
            'projections': [],
        }

        with patch(
            'nanoworks.engine.qe.run_scf',
            return_value=scf_workflow,
        ), patch(
            'nanoworks.engine.qe.run_bands_postprocess',
            return_value=bands_workflow,
        ), patch(
            'nanoworks.engine.qe.parse_bands_x_output',
            return_value=band_data,
        ), patch(
            'nanoworks.engine.qe.run_band_projections',
            return_value=projection_workflow,
        ) as run_projection, patch(
            'nanoworks.engine.qe.parse_projwfc_band_file',
            return_value=raw_projection,
        ) as parse_projection, patch(
            'nanoworks.engine.qe.prepare_qe_band_projection_data',
            return_value=prepared_projection,
        ) as prepare_projection:
            workflow = run_hybrid_bands(
                atoms=Atoms('Si'),
                scf_input_file='scf.in',
                scf_output_file='scf.out',
                bands_input_file='bands.in',
                bands_output_file='bands.out',
                state_dir='state',
                band_file='bands.dat',
                pseudopotentials={
                    'Si': 'Si.upf',
                },
                pseudo_dir='/tmp/pseudos',
                cutoff_ev=400.0,
                band_path=band_path,
                qpoint_grid=(2, 1, 1),
                xc_calc='HSE06',
                projected_band=True,
                projections=[
                    {
                        'atoms': [0],
                        'orbital': 'p',
                        'color': 'red',
                        'label': 'Si p',
                    },
                ],
                projection_input_file='projection.in',
                projection_output_file='projection.out',
                projection_prefix='projection',
            )

        run_projection.assert_called_once_with(
            input_file='projection.in',
            output_file='projection.out',
            state_dir='state',
            projection_prefix='projection',
            spinpol=False,
            parallel_cores=1,
            executable='projwfc.x',
            prefix='nanoworks',
        )
        parse_projection.assert_called_once_with(
            'projection.projwfc_up',
            kpoint_indices=[0, 1],
        )
        prepare_projection.assert_called_once_with(
            raw_projection,
            projections=[
                {
                    'atoms': [0],
                    'orbital': 'p',
                    'color': 'red',
                    'label': 'Si p',
                },
            ],
        )
        self.assertIs(
            workflow['band_projections']['up'],
            prepared_projection,
        )
        self.assertIsNone(
            workflow['band_projections']['down'],
        )

    def test_run_hybrid_bands_rejects_nonhybrid_xc(self):
        with self.assertRaisesRegex(
            ValueError,
            'require a hybrid functional',
        ):
            run_hybrid_bands(
                atoms=Atoms('Si'),
                scf_input_file='scf.in',
                scf_output_file='scf.out',
                bands_input_file='bands.in',
                bands_output_file='bands.out',
                state_dir='state',
                band_file='bands.dat',
                pseudopotentials={
                    'Si': 'Si.upf',
                },
                pseudo_dir='/tmp/pseudos',
                cutoff_ev=400.0,
                band_path={
                    'kpoints': [
                        (0.0, 0.0, 0.0),
                    ],
                },
                qpoint_grid=(1, 1, 1),
                xc_calc='PBE',
            )

    def test_run_hybrid_dos_composes_scf_dos_and_pdos(self):
        scf_workflow = {
            'result': {
                'job_done': True,
            },
        }
        dos_workflow = {
            'dos_file': 'si-hse.dos',
        }
        pdos_workflow = {
            'pdos_tot_file': 'si-hse.pdos_tot',
        }

        with patch(
            'nanoworks.engine.qe.run_scf',
            return_value=scf_workflow,
        ) as run_scf_mock, patch(
            'nanoworks.engine.qe.run_dos',
            return_value=dos_workflow,
        ) as run_dos_mock, patch(
            'nanoworks.engine.qe.run_projwfc',
            return_value=pdos_workflow,
        ) as run_projwfc_mock:
            workflow = run_hybrid_dos(
                atoms=Atoms('Si'),
                scf_input_file='scf.in',
                scf_output_file='scf.out',
                dos_input_file='dos.in',
                dos_output_file='dos.out',
                dos_file='dos.dat',
                pdos_input_file='pdos.in',
                pdos_output_file='pdos.out',
                pdos_prefix='pdos',
                state_dir='state',
                pseudopotentials={
                    'Si': 'Si.upf',
                },
                pseudo_dir='/tmp/pseudos',
                cutoff_ev=500.0,
                kpoint_density=2.5,
                kpoint_size=(4, 4, 4),
                gamma=True,
                total_charge=0.0,
                nbands=24,
                spinpol=False,
                setup_params=None,
                xc_calc='HSE03',
                exx_fraction=0.28,
                omega=0.15,
                occupation='tetrahedra',
                emin=-8.0,
                emax=8.0,
                delta_e=0.1,
                bz_sum='tetrahedra',
                parallel_cores=4,
            )

        run_scf_mock.assert_called_once()
        scf_call = run_scf_mock.call_args.kwargs
        self.assertEqual(
            scf_call['kpoint_density'],
            2.5,
        )
        self.assertEqual(
            scf_call['kpoint_size'],
            (4, 4, 4),
        )
        self.assertEqual(
            scf_call['xc_calc'],
            'HSE03',
        )
        self.assertEqual(
            scf_call['exx_fraction'],
            0.28,
        )
        self.assertEqual(
            scf_call['omega'],
            0.15,
        )
        run_dos_mock.assert_called_once_with(
            input_file='dos.in',
            output_file='dos.out',
            state_dir='state',
            dos_file='dos.dat',
            emin=-8.0,
            emax=8.0,
            delta_e=0.1,
            bz_sum='tetrahedra',
            degauss=None,
            ngauss=None,
            parallel_cores=4,
            executable='dos.x',
            prefix='nanoworks',
        )
        run_projwfc_mock.assert_called_once_with(
            input_file='pdos.in',
            output_file='pdos.out',
            state_dir='state',
            pdos_prefix='pdos',
            emin=-8.0,
            emax=8.0,
            delta_e=0.1,
            degauss=None,
            ngauss=None,
            parallel_cores=4,
            executable='projwfc.x',
            prefix='nanoworks',
        )
        self.assertIs(
            workflow['scf'],
            scf_workflow,
        )
        self.assertIs(
            workflow['dos'],
            dos_workflow,
        )
        self.assertIs(
            workflow['pdos'],
            pdos_workflow,
        )

    def test_run_hybrid_dos_can_shift_windows_from_fermi_reference(self):
        scf_workflow = {
            'result': {
                'fermi_energy_ev': 5.0,
            },
        }
        dos_workflow = {
            'dos_file': 'si-hse.dos',
        }
        pdos_workflow = {
            'pdos_tot_file': 'si-hse.pdos_tot',
        }

        with patch(
            'nanoworks.engine.qe.run_scf',
            return_value=scf_workflow,
        ), patch(
            'nanoworks.engine.qe.run_dos',
            return_value=dos_workflow,
        ) as run_dos_mock, patch(
            'nanoworks.engine.qe.run_projwfc',
            return_value=pdos_workflow,
        ) as run_projwfc_mock:
            workflow = run_hybrid_dos(
                atoms=Atoms('Si'),
                scf_input_file='scf.in',
                scf_output_file='scf.out',
                dos_input_file='dos.in',
                dos_output_file='dos.out',
                dos_file='dos.dat',
                pdos_input_file='pdos.in',
                pdos_output_file='pdos.out',
                pdos_prefix='pdos',
                state_dir='state',
                pseudopotentials={
                    'Si': 'Si.upf',
                },
                pseudo_dir='/tmp/pseudos',
                cutoff_ev=500.0,
                xc_calc='HSE06',
                emin=-8.0,
                emax=8.0,
                delta_e=0.1,
                bz_sum='tetrahedra',
                relative_to_fermi=True,
            )

        self.assertEqual(
            run_dos_mock.call_args.kwargs['emin'],
            -3.0,
        )
        self.assertEqual(
            run_dos_mock.call_args.kwargs['emax'],
            13.0,
        )
        self.assertEqual(
            run_projwfc_mock.call_args.kwargs['emin'],
            -3.0,
        )
        self.assertEqual(
            run_projwfc_mock.call_args.kwargs['emax'],
            13.0,
        )
        self.assertEqual(
            workflow['fermi_energy_ev'],
            5.0,
        )

    def test_run_hybrid_dos_rejects_nonhybrid_xc(self):
        with self.assertRaisesRegex(
            ValueError,
            'requires a hybrid functional',
        ):
            run_hybrid_dos(
                atoms=Atoms('Si'),
                scf_input_file='scf.in',
                scf_output_file='scf.out',
                dos_input_file='dos.in',
                dos_output_file='dos.out',
                dos_file='dos.dat',
                pdos_input_file='pdos.in',
                pdos_output_file='pdos.out',
                pdos_prefix='pdos',
                state_dir='state',
                pseudopotentials={
                    'Si': 'Si.upf',
                },
                pseudo_dir='/tmp/pseudos',
                cutoff_ev=500.0,
                xc_calc='PBE',
            )

    def test_run_bands_parses_eigenvalues(self):
        atoms = bulk(
            'Si',
            'diamond',
            a=5.43,
        )

        band_path = build_band_path(
            atoms=atoms,
            path='GX',
            npoints=2,
        )

        output_text = """
         Program PWSCF v.7.2 starts

         End of band structure calculation

              k = 0.0000 0.0000 0.0000 ( 123 PWs)   bands (ev):

            -5.0000  -1.0000   1.0000

              k = 0.5000 0.0000 0.5000 ( 120 PWs)   bands (ev):

            -4.5000  -0.5000   1.5000

         JOB DONE.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            state_dir = (
                tmpdir
                / 'state'
            )

            save_dir = (
                state_dir
                / 'nanoworks.save'
            )

            save_dir.mkdir(
                parents=True
            )

            (
                save_dir
                / 'data-file-schema.xml'
            ).write_text(
                '<qes/>',
                encoding='utf-8',
            )

            def fake_run_qe_program(**kwargs):
                Path(
                    kwargs['output_file']
                ).write_text(
                    output_text,
                    encoding='utf-8',
                )

                return {
                    'returncode': 0,
                }

            with patch(
                'nanoworks.engine.qe.run_qe_program',
                side_effect=fake_run_qe_program,
            ) as run:
                workflow = run_bands(
                    atoms=atoms,
                    input_file=tmpdir / 'bands.in',
                    output_file=tmpdir / 'bands.out',
                    state_dir=state_dir,
                    pseudopotentials={
                        'Si': 'Si.upf',
                    },
                    pseudo_dir='/tmp/pseudos',
                    cutoff_ev=400.0,
                    band_path=band_path,
                    nbands=3,
                )

            self.assertEqual(
                run.call_count,
                1,
            )

            self.assertTrue(
                workflow['input_file'].is_file()
            )

            input_text = (
                workflow['input_file']
                .read_text(
                    encoding='utf-8'
                )
            )

        self.assertIn(
            "calculation = 'bands'",
            input_text,
        )

        self.assertIn(
            "K_POINTS crystal",
            input_text,
        )

        self.assertEqual(
            workflow['bands']['nspins'],
            1,
        )

        self.assertEqual(
            workflow['bands']['nkpoints'],
            2,
        )

        self.assertEqual(
            workflow['bands']['nbands'],
            3,
        )

        self.assertEqual(
            workflow['bands']['eigenvalues_ev'],
            [[
                [-5.0, -1.0, 1.0],
                [-4.5, -0.5, 1.5],
            ]],
        )

    def test_parse_dos_output(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            dos_file = (
                Path(tmpdir)
                / 'nanoworks.dos'
            )

            dos_file.write_text(
                "# E (eV) DOS(E) Int DOS(E)\n"
                "-1.0000  0.1000  0.0100\n"
                " 0.0000  0.5000  0.2000\n"
                " 1.0000  0.2500  0.6000\n",
                encoding='utf-8',
            )

            result = parse_dos_output(
                dos_file
            )

        self.assertEqual(
            result['npoints'],
            3,
        )

        self.assertEqual(
            result['energies_ev'],
            [-1.0, 0.0, 1.0],
        )

        self.assertEqual(
            result['dos'],
            [0.1, 0.5, 0.25],
        )

        self.assertEqual(
            result['integrated_dos'],
            [0.01, 0.2, 0.6],
        )

        self.assertFalse(
            result['spin_polarized']
        )

        self.assertIsNone(
            result['dos_up']
        )

        self.assertIsNone(
            result['dos_down']
        )

    def test_parse_spin_polarized_dos_output(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            dos_file = (
                Path(tmpdir)
                / 'nanoworks.dos'
            )

            dos_file.write_text(
                "# E (eV) dosup(E) dosdw(E) Int dos(E)\n"
                "-1.0000  0.1000  0.0400  0.0100\n"
                " 0.0000  0.5000  0.2000  0.2500\n"
                " 1.0000  0.2500  0.1500  0.7000\n",
                encoding='utf-8',
            )

            result = parse_dos_output(
                dos_file
            )

        self.assertTrue(
            result['spin_polarized']
        )

        self.assertEqual(
            result['energies_ev'],
            [-1.0, 0.0, 1.0],
        )

        self.assertEqual(
            result['dos_up'],
            [0.1, 0.5, 0.25],
        )

        self.assertEqual(
            result['dos_down'],
            [0.04, 0.2, 0.15],
        )

        self.assertEqual(
            result['dos'],
            [0.14, 0.7, 0.4],
        )

        self.assertEqual(
            result['integrated_dos'],
            [0.01, 0.25, 0.7],
        )

    def test_parse_dos_output_rejects_mixed_spin_columns(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            dos_file = (
                Path(tmpdir)
                / 'nanoworks.dos'
            )

            dos_file.write_text(
                "# inconsistent DOS data\n"
                "-1.0 0.1 0.01\n"
                " 0.0 0.5 0.2 0.25\n",
                encoding='utf-8',
            )

            with self.assertRaisesRegex(
                ValueError,
                'inconsistent spin column counts',
            ):
                parse_dos_output(
                    dos_file
                )

    def test_parse_dos_output_supports_fortran_exponents(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            dos_file = (
                Path(tmpdir)
                / 'nanoworks.dos'
            )

            dos_file.write_text(
                "# E DOS IntDOS\n"
                "1.000D+00 2.500D-01 6.000D-01\n",
                encoding='utf-8',
            )

            result = parse_dos_output(
                dos_file
            )

        self.assertEqual(
            result['energies_ev'],
            [1.0],
        )

        self.assertEqual(
            result['dos'],
            [0.25],
        )

    def test_parse_dos_output_rejects_missing_data(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            dos_file = (
                Path(tmpdir)
                / 'nanoworks.dos'
            )

            dos_file.write_text(
                "# no DOS data\n",
                encoding='utf-8',
            )

            with self.assertRaisesRegex(
                ValueError,
                'No DOS data could be parsed',
            ):
                parse_dos_output(
                    dos_file
                )

def test_render_projwfc_input(self):
    text = render_projwfc_input(
        prefix='nanoworks',
        outdir='/tmp/qe-state',
        filpdos='/tmp/gaas-pdos',
        filproj='/tmp/gaas-projections.dat',
        emin=1.0,
        emax=10.0,
        delta_e=0.02,
    )

    self.assertIn(
        '&PROJWFC',
        text,
    )

    self.assertIn(
        "prefix = 'nanoworks'",
        text,
    )

    self.assertIn(
        "outdir = '/tmp/qe-state'",
        text,
    )

    self.assertIn(
        "filpdos = '/tmp/gaas-pdos'",
        text,
    )

    self.assertIn(
        "filproj = '/tmp/gaas-projections.dat'",
        text,
    )

    self.assertIn(
        'Emin = 1',
        text,
    )

    self.assertIn(
        'Emax = 10',
        text,
    )

    self.assertIn(
        'DeltaE = 0.02',
        text,
    )

    self.assertIn(
        'lsym = .true.',
        text,
    )

    self.assertNotIn(
        'degauss',
        text,
    )

def test_render_projwfc_input_omits_optional_projection_file(self):
    text = render_projwfc_input(
        prefix='nanoworks',
        outdir='/tmp/qe-state',
        filpdos='/tmp/gaas-pdos',
    )

    self.assertNotIn(
        'filproj',
        text,
    )

def test_render_projwfc_input_rejects_invalid_energy_range(self):
    with self.assertRaisesRegex(
        ValueError,
        'Emax must be greater than Emin',
    ):
        render_projwfc_input(
            emin=5.0,
            emax=-5.0,
        )

def test_run_projwfc_requires_qe_state(self):
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(
            tmpdir
        )

        with self.assertRaisesRegex(
            FileNotFoundError,
            'valid QE electronic state',
        ):
            run_projwfc(
                input_file=tmpdir / 'pdos.in',
                output_file=tmpdir / 'pdos.out',
                state_dir=tmpdir / 'state',
                pdos_prefix=tmpdir / 'pdos',
            )

def test_parse_projwfc_pdos_p_file(self):
    with tempfile.TemporaryDirectory() as tmpdir:
        pdos_file = (
            Path(tmpdir)
            / 'test.pdos_atm#2(As)_wfc#3(p)'
        )

        pdos_file.write_text(
            "# E (eV) ldos(E) pdos(E) pdos(E) pdos(E)\n"
            "2.048 0.186E-01 0.621E-02 0.622E-02 0.623E-02\n"
            "2.078 0.175E-01 0.583E-02 0.584E-02 0.585E-02\n",
            encoding='utf-8',
        )

        result = parse_projwfc_pdos_file(
            pdos_file
        )

    self.assertEqual(
        result['atom_index'],
        2,
    )

    self.assertEqual(
        result['symbol'],
        'As',
    )

    self.assertEqual(
        result['wfc_index'],
        3,
    )

    self.assertEqual(
        result['orbital'],
        'p',
    )

    self.assertEqual(
        result['components']['pz'],
        [0.00621, 0.00583],
    )

    self.assertEqual(
        result['components']['px'],
        [0.00622, 0.00584],
    )

    self.assertEqual(
        result['components']['py'],
        [0.00623, 0.00585],
    )
    
    self.assertFalse(
        result['spin_polarized']
    )

    self.assertIsNone(
        result['ldos_up']
    )

    self.assertIsNone(
        result['components_up']
    )

def test_parse_spin_polarized_projwfc_pdos_p_file(self):
    with tempfile.TemporaryDirectory() as tmpdir:
        pdos_file = (
            Path(tmpdir)
            / 'test.pdos_atm#2(As)_wfc#3(p)'
        )

        pdos_file.write_text(
            "# E ldosup ldosdw "
            "pzup pzdw pxup pxdw pyup pydw\n"
            "1.0 0.60 0.30 "
            "0.10 0.05 0.20 0.10 0.30 0.15\n"
            "2.0 0.90 0.60 "
            "0.20 0.10 0.30 0.20 0.40 0.30\n",
            encoding='utf-8',
        )

        result = parse_projwfc_pdos_file(
            pdos_file
        )

    self.assertTrue(
        result['spin_polarized']
    )

    self.assertEqual(
        result['ldos_up'],
        [0.6, 0.9],
    )

    self.assertEqual(
        result['ldos_down'],
        [0.3, 0.6],
    )

    self.assertEqual(
        result['components_up']['pz'],
        [0.1, 0.2],
    )

    self.assertEqual(
        result['components_down']['pz'],
        [0.05, 0.1],
    )

    self.assertEqual(
        result['components_up']['px'],
        [0.2, 0.3],
    )

    self.assertEqual(
        result['components_down']['py'],
        [0.15, 0.3],
    )

    self.assertAlmostEqual(
        result['components']['pz'][0],
        0.15,
    )

    self.assertAlmostEqual(
        result['ldos'][1],
        1.5,
    )

def test_parse_projwfc_pdos_d_file(self):
    with tempfile.TemporaryDirectory() as tmpdir:
        pdos_file = (
            Path(tmpdir)
            / 'test.pdos_atm#1(Ga)_wfc#1(d)'
        )

        pdos_file.write_text(
            "# E ldos d1 d2 d3 d4 d5\n"
            "2.048 0.198E-02 "
            "0.654E-03 0.223E-03 0.224E-03 "
            "0.438E-03 0.439E-03\n",
            encoding='utf-8',
        )

        result = parse_projwfc_pdos_file(
            pdos_file
        )

    self.assertEqual(
        result['orbital'],
        'd',
    )

    self.assertEqual(
        result['components']['d3z2_r2'],
        [0.000654],
    )

    self.assertEqual(
        result['components']['dxz'],
        [0.000223],
    )

    self.assertEqual(
        result['components']['dyz'],
        [0.000224],
    )

    self.assertEqual(
        result['components']['dx2_y2'],
        [0.000438],
    )

    self.assertEqual(
        result['components']['dxy'],
        [0.000439],
    )

def test_parse_projwfc_pdos_s_file(self):
    with tempfile.TemporaryDirectory() as tmpdir:
        pdos_file = (
            Path(tmpdir)
            / 'test.pdos_atm#1(Ga)_wfc#2(s)'
        )

        pdos_file.write_text(
            "# E ldos pdos\n"
            "2.048 0.842E-01 0.842E-01\n",
            encoding='utf-8',
        )

        result = parse_projwfc_pdos_file(
            pdos_file
        )

    self.assertEqual(
        result['components']['s'],
        [0.0842],
    )

def test_aggregate_projwfc_pdos(self):
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(
            tmpdir
        )

        prefix = (
            tmpdir
            / 'nanoworks-pdos'
        )

        ga_s = Path(
            str(prefix)
            + '.pdos_atm#1(Ga)_wfc#2(s)'
        )

        ga_p = Path(
            str(prefix)
            + '.pdos_atm#1(Ga)_wfc#3(p)'
        )

        as_p = Path(
            str(prefix)
            + '.pdos_atm#2(As)_wfc#3(p)'
        )

        ga_s.write_text(
            "# E ldos s\n"
            "1.0 0.10 0.10\n"
            "2.0 0.20 0.20\n",
            encoding='utf-8',
        )

        ga_p.write_text(
            "# E ldos pz px py\n"
            "1.0 0.60 0.10 0.20 0.30\n"
            "2.0 0.90 0.20 0.30 0.40\n",
            encoding='utf-8',
        )

        as_p.write_text(
            "# E ldos pz px py\n"
            "1.0 0.30 0.05 0.10 0.15\n"
            "2.0 0.60 0.10 0.20 0.30\n",
            encoding='utf-8',
        )

        result = aggregate_projwfc_pdos(
            prefix
        )

    self.assertEqual(
        result['energies_ev'],
        [1.0, 2.0],
    )

    self.assertEqual(
        result['s_total'],
        [0.1, 0.2],
    )

    self.assertEqual(
        result['p_total'],
        [0.9, 1.5],
    )

    self.assertAlmostEqual(
        result['total'][0],
        1.0,
    )

    self.assertAlmostEqual(
        result['total'][1],
        1.7,
    )

    self.assertEqual(
        result['pz'],
        [0.15, 0.30],
    )

    self.assertEqual(
        result['px'],
        [0.30, 0.50],
    )

    self.assertEqual(
        result['py'],
        [0.45, 0.70],
    )

    self.assertEqual(
        result['d_total'],
        [0.0, 0.0],
    )

    self.assertEqual(
        result['f_total'],
        [0.0, 0.0],
    )

    self.assertFalse(
        result['spin_polarized']
    )

    self.assertIsNone(
        result['spin_up']
    )

    self.assertIsNone(
        result['spin_down']
    )

def test_aggregate_spin_polarized_projwfc_pdos(self):
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(
            tmpdir
        )

        prefix = (
            tmpdir
            / 'nanoworks-pdos'
        )

        s_file = Path(
            str(prefix)
            + '.pdos_atm#1(Ga)_wfc#1(s)'
        )

        p_file = Path(
            str(prefix)
            + '.pdos_atm#1(Ga)_wfc#2(p)'
        )

        s_file.write_text(
            "# E ldosup ldosdw sup sdw\n"
            "1.0 0.10 0.04 0.10 0.04\n"
            "2.0 0.20 0.08 0.20 0.08\n",
            encoding='utf-8',
        )

        p_file.write_text(
            "# E ldosup ldosdw "
            "pzup pzdw pxup pxdw pyup pydw\n"
            "1.0 0.60 0.30 "
            "0.10 0.05 0.20 0.10 0.30 0.15\n"
            "2.0 0.90 0.60 "
            "0.20 0.10 0.30 0.20 0.40 0.30\n",
            encoding='utf-8',
        )

        result = aggregate_projwfc_pdos(
            prefix
        )

    self.assertTrue(
        result['spin_polarized']
    )

    self.assertEqual(
        result['spin_up']['s_total'],
        [0.1, 0.2],
    )

    self.assertEqual(
        result['spin_down']['s_total'],
        [0.04, 0.08],
    )

    self.assertEqual(
        result['spin_up']['p_total'],
        [0.6, 0.9],
    )

    self.assertEqual(
        result['spin_down']['p_total'],
        [0.3, 0.6],
    )

    self.assertEqual(
        result['spin_up']['pz'],
        [0.1, 0.2],
    )

    self.assertEqual(
        result['spin_down']['py'],
        [0.15, 0.3],
    )

    self.assertAlmostEqual(
        result['spin_up']['total'][0],
        0.7,
    )

    self.assertAlmostEqual(
        result['spin_down']['total'][0],
        0.34,
    )

    self.assertAlmostEqual(
        result['total'][0],
        1.04,
    )

def test_aggregate_projwfc_pdos_rejects_mismatched_energy_grid(self):
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(
            tmpdir
        )

        prefix = (
            tmpdir
            / 'nanoworks-pdos'
        )

        first = Path(
            str(prefix)
            + '.pdos_atm#1(Ga)_wfc#2(s)'
        )

        second = Path(
            str(prefix)
            + '.pdos_atm#2(As)_wfc#2(s)'
        )

        first.write_text(
            "# E ldos s\n"
            "1.0 0.10 0.10\n"
            "2.0 0.20 0.20\n",
            encoding='utf-8',
        )

        second.write_text(
            "# E ldos s\n"
            "1.0 0.10 0.10\n"
            "2.1 0.20 0.20\n",
            encoding='utf-8',
        )

        with self.assertRaisesRegex(
            ValueError,
            'same energy grid',
        ):
            aggregate_projwfc_pdos(
                prefix
            )

    def test_parse_pw_relaxed_atomic_positions(self):
        atoms = bulk(
            'Si',
            'diamond',
            a=5.43,
        )

        output_text = """
        Begin final coordinates
        ATOMIC_POSITIONS (angstrom)
        Si 0.100000 0.200000 0.300000
        Si 1.400000 1.500000 1.600000
        End final coordinates
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = (
                Path(tmpdir)
                / 'relax.out'
            )

            output_file.write_text(
                output_text,
                encoding='utf-8',
            )

            relaxed = parse_pw_relaxed_structure(
                output_file,
                atoms,
            )

        self.assertAlmostEqual(
            relaxed.positions[0, 0],
            0.1,
        )
        self.assertAlmostEqual(
            relaxed.positions[0, 1],
            0.2,
        )
        self.assertAlmostEqual(
            relaxed.positions[1, 2],
            1.6,
        )
        self.assertTrue(
            relaxed.cell == atoms.cell
        )

    def test_parse_pw_relaxed_cell_and_crystal_positions(self):
        atoms = bulk(
            'Si',
            'diamond',
            a=5.43,
        )

        output_text = """
        Begin final coordinates
        CELL_PARAMETERS (angstrom)
        3.000000 0.000000 0.000000
        0.000000 3.000000 0.000000
        0.000000 0.000000 3.000000
        ATOMIC_POSITIONS (crystal)
        Si 0.000000 0.000000 0.000000
        Si 0.250000 0.250000 0.250000
        End final coordinates
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = (
                Path(tmpdir)
                / 'vc-relax.out'
            )

            output_file.write_text(
                output_text,
                encoding='utf-8',
            )

            relaxed = parse_pw_relaxed_structure(
                output_file,
                atoms,
            )

        self.assertAlmostEqual(
            relaxed.cell[0, 0],
            3.0,
        )
        self.assertAlmostEqual(
            relaxed.cell[1, 1],
            3.0,
        )
        self.assertAlmostEqual(
            relaxed.positions[1, 0],
            0.75,
        )
        self.assertAlmostEqual(
            relaxed.positions[1, 1],
            0.75,
        )
        self.assertAlmostEqual(
            relaxed.positions[1, 2],
            0.75,
        )

    def test_parse_pw_relaxed_structure_requires_final_coordinates(self):
        atoms = bulk(
            'Si',
            'diamond',
            a=5.43,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = (
                Path(tmpdir)
                / 'relax.out'
            )

            output_file.write_text(
                "Program PWSCF v.7.2\nJOB DONE.\n",
                encoding='utf-8',
            )

            with self.assertRaisesRegex(
                ValueError,
                'No final relaxed coordinates',
            ):
                parse_pw_relaxed_structure(
                    output_file,
                    atoms,
                )

    def test_run_relax_returns_optimized_structure(self):
        atoms = bulk(
            'Si',
            'diamond',
            a=5.43,
        )

        output_text = """
        Program PWSCF v.7.2 starts

        !    total energy              =    -15.00000000 Ry

        bfgs converged in 5 scf cycles and 4 bfgs steps

        Begin final coordinates
        ATOMIC_POSITIONS (angstrom)
        Si 0.100000 0.200000 0.300000
        Si 1.400000 1.500000 1.600000
        End final coordinates

        JOB DONE.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(
                tmpdir
            )

            def fake_run_qe_program(**kwargs):
                Path(
                    kwargs['output_file']
                ).write_text(
                    output_text,
                    encoding='utf-8',
                )

                return {
                    'returncode': 0,
                }

            with patch(
                'nanoworks.engine.qe.run_qe_program',
                side_effect=fake_run_qe_program,
            ) as run:
                workflow = run_relax(
                    atoms=atoms,
                    input_file=tmpdir / 'relax.in',
                    output_file=tmpdir / 'relax.out',
                    state_dir=tmpdir / 'state',
                    pseudopotentials={
                        'Si': 'Si.upf',
                    },
                    pseudo_dir='/tmp/pseudos',
                    cutoff_ev=400.0,
                    optimizer='LBFGS',
                    max_force=0.05,
                    max_step=0.20,
                    relax_cell=[
                        False,
                        False,
                        False,
                        False,
                        False,
                        False,
                    ],
                )

            self.assertEqual(
                run.call_count,
                1,
            )

            self.assertTrue(
                workflow['input_file'].is_file()
            )

            input_text = (
                workflow['input_file']
                .read_text(
                    encoding='utf-8',
                )
            )

        self.assertIn(
            "calculation = 'relax'",
            input_text,
        )
        self.assertEqual(
            workflow['calculation'],
            'relax',
        )
        self.assertTrue(
            workflow['geometry_converged']
        )
        self.assertAlmostEqual(
            workflow['result']['total_energy_ry'],
            -15.0,
        )
        self.assertAlmostEqual(
            workflow['atoms'].positions[0, 0],
            0.1,
        )
        self.assertAlmostEqual(
            workflow['atoms'].positions[1, 2],
            1.6,
        )

    def test_run_relax_rejects_unconverged_geometry(self):
        atoms = bulk(
            'Si',
            'diamond',
            a=5.43,
        )

        output_text = """
        Program PWSCF v.7.2 starts

        !    total energy              =    -15.00000000 Ry

        The maximum number of steps has been reached.

        Begin final coordinates
        ATOMIC_POSITIONS (angstrom)
        Si 0.100000 0.200000 0.300000
        Si 1.400000 1.500000 1.600000
        End final coordinates

        JOB DONE.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(
                tmpdir
            )

            def fake_run_qe_program(**kwargs):
                Path(
                    kwargs['output_file']
                ).write_text(
                    output_text,
                    encoding='utf-8',
                )

                return {
                    'returncode': 0,
                }

            with patch(
                'nanoworks.engine.qe.run_qe_program',
                side_effect=fake_run_qe_program,
            ):
                with self.assertRaisesRegex(
                    RuntimeError,
                    'did not report BFGS convergence',
                ):
                    run_relax(
                        atoms=atoms,
                        input_file=tmpdir / 'relax.in',
                        output_file=tmpdir / 'relax.out',
                        state_dir=tmpdir / 'state',
                        pseudopotentials={
                            'Si': 'Si.upf',
                        },
                        pseudo_dir='/tmp/pseudos',
                        cutoff_ev=400.0,
                        optimizer='LBFGS',
                        max_force=0.05,
                        max_step=0.20,
                        relax_cell=[
                            False,
                            False,
                            False,
                            False,
                            False,
                            False,
                        ],
                    )

    def test_qe_magnetic_species_support_uniform_moment(self):
        atoms = Atoms(
            'Fe2',
            positions=[
                [0.0, 0.0, 0.0],
                [1.0, 1.0, 1.0],
            ],
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(
                tmpdir
            )

            (
                tmpdir
                / 'Fe.upf'
            ).write_text(
                '<PP_HEADER z_valence="8.0" />\n',
                encoding='utf-8',
            )

            model = build_qe_magnetic_species(
                atoms=atoms,
                pseudopotentials={
                    'Fe': 'Fe.upf',
                },
                pseudo_dir=tmpdir,
                magnetic_moments=[
                    2.0,
                    2.0,
                ],
            )

        self.assertEqual(
            model['ntyp'],
            1,
        )

        self.assertEqual(
            model['species_labels'],
            [
                'Fe',
                'Fe',
            ],
        )

        self.assertAlmostEqual(
            model[
                'starting_magnetizations'
            ][
                'starting_magnetization(1)'
            ],
            0.25,
        )

    def test_qe_magnetic_species_split_antiferromagnetic_atoms(self):
        atoms = Atoms(
            'Fe2O3',
            positions=[
                [0.0, 0.0, 0.0],
                [1.0, 1.0, 1.0],
                [2.0, 0.0, 0.0],
                [0.0, 2.0, 0.0],
                [0.0, 0.0, 2.0],
            ],
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(
                tmpdir
            )

            (
                tmpdir
                / 'Fe.upf'
            ).write_text(
                '<PP_HEADER z_valence="8.0" />\n',
                encoding='utf-8',
            )

            (
                tmpdir
                / 'O.upf'
            ).write_text(
                '<PP_HEADER z_valence="6.0" />\n',
                encoding='utf-8',
            )

            model = build_qe_magnetic_species(
                atoms=atoms,
                pseudopotentials={
                    'Fe': 'Fe.upf',
                    'O': 'O.upf',
                },
                pseudo_dir=tmpdir,
                magnetic_moments=[
                    4.0,
                    -4.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            )

        self.assertEqual(
            model['ntyp'],
            3,
        )

        self.assertEqual(
            model['species_labels'],
            [
                'Fe1',
                'Fe2',
                'O',
                'O',
                'O',
            ],
        )

        self.assertEqual(
            [
                species[0]
                for species in model['species']
            ],
            [
                'Fe1',
                'Fe2',
                'O',
            ],
        )

        self.assertAlmostEqual(
            model[
                'starting_magnetizations'
            ][
                'starting_magnetization(1)'
            ],
            0.5,
        )

        self.assertAlmostEqual(
            model[
                'starting_magnetizations'
            ][
                'starting_magnetization(2)'
            ],
            -0.5,
        )

        self.assertEqual(
            model[
                'starting_magnetizations'
            ][
                'starting_magnetization(3)'
            ],
            0.0,
        )

        self.assertEqual(
            model[
                'positions'
            ][
                'positions'
            ][0][0],
            'Fe1',
        )

        self.assertEqual(
            model[
                'positions'
            ][
                'positions'
            ][1][0],
            'Fe2',
        )

    def test_qe_magnetic_species_rejects_zero_moments(self):
        atoms = Atoms(
            'Fe',
            positions=[
                [0.0, 0.0, 0.0],
            ],
        )

        with self.assertRaisesRegex(
            ValueError,
            'at least one non-zero',
        ):
            build_qe_magnetic_species(
                atoms=atoms,
                pseudopotentials={
                    'Fe': 'Fe.upf',
                },
                pseudo_dir='/tmp',
                magnetic_moments=[
                    0.0,
                ],
            )

    def test_render_spin_polarized_scf_input(self):
        atoms = Atoms(
            'Fe2O3',
            positions=[
                [0.0, 0.0, 0.0],
                [1.0, 1.0, 1.0],
                [2.0, 0.0, 0.0],
                [0.0, 2.0, 0.0],
                [0.0, 0.0, 2.0],
            ],
            cell=[
                [4.0, 0.0, 0.0],
                [0.0, 4.0, 0.0],
                [0.0, 0.0, 4.0],
            ],
            pbc=True,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(
                tmpdir
            )

            (
                tmpdir
                / 'Fe.upf'
            ).write_text(
                '<PP_HEADER z_valence="8.0" />\n',
                encoding='utf-8',
            )

            (
                tmpdir
                / 'O.upf'
            ).write_text(
                '<PP_HEADER z_valence="6.0" />\n',
                encoding='utf-8',
            )

            text = render_scf_input(
                atoms=atoms,
                pseudopotentials={
                    'Fe': 'Fe.upf',
                    'O': 'O.upf',
                },
                pseudo_dir=tmpdir,
                cutoff_ev=500.0,
                kpoint_size=(
                    4,
                    4,
                    4,
                ),
                spinpol=True,
                magnetic_moments=[
                    4.0,
                    -4.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            )

        self.assertIn(
            'ntyp = 3',
            text,
        )

        self.assertIn(
            'nspin = 2',
            text,
        )

        self.assertIn(
            'starting_magnetization(1) = 0.5',
            text,
        )

        self.assertIn(
            'starting_magnetization(2) = -0.5',
            text,
        )

        self.assertIn(
            'starting_magnetization(3) = 0',
            text,
        )

        self.assertIn(
            'Fe1 55.84500000 Fe.upf',
            text,
        )

        self.assertIn(
            'Fe2 55.84500000 Fe.upf',
            text,
        )

        self.assertIn(
            'O 15.99900000 O.upf',
            text,
        )

        self.assertIn(
            'Fe1 0.000000000000 0.000000000000 0.000000000000',
            text,
        )

        self.assertIn(
            'Fe2 1.000000000000 1.000000000000 1.000000000000',
            text,
        )

    def test_render_spin_input_requires_magnetic_moments(self):
        atoms = bulk(
            'Fe',
            'bcc',
            a=2.87,
        )

        with self.assertRaisesRegex(
            ValueError,
            'requires initial magnetic moments',
        ):
            render_scf_input(
                atoms=atoms,
                pseudopotentials={
                    'Fe': 'Fe.upf',
                },
                pseudo_dir='/tmp',
                cutoff_ev=500.0,
                kpoint_size=(
                    4,
                    4,
                    4,
                ),
                spinpol=True,
            )

    def test_render_nonspin_input_rejects_magnetic_moments(self):
        atoms = bulk(
            'Fe',
            'bcc',
            a=2.87,
        )

        with self.assertRaisesRegex(
            ValueError,
            'require spinpol=True',
        ):
            render_scf_input(
                atoms=atoms,
                pseudopotentials={
                    'Fe': 'Fe.upf',
                },
                pseudo_dir='/tmp',
                cutoff_ev=500.0,
                kpoint_size=(
                    4,
                    4,
                    4,
                ),
                spinpol=False,
                magnetic_moments=[
                    2.0,
                ],
            )

    def test_parse_spin_polarized_pw_bands_output(self):
        output_text = """
        ------ SPIN UP ------------

        k = 0.0000 0.0000 0.0000 (100 PWs) bands (ev):
        -5.0 -1.0 1.0

        k = 0.5000 0.0000 0.5000 (90 PWs) bands (ev):
        -4.5 -0.5 1.5

        ------ SPIN DOWN ----------

        k = 0.0000 0.0000 0.0000 (100 PWs) bands (ev):
        -4.8 -0.8 1.2

        k = 0.5000 0.0000 0.5000 (90 PWs) bands (ev):
        -4.3 -0.3 1.7

        JOB DONE.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = (
                Path(tmpdir)
                / 'spin-bands.out'
            )

            output_file.write_text(
                output_text,
                encoding='utf-8',
            )

            result = parse_pw_bands_output(
                output_file
            )

        self.assertTrue(
            result['spin_polarized']
        )
        self.assertEqual(
            result['nspins'],
            2,
        )
        self.assertEqual(
            result['nkpoints'],
            2,
        )
        self.assertEqual(
            result['nbands'],
            3,
        )
        self.assertEqual(
            result['eigenvalues_ev'][0][1],
            [-4.5, -0.5, 1.5],
        )
        self.assertEqual(
            result['eigenvalues_ev'][1][1],
            [-4.3, -0.3, 1.7],
        )

    def test_prepare_spin_polarized_qe_band_data(self):
        bands = {
            'spin_polarized': True,
            'nspins': 2,
            'nkpoints': 2,
            'nbands': 2,
            'eigenvalues_ev': [
                [
                    [0.0, 2.0],
                    [0.5, 2.5],
                ],
                [
                    [0.2, 2.2],
                    [0.7, 2.7],
                ],
            ],
        }

        band_path = {
            'distances': [
                0.0,
                0.5,
            ],
            'special_distances': [
                0.0,
                0.5,
            ],
            'labels': [
                'G',
                'X',
            ],
        }

        result = prepare_qe_band_data(
            bands=bands,
            band_path=band_path,
            reference_energy=1.0,
        )

        self.assertTrue(
            result['spin_polarized']
        )
        self.assertEqual(
            result['nspins'],
            2,
        )
        self.assertIsNone(
            result['eigenvalues_ev']
        )
        self.assertEqual(
            result['eigenvalues_up_ev'],
            [
                [-1.0, 1.0],
                [-0.5, 1.5],
            ],
        )
        self.assertEqual(
            result['eigenvalues_down_ev'],
            [
                [-0.8, 1.2],
                [-0.3, 1.7],
            ],
        )

def test_run_band_projections_requires_qe_state(self):
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(
            tmpdir
        )

        with self.assertRaisesRegex(
            FileNotFoundError,
            'valid QE bands state',
        ):
            run_band_projections(
                input_file=tmpdir / 'proj.in',
                output_file=tmpdir / 'proj.out',
                state_dir=tmpdir / 'state',
                projection_prefix=(
                    tmpdir
                    / 'bands-proj'
                ),
            )


def test_run_spin_polarized_band_projections(self):
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(
            tmpdir
        )

        state_dir = (
            tmpdir
            / 'state'
        )

        save_dir = (
            state_dir
            / 'nanoworks.save'
        )

        save_dir.mkdir(
            parents=True
        )

        (
            save_dir
            / 'data-file-schema.xml'
        ).write_text(
            '<espresso/>',
            encoding='utf-8',
        )

        projection_prefix = (
            tmpdir
            / 'bands-proj'
        )

        def fake_run_qe_program(**kwargs):
            Path(
                kwargs['output_file']
            ).write_text(
                'JOB DONE.\n',
                encoding='utf-8',
            )

            Path(
                str(projection_prefix)
                + '.projwfc_up'
            ).write_text(
                'up projections\n',
                encoding='utf-8',
            )

            Path(
                str(projection_prefix)
                + '.projwfc_down'
            ).write_text(
                'down projections\n',
                encoding='utf-8',
            )

            return {
                'returncode': 0,
            }

        with patch(
            'nanoworks.engine.qe.run_qe_program',
            side_effect=fake_run_qe_program,
        ):
            workflow = run_band_projections(
                input_file=tmpdir / 'proj.in',
                output_file=tmpdir / 'proj.out',
                state_dir=state_dir,
                projection_prefix=projection_prefix,
                spinpol=True,
            )

        input_text = (
            workflow['input_file']
            .read_text(
                encoding='utf-8',
            )
        )

    self.assertIn(
        "filproj = "
        f"'{projection_prefix}'",
        input_text,
    )

    self.assertIn(
        'lsym = .false.',
        input_text,
    )

    self.assertTrue(
        workflow[
            'projection_up_file'
        ].is_file()
    )

    self.assertTrue(
        workflow[
            'projection_down_file'
        ].is_file()
    )

    self.assertEqual(
        len(
            workflow[
                'projection_files'
            ]
        ),
        2,
    )

    def test_parse_projwfc_band_file(self):
        projection_text = """
        8 8 8 8 8 8 1 1
        0 5.0 0.0 0.0 0.0 0.0 0.0
        1.0 0.0 0.0
        0.0 1.0 0.0
        0.0 0.0 1.0
        10.0 4.0 30.0 9
        1 Fe 16.0
        1 0.0 0.0 0.0 1
        2 2 2
        F F
        1 1 Fe 4S 1 0 1
        1 1 0.10
        1 2 0.20
        2 1 0.30
        2 2 0.40
        2 1 Fe 3D 2 2 1
        1 1 0.50
        1 2 0.60
        2 1 0.70
        2 2 0.80
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            projection_file = (
                Path(tmpdir)
                / 'bands.projwfc_up'
            )

            projection_file.write_text(
                projection_text,
                encoding='utf-8',
            )

            result = parse_projwfc_band_file(
                projection_file
            )

        self.assertEqual(
            result['natoms'],
            1,
        )
        self.assertEqual(
            result['natomic_states'],
            2,
        )
        self.assertEqual(
            result['nkpoints'],
            2,
        )
        self.assertEqual(
            result['nbands'],
            2,
        )
        self.assertEqual(
            result['states'][0]['atom_index'],
            0,
        )
        self.assertEqual(
            result['states'][0]['orbital'],
            's',
        )
        self.assertEqual(
            result['states'][1]['orbital'],
            'd',
        )
        self.assertEqual(
            result['states'][1]['weights'],
            [
                [0.5, 0.6],
                [0.7, 0.8],
            ],
        )

    def test_parse_projwfc_band_file_normalizes_spin_down_indices(
        self,
    ):
        projection_text = """
        8 8 8 8 8 8 1 1
        1 5.0 0.0 0.0 0.0 0.0 0.0
        10.0 4.0 30.0 9
        1 Fe 16.0
        1 0.0 0.0 0.0 1
        1 2 1
        F F
        1 1 Fe 3D 1 2 1
        3 1 0.25
        4 1 0.75
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            projection_file = (
                Path(tmpdir)
                / 'bands.projwfc_down'
            )

            projection_file.write_text(
                projection_text,
                encoding='utf-8',
            )

            result = parse_projwfc_band_file(
                projection_file
            )

        self.assertEqual(
            result['states'][0]['weights'],
            [
                [0.25],
                [0.75],
            ],
        )

    def test_parse_projwfc_band_file_selects_requested_kpoints(self):
        projection_text = """
        8 8 8 8 8 8 1 1
        0 5.0 0.0 0.0 0.0 0.0 0.0
        10.0 4.0 30.0 9
        1 Fe 16.0
        1 0.0 0.0 0.0 1
        3 1 1
        F F
        1 1 Fe 4S 1 0 1
        1 1 0.10
        2 1 0.20
        3 1 0.30
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            projection_file = (
                Path(tmpdir)
                / 'bands.projwfc_up'
            )

            projection_file.write_text(
                projection_text,
                encoding='utf-8',
            )

            result = parse_projwfc_band_file(
                projection_file,
                kpoint_indices=[2, 0],
            )

        self.assertEqual(
            result['nkpoints'],
            2,
        )
        self.assertEqual(
            result['kpoint_indices'],
            [2, 0],
        )
        self.assertEqual(
            result['states'][0]['weights'],
            [
                [0.30],
                [0.10],
            ],
        )

    def test_parse_projwfc_band_file_rejects_spin_orbit(
        self,
    ):
        projection_text = """
        8 8 8 8 8 8 1 1
        1 5.0 0.0 0.0 0.0 0.0 0.0
        10.0 4.0 30.0 9
        1 Fe 16.0
        1 0.0 0.0 0.0 1
        1 1 1
        F T
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            projection_file = (
                Path(tmpdir)
                / 'bands.projwfc_up'
            )

            projection_file.write_text(
                projection_text,
                encoding='utf-8',
            )

            with self.assertRaisesRegex(
                NotImplementedError,
                'spin-orbit',
            ):
                parse_projwfc_band_file(
                    projection_file
                )

    def test_prepare_qe_band_projection_data(self):
        projection_result = {
            'natoms': 2,
            'nkpoints': 2,
            'nbands': 2,
            'states': [
                {
                    'atom_index': 0,
                    'orbital': 's',
                    'weights': [
                        [0.10, 0.20],
                        [0.30, 0.40],
                    ],
                },
                {
                    'atom_index': 0,
                    'orbital': 'd',
                    'weights': [
                        [0.50, 0.60],
                        [0.70, 0.80],
                    ],
                },
                {
                    'atom_index': 1,
                    'orbital': 'd',
                    'weights': [
                        [0.01, 0.02],
                        [0.03, 0.04],
                    ],
                },
            ],
        }

        result = prepare_qe_band_projection_data(
            projection_result,
            projections=[
                {
                    'atoms': [0],
                    'orbital': 'd',
                    'color': 'red',
                    'label': 'Fe d',
                },
                {
                    'atoms': [0, 1],
                    'orbital': 'd',
                    'color': 'green',
                    'label': 'Total d',
                },
            ],
        )

        first = result[
            'projections'
        ][0]

        second = result[
            'projections'
        ][1]

        self.assertEqual(
            first['selected_state_count'],
            1,
        )
        self.assertEqual(
            first['weights'],
            [
                [0.50, 0.60],
                [0.70, 0.80],
            ],
        )
        self.assertEqual(
            second['selected_state_count'],
            2,
        )
        self.assertEqual(
            second['weights'],
            [
                [0.51, 0.62],
                [0.73, 0.84],
            ],
        )

    def test_prepare_qe_band_projection_data_builds_total(
        self,
    ):
        projection_result = {
            'natoms': 2,
            'nkpoints': 1,
            'nbands': 2,
            'states': [
                {
                    'atom_index': 0,
                    'orbital': 's',
                    'weights': [
                        [0.10, 0.20],
                    ],
                },
                {
                    'atom_index': 1,
                    'orbital': 'p',
                    'weights': [
                        [0.30, 0.40],
                    ],
                },
            ],
        }

        result = prepare_qe_band_projection_data(
            projection_result,
            projections=[],
        )

        total = result[
            'projections'
        ][0]

        self.assertEqual(
            total['atoms'],
            [0, 1],
        )
        self.assertIsNone(
            total['orbital']
        )
        self.assertEqual(
            total['color'],
            'blue',
        )
        self.assertEqual(
            total['label'],
            'Total Contribution',
        )
        self.assertEqual(
            total['weights'],
            [
                [0.40, 0.60],
            ],
        )

    def test_prepare_qe_band_projection_data_rejects_atom_index(
        self,
    ):
        projection_result = {
            'natoms': 2,
            'nkpoints': 1,
            'nbands': 1,
            'states': [],
        }

        with self.assertRaisesRegex(
            ValueError,
            'atom index 2 is out of range',
        ):
            prepare_qe_band_projection_data(
                projection_result,
                projections=[
                    {
                        'atoms': [2],
                        'orbital': 'd',
                    },
                ],
            )

    def test_prepare_qe_band_projection_data_rejects_orbital(
        self,
    ):
        projection_result = {
            'natoms': 1,
            'nkpoints': 1,
            'nbands': 1,
            'states': [],
        }

        with self.assertRaisesRegex(
            ValueError,
            'Unsupported QE band projection orbital',
        ):
            prepare_qe_band_projection_data(
                projection_result,
                projections=[
                    {
                        'atoms': [0],
                        'orbital': 'g',
                    },
                ],
            )

    def test_run_bands_prepares_projected_band_data(self):
        atoms = bulk(
            'Si',
            'diamond',
            a=5.43,
        )

        band_path = build_band_path(
            atoms=atoms,
            path='GX',
            npoints=2,
        )

        output_text = """
         Program PWSCF v.7.2 starts

         End of band structure calculation

              k = 0.0000 0.0000 0.0000 ( 123 PWs)   bands (ev):

            -5.0000  -1.0000

              k = 0.5000 0.0000 0.5000 ( 120 PWs)   bands (ev):

            -4.5000  -0.5000

         JOB DONE.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(
                tmpdir
            )

            state_dir = (
                tmpdir
                / 'state'
            )

            save_dir = (
                state_dir
                / 'nanoworks.save'
            )

            save_dir.mkdir(
                parents=True
            )

            (
                save_dir
                / 'data-file-schema.xml'
            ).write_text(
                '<qes/>',
                encoding='utf-8',
            )

            def fake_run_qe_program(**kwargs):
                Path(
                    kwargs['output_file']
                ).write_text(
                    output_text,
                    encoding='utf-8',
                )

                return {
                    'returncode': 0,
                }

            projection_workflow = {
                'projection_up_file': (
                    tmpdir
                    / 'projection.projwfc_up'
                ),
                'projection_down_file': None,
            }

            raw_projection = {
                'natoms': 2,
                'nkpoints': 2,
                'nbands': 2,
                'states': [],
            }

            prepared_projection = {
                'natoms': 2,
                'nkpoints': 2,
                'nbands': 2,
                'projections': [],
            }

            with patch(
                'nanoworks.engine.qe.run_qe_program',
                side_effect=fake_run_qe_program,
            ), patch(
                'nanoworks.engine.qe.run_band_projections',
                return_value=projection_workflow,
            ) as run_projection, patch(
                'nanoworks.engine.qe.parse_projwfc_band_file',
                return_value=raw_projection,
            ) as parse_projection, patch(
                'nanoworks.engine.qe.prepare_qe_band_projection_data',
                return_value=prepared_projection,
            ) as prepare_projection:
                workflow = run_bands(
                    atoms=atoms,
                    input_file=tmpdir / 'bands.in',
                    output_file=tmpdir / 'bands.out',
                    state_dir=state_dir,
                    pseudopotentials={
                        'Si': 'Si.upf',
                    },
                    pseudo_dir='/tmp/pseudos',
                    cutoff_ev=400.0,
                    band_path=band_path,
                    nbands=2,
                    projected_band=True,
                    projections=[
                        {
                            'atoms': [0],
                            'orbital': 'p',
                            'color': 'red',
                            'label': 'Si p',
                        },
                    ],
                    projection_input_file=(
                        tmpdir
                        / 'projection.in'
                    ),
                    projection_output_file=(
                        tmpdir
                        / 'projection.out'
                    ),
                    projection_prefix=(
                        tmpdir
                        / 'projection'
                    ),
                )

        self.assertEqual(
            run_projection.call_count,
            1,
        )

        self.assertEqual(
            parse_projection.call_count,
            1,
        )

        prepare_projection.assert_called_once_with(
            raw_projection,
            projections=[
                {
                    'atoms': [0],
                    'orbital': 'p',
                    'color': 'red',
                    'label': 'Si p',
                },
            ],
        )

        self.assertFalse(
            workflow[
                'band_projections'
            ][
                'spin_polarized'
            ]
        )

        self.assertEqual(
            workflow[
                'band_projections'
            ][
                'up'
            ],
            prepared_projection,
        )

        self.assertIsNone(
            workflow[
                'band_projections'
            ][
                'down'
            ]
        )

    def test_run_bands_requires_projected_band_paths(self):
        atoms = bulk(
            'Si',
            'diamond',
            a=5.43,
        )

        band_path = build_band_path(
            atoms=atoms,
            path='GX',
            npoints=1,
        )

        output_text = """
         Program PWSCF v.7.2 starts

              k = 0.0000 0.0000 0.0000 ( 123 PWs)   bands (ev):

            -5.0000

         JOB DONE.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(
                tmpdir
            )

            state_dir = (
                tmpdir
                / 'state'
            )

            save_dir = (
                state_dir
                / 'nanoworks.save'
            )

            save_dir.mkdir(
                parents=True
            )

            (
                save_dir
                / 'data-file-schema.xml'
            ).write_text(
                '<qes/>',
                encoding='utf-8',
            )

            def fake_run_qe_program(**kwargs):
                Path(
                    kwargs['output_file']
                ).write_text(
                    output_text,
                    encoding='utf-8',
                )

                return {
                    'returncode': 0,
                }

            with patch(
                'nanoworks.engine.qe.run_qe_program',
                side_effect=fake_run_qe_program,
            ):
                with self.assertRaisesRegex(
                    ValueError,
                    'projection_input_file',
                ):
                    run_bands(
                        atoms=atoms,
                        input_file=(
                            tmpdir
                            / 'bands.in'
                        ),
                        output_file=(
                            tmpdir
                            / 'bands.out'
                        ),
                        state_dir=state_dir,
                        pseudopotentials={
                            'Si': 'Si.upf',
                        },
                        pseudo_dir='/tmp/pseudos',
                        cutoff_ev=400.0,
                        band_path=band_path,
                        nbands=1,
                        projected_band=True,
                    )

    def test_render_pp_total_density_input(self):
        text = render_pp_input(
            prefix='nanoworks',
            outdir='/tmp/qe-state',
            filplot='/tmp/total.pp',
            fileout='/tmp/total.cube',
            plot_num=0,
            spin_component=0,
        )

        self.assertIn(
            '&INPUTPP',
            text,
        )
        self.assertIn(
            "prefix = 'nanoworks'",
            text,
        )
        self.assertIn(
            "outdir = '/tmp/qe-state'",
            text,
        )
        self.assertIn(
            "filplot = '/tmp/total.pp'",
            text,
        )
        self.assertIn(
            'plot_num = 0',
            text,
        )
        self.assertIn(
            'spin_component = 0',
            text,
        )
        self.assertIn(
            '&PLOT',
            text,
        )
        self.assertIn(
            'iflag = 3',
            text,
        )
        self.assertIn(
            'output_format = 6',
            text,
        )
        self.assertIn(
            "fileout = '/tmp/total.cube'",
            text,
        )

    def test_render_pp_spin_density_input(self):
        text = render_pp_input(
            filplot='/tmp/spin.pp',
            fileout='/tmp/spin.cube',
            plot_num=6,
        )

        self.assertIn(
            'plot_num = 6',
            text,
        )
        self.assertNotIn(
            'spin_component',
            text,
        )

    def test_render_pp_input_rejects_unsupported_plot(self):
        with self.assertRaisesRegex(
            ValueError,
            'Unsupported QE density plot number',
        ):
            render_pp_input(
                plot_num=9,
            )

    def test_render_pp_input_rejects_invalid_spin_component(
        self,
    ):
        with self.assertRaisesRegex(
            ValueError,
            'must be 0, 1, or 2',
        ):
            render_pp_input(
                plot_num=0,
                spin_component=3,
            )

    def test_render_pp_input_rejects_spin_component_for_spin_density(
        self,
    ):
        with self.assertRaisesRegex(
            ValueError,
            'supported only',
        ):
            render_pp_input(
                plot_num=6,
                spin_component=1,
            )

    def test_run_pp_density_requires_ground_state(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(
                tmpdir
            )

            with self.assertRaisesRegex(
                FileNotFoundError,
                'valid QE ground-state result',
            ):
                run_pp_density(
                    input_file=(
                        tmpdir
                        / 'density.in'
                    ),
                    output_file=(
                        tmpdir
                        / 'density.out'
                    ),
                    state_dir=(
                        tmpdir
                        / 'state'
                    ),
                    filplot=(
                        tmpdir
                        / 'density.pp'
                    ),
                    cube_file=(
                        tmpdir
                        / 'density.cube'
                    ),
                )

    def test_run_pp_density_creates_cube_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(
                tmpdir
            )

            state_dir = (
                tmpdir
                / 'state'
            )

            save_dir = (
                state_dir
                / 'nanoworks.save'
            )

            save_dir.mkdir(
                parents=True
            )

            (
                save_dir
                / 'data-file-schema.xml'
            ).write_text(
                '<qes/>',
                encoding='utf-8',
            )

            cube_file = (
                tmpdir
                / 'density.cube'
            )

            def fake_run_qe_program(**kwargs):
                Path(
                    kwargs['output_file']
                ).write_text(
                    "Program PP\nJOB DONE.\n",
                    encoding='utf-8',
                )

                cube_file.write_text(
                    "Cube density data\n",
                    encoding='utf-8',
                )

                return {
                    'returncode': 0,
                }

            with patch(
                'nanoworks.engine.qe.run_qe_program',
                side_effect=fake_run_qe_program,
            ) as run:
                workflow = run_pp_density(
                    input_file=(
                        tmpdir
                        / 'density.in'
                    ),
                    output_file=(
                        tmpdir
                        / 'density.out'
                    ),
                    state_dir=state_dir,
                    filplot=(
                        tmpdir
                        / 'density.pp'
                    ),
                    cube_file=cube_file,
                    plot_num=0,
                    spin_component=1,
                    parallel_cores=4,
                )

            self.assertEqual(
                run.call_count,
                1,
            )

            self.assertEqual(
                run.call_args.kwargs[
                    'executable'
                ],
                'pp.x',
            )

            input_text = (
                workflow[
                    'input_file'
                ]
                .read_text(
                    encoding='utf-8'
                )
            )

            self.assertIn(
                'plot_num = 0',
                input_text,
            )

            self.assertIn(
                'spin_component = 1',
                input_text,
            )

            self.assertTrue(
                workflow[
                    'cube_file'
                ].is_file()
            )

            self.assertEqual(
                workflow[
                    'spin_component'
                ],
                1,
            )

    def test_run_pp_density_requires_cube_output(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(
                tmpdir
            )

            state_dir = (
                tmpdir
                / 'state'
            )

            save_dir = (
                state_dir
                / 'nanoworks.save'
            )

            save_dir.mkdir(
                parents=True
            )

            (
                save_dir
                / 'data-file-schema.xml'
            ).write_text(
                '<qes/>',
                encoding='utf-8',
            )

            def fake_run_qe_program(**kwargs):
                Path(
                    kwargs['output_file']
                ).write_text(
                    "Program PP\nJOB DONE.\n",
                    encoding='utf-8',
                )

                return {
                    'returncode': 0,
                }

            with patch(
                'nanoworks.engine.qe.run_qe_program',
                side_effect=fake_run_qe_program,
            ):
                with self.assertRaisesRegex(
                    RuntimeError,
                    'Cube file was not created',
                ):
                    run_pp_density(
                        input_file=(
                            tmpdir
                            / 'density.in'
                        ),
                        output_file=(
                            tmpdir
                            / 'density.out'
                        ),
                        state_dir=state_dir,
                        filplot=(
                            tmpdir
                            / 'density.pp'
                        ),
                        cube_file=(
                            tmpdir
                            / 'density.cube'
                        ),
                    )

    def test_resolve_and_render_qe_hubbard(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            pseudo_dir = Path(
                tmpdir
            )

            (
                pseudo_dir
                / 'O.upf'
            ).write_text(
                """
                <UPF version="2.0.1">
                <PP_PSWFC>
                  <PP_CHI.1 label="2S"></PP_CHI.1>
                  <PP_CHI.2 label="2P"></PP_CHI.2>
                </PP_PSWFC>
                </UPF>
                """,
                encoding='utf-8',
            )

            (
                pseudo_dir
                / 'Zn.upf'
            ).write_text(
                """
                <UPF version="2.0.1">
                <PP_PSWFC>
                  <PP_CHI.1 label="4S"></PP_CHI.1>
                  <PP_CHI.2 label="3D"></PP_CHI.2>
                </PP_PSWFC>
                </UPF>
                """,
                encoding='utf-8',
            )

            settings = resolve_qe_hubbard(
                setup_params={
                    'O': ':p,7.0',
                    'Zn': ':d,10.0',
                },
                pseudopotentials={
                    'O': 'O.upf',
                    'Zn': 'Zn.upf',
                },
                pseudo_dir=pseudo_dir,
            )

        text = render_qe_hubbard_card(
            settings
        )

        self.assertEqual(
            settings['projector'],
            'ortho-atomic',
        )
        self.assertIn(
            'HUBBARD (ortho-atomic)',
            text,
        )
        self.assertIn(
            'U O-2p 7',
            text,
        )
        self.assertIn(
            'U Zn-3d 10',
            text,
        )

    def test_render_qe_hubbard_expands_magnetic_species(self):
        settings = {
            'projector': 'ortho-atomic',
            'parameters': [
                {
                    'symbol': 'Fe',
                    'manifold': '3d',
                    'value_ev': 4.0,
                },
            ],
            'notices': [],
        }

        text = render_qe_hubbard_card(
            settings,
            species_by_element={
                'Fe': [
                    'Fe1',
                    'Fe2',
                ],
            },
        )

        self.assertIn(
            'U Fe1-3d 4',
            text,
        )
        self.assertIn(
            'U Fe2-3d 4',
            text,
        )

    def test_resolve_qe_hubbard_rejects_ambiguous_orbital(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            pseudo_dir = Path(
                tmpdir
            )

            (
                pseudo_dir
                / 'X.upf'
            ).write_text(
                """
                <UPF version="2.0.1">
                <PP_PSWFC>
                  <PP_CHI.1 label="2P"></PP_CHI.1>
                  <PP_CHI.2 label="3P"></PP_CHI.2>
                </PP_PSWFC>
                </UPF>
                """,
                encoding='utf-8',
            )

            with self.assertRaisesRegex(
                ValueError,
                'matches multiple pseudopotential manifolds',
            ):
                resolve_qe_hubbard(
                    setup_params={
                        'X': ':p,6.0',
                    },
                    pseudopotentials={
                        'X': 'X.upf',
                    },
                    pseudo_dir=pseudo_dir,
                )

    def test_render_pw_input_with_spin_and_hubbard(self):
        atoms = Atoms(
            'Fe2',
            positions=[
                [0.0, 0.0, 0.0],
                [1.5, 1.5, 1.5],
            ],
            cell=[
                [3.0, 0.0, 0.0],
                [0.0, 3.0, 0.0],
                [0.0, 0.0, 3.0],
            ],
            pbc=True,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            pseudo_dir = Path(
                tmpdir
            )

            (
                pseudo_dir
                / 'Fe.upf'
            ).write_text(
                """
                <UPF version="2.0.1">
                <PP_HEADER
                    element="Fe"
                    z_valence="8.0"
                />
                <PP_PSWFC>
                  <PP_CHI.1 label="4S"></PP_CHI.1>
                  <PP_CHI.2 label="3D"></PP_CHI.2>
                </PP_PSWFC>
                </UPF>
                """,
                encoding='utf-8',
            )

            text = render_pw_input(
                calculation='scf',
                atoms=atoms,
                pseudopotentials={
                    'Fe': 'Fe.upf',
                },
                pseudo_dir=pseudo_dir,
                cutoff_ev=500.0,
                kpoint_size=(4, 4, 4),
                spinpol=True,
                magnetic_moments=[
                    2.0,
                    -2.0,
                ],
                setup_params={
                    'Fe': ':d,4.0',
                },
            )

        self.assertIn(
            'Fe1 ',
            text,
        )
        self.assertIn(
            'Fe2 ',
            text,
        )
        self.assertIn(
            'HUBBARD (ortho-atomic)',
            text,
        )
        self.assertIn(
            'U Fe1-3d 4',
            text,
        )
        self.assertIn(
            'U Fe2-3d 4',
            text,
        )

    def test_pw_input_wrappers_forward_setup_params(self):
        atoms = bulk(
            'Si',
            'diamond',
            a=5.43,
        )

        setup_params = {
            'Si': ':p,4.0',
        }

        common = {
            'atoms': atoms,
            'pseudopotentials': {
                'Si': 'Si.upf',
            },
            'cutoff_ev': 400.0,
            'setup_params': setup_params,
        }

        with patch(
            'nanoworks.engine.qe.render_pw_input',
            return_value='rendered',
        ) as render:
            render_scf_input(
                **common,
                kpoint_size=(4, 4, 4),
            )

            self.assertEqual(
                render.call_args.kwargs[
                    'setup_params'
                ],
                setup_params,
            )

            render_nscf_input(
                **common,
                kpoint_size=(6, 6, 6),
            )

            self.assertEqual(
                render.call_args.kwargs[
                    'setup_params'
                ],
                setup_params,
            )

            render_relax_input(
                **common,
                kpoint_size=(4, 4, 4),
                optimizer='LBFGS',
                max_force=0.05,
                max_step=0.1,
                relax_cell=[
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                ],
            )

            self.assertEqual(
                render.call_args.kwargs[
                    'setup_params'
                ],
                setup_params,
            )

            render_bands_input(
                **common,
                band_path={
                    'kpoints': [
                        [
                            0.0,
                            0.0,
                            0.0,
                        ],
                    ],
                },
            )

            self.assertEqual(
                render.call_args.kwargs[
                    'setup_params'
                ],
                setup_params,
            )

    def test_run_scf_writes_hubbard_card(self):
        atoms = bulk(
            'Si',
            'diamond',
            a=5.43,
        )

        output_text = """
        Program PWSCF v.7.2 starts

        !    total energy = -15.00000000 Ry

        JOB DONE.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(
                tmpdir
            )

            pseudo_file = (
                tmpdir
                / 'Si.upf'
            )

            pseudo_file.write_text(
                """
                <UPF version="2.0.1">
                <PP_HEADER
                    element="Si"
                    z_valence="4.0"
                />
                <PP_PSWFC>
                  <PP_CHI.1 label="3S"></PP_CHI.1>
                  <PP_CHI.2 label="3P"></PP_CHI.2>
                </PP_PSWFC>
                </UPF>
                """,
                encoding='utf-8',
            )

            def fake_run_qe_program(**kwargs):
                Path(
                    kwargs['output_file']
                ).write_text(
                    output_text,
                    encoding='utf-8',
                )

                return {
                    'returncode': 0,
                }

            with patch(
                'nanoworks.engine.qe.run_qe_program',
                side_effect=fake_run_qe_program,
            ):
                workflow = run_scf(
                    atoms=atoms,
                    input_file=tmpdir / 'scf.in',
                    output_file=tmpdir / 'scf.out',
                    state_dir=tmpdir / 'state',
                    pseudopotentials={
                        'Si': 'Si.upf',
                    },
                    pseudo_dir=tmpdir,
                    cutoff_ev=400.0,
                    setup_params={
                        'Si': ':p,4.0',
                    },
                )

            input_text = (
                workflow['input_file']
                .read_text(
                    encoding='utf-8',
                )
            )

        self.assertIn(
            'HUBBARD (ortho-atomic)',
            input_text,
        )
        self.assertIn(
            'U Si-3p 4',
            input_text,
        )
        
    def test_parse_pw_bands_output_supports_adjacent_values(self):
        output_text = """
        Program PWSCF v.7.2 starts

             k = 0.0000 0.0000 0.0000 (123 PWs) bands (ev):

          -114.3470-114.3470  -20.5000   1.2500

        JOB DONE.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = (
                Path(tmpdir)
                / 'bands.out'
            )

            output_file.write_text(
                output_text,
                encoding='utf-8',
            )

            result = parse_pw_bands_output(
                output_file
            )

        self.assertEqual(
            result['nkpoints'],
            1,
        )
        self.assertEqual(
            result['nbands'],
            4,
        )
        self.assertEqual(
            result['eigenvalues_ev'][0][0],
            [
                -114.347,
                -114.347,
                -20.5,
                1.25,
            ],
        )

    def test_build_ph_settings_for_regular_grid(self):
        settings = build_ph_settings(
            prefix='nanoworks',
            outdir='/tmp/qe-state',
            fildyn='Si-PHONON-QE-Result-Dynamical-Matrix',
            qpoint_grid=(2, 3, 4),
            tr2_ph=1.0e-14,
        )

        self.assertEqual(settings['prefix'], 'nanoworks')
        self.assertEqual(settings['outdir'], '/tmp/qe-state')
        self.assertEqual(
            settings['fildyn'],
            'Si-PHONON-QE-Result-Dynamical-Matrix',
        )
        self.assertEqual(settings['tr2_ph'], 1.0e-14)
        self.assertTrue(settings['ldisp'])
        self.assertEqual(
            (
                settings['nq1'],
                settings['nq2'],
                settings['nq3'],
            ),
            (2, 3, 4),
        )

    def test_render_ph_input_for_regular_grid(self):
        text = render_ph_input(
            prefix='nanoworks',
            outdir='/tmp/qe-state',
            fildyn='Si-PHONON-QE-Result-Dynamical-Matrix',
            qpoint_grid=(2, 2, 2),
        )

        self.assertTrue(
            text.startswith(
                'Nanoworks native QE phonon calculation\n'
                '&INPUTPH\n'
            )
        )
        self.assertIn(
            "  prefix = 'nanoworks',",
            text,
        )
        self.assertIn(
            "  outdir = '/tmp/qe-state',",
            text,
        )
        self.assertIn(
            "  fildyn = "
            "'Si-PHONON-QE-Result-Dynamical-Matrix',",
            text,
        )
        self.assertIn(
            '  tr2_ph = 1e-12,',
            text,
        )
        self.assertIn(
            '  ldisp = .true.,',
            text,
        )
        self.assertIn('  nq1 = 2,', text)
        self.assertIn('  nq2 = 2,', text)
        self.assertIn('  nq3 = 2,', text)
        self.assertTrue(
            text.endswith('/\n')
        )

    def test_build_ph_settings_rejects_invalid_grid(self):
        for qpoint_grid in (
            (2, 2),
            (2, 0, 2),
            (2, 2.5, 2),
            (True, 2, 2),
        ):
            with self.subTest(
                qpoint_grid=qpoint_grid
            ):
                with self.assertRaises(
                    (TypeError, ValueError)
                ):
                    build_ph_settings(
                        prefix='nanoworks',
                        outdir='/tmp/qe-state',
                        fildyn='si.dyn',
                        qpoint_grid=qpoint_grid,
                    )

    def test_build_ph_settings_rejects_invalid_threshold(self):
        for tr2_ph in (
            0.0,
            -1.0e-12,
            float('inf'),
        ):
            with self.subTest(
                tr2_ph=tr2_ph
            ):
                with self.assertRaisesRegex(
                    ValueError,
                    'positive finite value',
                ):
                    build_ph_settings(
                        prefix='nanoworks',
                        outdir='/tmp/qe-state',
                        fildyn='si.dyn',
                        qpoint_grid=(2, 2, 2),
                        tr2_ph=tr2_ph,
                    )

    def test_build_ph_settings_rejects_empty_identifiers(self):
        for keyword, value in (
            ('prefix', None),
            ('outdir', ''),
            ('fildyn', '   '),
        ):
            arguments = {
                'prefix': 'nanoworks',
                'outdir': '/tmp/qe-state',
                'fildyn': 'si.dyn',
                'qpoint_grid': (2, 2, 2),
            }

            arguments[keyword] = value

            with self.subTest(keyword=keyword):
                with self.assertRaisesRegex(
                    ValueError,
                    f'{keyword} must not be empty',
                ):
                    build_ph_settings(
                        **arguments
                    )

    def test_render_ph_input_rejects_invalid_title(self):
        for title in (
            None,
            '',
            'line one\nline two',
        ):
            with self.subTest(title=title):
                with self.assertRaisesRegex(
                    ValueError,
                    'one non-empty line',
                ):
                    render_ph_input(
                        prefix='nanoworks',
                        outdir='/tmp/qe-state',
                        fildyn='si.dyn',
                        qpoint_grid=(2, 2, 2),
                        title=title,
                    )

    def test_build_q2r_settings(self):
        settings = build_q2r_settings(
            fildyn=(
                'Si-PHONON-QE-Result-'
                'Dynamical-Matrix'
            ),
            flfrc=(
                'Si-PHONON-QE-Result-'
                'Force-Constants.fc'
            ),
            zasr='crystal',
        )

        self.assertEqual(
            settings['fildyn'],
            (
                'Si-PHONON-QE-Result-'
                'Dynamical-Matrix'
            ),
        )
        self.assertEqual(
            settings['flfrc'],
            (
                'Si-PHONON-QE-Result-'
                'Force-Constants.fc'
            ),
        )
        self.assertEqual(
            settings['zasr'],
            'crystal',
        )

    def test_render_q2r_input(self):
        text = render_q2r_input(
            fildyn=(
                'Si-PHONON-QE-Result-'
                'Dynamical-Matrix'
            ),
            flfrc=(
                'Si-PHONON-QE-Result-'
                'Force-Constants.fc'
            ),
        )

        self.assertTrue(
            text.startswith('&INPUT\n')
        )
        self.assertIn(
            "  fildyn = "
            "'Si-PHONON-QE-Result-"
            "Dynamical-Matrix',",
            text,
        )
        self.assertIn(
            "  flfrc = "
            "'Si-PHONON-QE-Result-"
            "Force-Constants.fc',",
            text,
        )
        self.assertIn(
            "  zasr = 'no',",
            text,
        )
        self.assertTrue(
            text.endswith('/\n')
        )

    def test_build_q2r_settings_normalizes_zasr(self):
        settings = build_q2r_settings(
            fildyn='si.dyn',
            flfrc='si.fc',
            zasr=' CRYSTAL ',
        )

        self.assertEqual(
            settings['zasr'],
            'crystal',
        )

    def test_build_q2r_settings_rejects_empty_filenames(self):
        for keyword, value in (
            ('fildyn', None),
            ('fildyn', ''),
            ('flfrc', None),
            ('flfrc', '   '),
        ):
            arguments = {
                'fildyn': 'si.dyn',
                'flfrc': 'si.fc',
            }

            arguments[keyword] = value

            with self.subTest(
                keyword=keyword,
                value=value,
            ):
                with self.assertRaisesRegex(
                    ValueError,
                    f'{keyword} must not be empty',
                ):
                    build_q2r_settings(
                        **arguments
                    )

    def test_build_q2r_settings_rejects_invalid_zasr(self):
        for zasr in (
            None,
            '',
            'all',
            'invalid',
        ):
            with self.subTest(zasr=zasr):
                with self.assertRaisesRegex(
                    ValueError,
                    'zasr|Unsupported',
                ):
                    build_q2r_settings(
                        fildyn='si.dyn',
                        flfrc='si.fc',
                        zasr=zasr,
                    )

    def test_build_matdyn_band_settings(self):
        settings = build_matdyn_band_settings(
            flfrc='si.fc',
            flfrq='si.freq',
            acoustic_sum_rule=True,
        )

        self.assertEqual(
            settings['flfrc'],
            'si.fc',
        )
        self.assertEqual(
            settings['flfrq'],
            'si.freq',
        )
        self.assertEqual(
            settings['asr'],
            'crystal',
        )
        self.assertFalse(settings['dos'])
        self.assertFalse(
            settings['q_in_band_form']
        )
        self.assertTrue(
            settings['q_in_cryst_coord']
        )

    def test_build_matdyn_band_settings_can_disable_asr(self):
        settings = build_matdyn_band_settings(
            flfrc='si.fc',
            flfrq='si.freq',
            acoustic_sum_rule=False,
        )

        self.assertEqual(
            settings['asr'],
            'no',
        )

    def test_render_matdyn_qpoints(self):
        text = render_matdyn_qpoints(
            {
                'option': 'crystal',
                'kpoints': [
                    (0.0, 0.0, 0.0),
                    (0.0, 0.25, 0.0),
                    (0.0, 0.5, 0.0),
                ],
                'npoints': 3,
            }
        )

        self.assertEqual(
            text.splitlines(),
            [
                '3',
                '0.000000000000 0.000000000000 0.000000000000',
                '0.000000000000 0.250000000000 0.000000000000',
                '0.000000000000 0.500000000000 0.000000000000',
            ],
        )

    def test_render_matdyn_band_input(self):
        band_path = build_band_path(
            atoms=Atoms(
                'Si',
                cell=[4.0, 4.0, 4.0],
                pbc=True,
            ),
            path='GX',
            npoints=3,
        )

        text = render_matdyn_band_input(
            flfrc='si.fc',
            flfrq='si.freq',
            band_path=band_path,
        )

        self.assertIn(
            "  flfrc = 'si.fc',",
            text,
        )
        self.assertIn(
            "  asr = 'crystal',",
            text,
        )
        self.assertIn(
            '  dos = .false.,',
            text,
        )
        self.assertIn(
            "  flfrq = 'si.freq',",
            text,
        )
        self.assertIn(
            '  q_in_band_form = .false.,',
            text,
        )
        self.assertIn(
            '  q_in_cryst_coord = .true.,',
            text,
        )
        self.assertIn(
            '\n3\n',
            text,
        )
        self.assertTrue(
            text.endswith('\n')
        )

    def test_render_matdyn_qpoints_rejects_invalid_path(self):
        invalid_paths = (
            None,
            {},
            {
                'option': 'cartesian',
                'kpoints': [(0.0, 0.0, 0.0)],
            },
            {
                'option': 'crystal',
                'kpoints': [(0.0, 0.0)],
            },
            {
                'option': 'crystal',
                'kpoints': [(0.0, 0.0, 0.0)],
                'npoints': 2,
            },
        )

        for band_path in invalid_paths:
            with self.subTest(
                band_path=band_path
            ):
                with self.assertRaises(
                    (TypeError, ValueError)
                ):
                    render_matdyn_qpoints(
                        band_path
                    )

    def test_build_matdyn_band_settings_rejects_invalid_values(self):
        for arguments in (
            {
                'flfrc': '',
                'flfrq': 'si.freq',
            },
            {
                'flfrc': 'si.fc',
                'flfrq': None,
            },
            {
                'flfrc': 'si.fc',
                'flfrq': 'si.freq',
                'acoustic_sum_rule': 'yes',
            },
        ):
            with self.subTest(arguments=arguments):
                with self.assertRaises(
                    (TypeError, ValueError)
                ):
                    build_matdyn_band_settings(
                        **arguments
                    )

    def test_build_matdyn_dos_settings(self):
        settings = build_matdyn_dos_settings(
            flfrc='si.fc',
            fldos='si.dos',
            qpoint_grid=(20, 18, 16),
            acoustic_sum_rule=True,
        )

        self.assertEqual(
            settings['flfrc'],
            'si.fc',
        )
        self.assertEqual(
            settings['fldos'],
            'si.dos',
        )
        self.assertEqual(
            settings['asr'],
            'crystal',
        )
        self.assertTrue(settings['dos'])
        self.assertEqual(
            (
                settings['nk1'],
                settings['nk2'],
                settings['nk3'],
            ),
            (20, 18, 16),
        )

    def test_render_matdyn_dos_input(self):
        text = render_matdyn_dos_input(
            flfrc='si.fc',
            fldos='si.dos',
            qpoint_grid=(20, 20, 20),
        )

        self.assertIn(
            "  flfrc = 'si.fc',",
            text,
        )
        self.assertIn(
            "  asr = 'crystal',",
            text,
        )
        self.assertIn(
            '  dos = .true.,',
            text,
        )
        self.assertIn('  nk1 = 20,', text)
        self.assertIn('  nk2 = 20,', text)
        self.assertIn('  nk3 = 20,', text)
        self.assertIn(
            "  fldos = 'si.dos',",
            text,
        )
        self.assertTrue(
            text.endswith('/\n')
        )

    def test_build_matdyn_dos_settings_can_disable_asr(self):
        settings = build_matdyn_dos_settings(
            flfrc='si.fc',
            fldos='si.dos',
            qpoint_grid=(8, 8, 8),
            acoustic_sum_rule=False,
        )

        self.assertEqual(
            settings['asr'],
            'no',
        )

    def test_build_matdyn_dos_settings_rejects_invalid_values(self):
        invalid_arguments = (
            {
                'flfrc': '',
                'fldos': 'si.dos',
                'qpoint_grid': (20, 20, 20),
            },
            {
                'flfrc': 'si.fc',
                'fldos': None,
                'qpoint_grid': (20, 20, 20),
            },
            {
                'flfrc': 'si.fc',
                'fldos': 'si.dos',
                'qpoint_grid': (20, 20),
            },
            {
                'flfrc': 'si.fc',
                'fldos': 'si.dos',
                'qpoint_grid': (20, 0, 20),
            },
            {
                'flfrc': 'si.fc',
                'fldos': 'si.dos',
                'qpoint_grid': (20, 2.5, 20),
            },
            {
                'flfrc': 'si.fc',
                'fldos': 'si.dos',
                'qpoint_grid': (20, 20, 20),
                'acoustic_sum_rule': 'yes',
            },
        )

        for arguments in invalid_arguments:
            with self.subTest(arguments=arguments):
                with self.assertRaises(
                    (TypeError, ValueError)
                ):
                    build_matdyn_dos_settings(
                        **arguments
                    )

    def test_parse_qe_auxiliary_output(self):
        output_text = """
        Program PHONON v.7.2 starts

        JOB DONE.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = (
                Path(tmpdir)
                / 'ph.out'
            )

            output_file.write_text(
                output_text,
                encoding='utf-8',
            )

            result = parse_qe_auxiliary_output(
                output_file,
                expected_program='PHONON',
            )

        self.assertEqual(
            result['program'],
            'PHONON',
        )
        self.assertEqual(
            result['qe_version'],
            (7, 2),
        )
        self.assertTrue(
            result['job_done']
        )

    def test_parse_qe_auxiliary_output_supports_patch_version(self):
        output_text = """
        Program MATDYN v.7.2.1 starts

        JOB DONE.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = (
                Path(tmpdir)
                / 'matdyn.out'
            )

            output_file.write_text(
                output_text,
                encoding='utf-8',
            )

            result = parse_qe_auxiliary_output(
                output_file,
                expected_program='matdyn',
            )

        self.assertEqual(
            result['qe_version'],
            (7, 2, 1),
        )

    def test_parse_qe_auxiliary_output_preserves_incomplete_job(self):
        output_text = """
        Program Q2R v.7.2 starts
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = (
                Path(tmpdir)
                / 'q2r.out'
            )

            output_file.write_text(
                output_text,
                encoding='utf-8',
            )

            result = parse_qe_auxiliary_output(
                output_file,
                expected_program='Q2R',
            )

        self.assertFalse(
            result['job_done']
        )

    def test_parse_qe_auxiliary_output_rejects_wrong_program(self):
        output_text = """
        Program Q2R v.7.2 starts

        JOB DONE.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = (
                Path(tmpdir)
                / 'q2r.out'
            )

            output_file.write_text(
                output_text,
                encoding='utf-8',
            )

            with self.assertRaisesRegex(
                ValueError,
                'expected MATDYN, found Q2R',
            ):
                parse_qe_auxiliary_output(
                    output_file,
                    expected_program='MATDYN',
                )

    def test_parse_qe_auxiliary_output_requires_output_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = (
                Path(tmpdir)
                / 'missing.out'
            )

            with self.assertRaises(
                FileNotFoundError
            ):
                parse_qe_auxiliary_output(
                    output_file,
                    expected_program='PHONON',
                )

    def test_parse_matdyn_frequency_file(self):
        frequency_text = """
 &plot nbnd=   6, nks=   2 /
          0.000000  0.000000  0.000000
  -12.5000-11.2500  -0.5000   5.0000  10.0000  15.0000
          0.000000  0.500000  0.000000
   1.0000D+01  20.0000  30.0000  40.0000  50.0000  60.0000
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            frequency_file = (
                Path(tmpdir)
                / 'si.freq'
            )
            frequency_file.write_text(
                frequency_text,
                encoding='utf-8',
            )

            result = parse_matdyn_frequency_file(
                frequency_file
            )

        self.assertEqual(
            result['nqpoints'],
            2,
        )
        self.assertEqual(
            result['nmodes'],
            6,
        )
        self.assertEqual(
            result['qpoints'][1],
            (0.0, 0.5, 0.0),
        )
        self.assertEqual(
            result['frequencies_cm1'][0][:3],
            [-12.5, -11.25, -0.5],
        )
        self.assertEqual(
            result['frequencies_cm1'][1][0],
            10.0,
        )
        self.assertAlmostEqual(
            result['frequencies_thz'][0][0],
            -12.5 * THZ_PER_CM_MINUS_ONE,
        )

    def test_parse_matdyn_frequency_file_rejects_bad_count(self):
        frequency_text = """
 &plot nbnd=   3, nks=   1 /
          0.000000  0.000000  0.000000
   1.0000   2.0000
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            frequency_file = (
                Path(tmpdir)
                / 'bad.freq'
            )
            frequency_file.write_text(
                frequency_text,
                encoding='utf-8',
            )

            with self.assertRaisesRegex(
                ValueError,
                'data count does not match',
            ):
                parse_matdyn_frequency_file(
                    frequency_file
                )

    def test_parse_matdyn_frequency_file_requires_header(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            frequency_file = (
                Path(tmpdir)
                / 'bad.freq'
            )
            frequency_file.write_text(
                'not a matdyn frequency file',
                encoding='utf-8',
            )

            with self.assertRaisesRegex(
                ValueError,
                'header could not be parsed',
            ):
                parse_matdyn_frequency_file(
                    frequency_file
                )

    def test_parse_matdyn_frequency_file_requires_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaises(
                FileNotFoundError
            ):
                parse_matdyn_frequency_file(
                    Path(tmpdir)
                    / 'missing.freq'
                )

    def test_parse_matdyn_dos_file(self):
        dos_text = """
# Frequency[cm^-1] DOS PDOS
-1.0000000000E+01 2.5000000000E-03 1.0000E-03 1.5000E-03
 0.0000000000D+00 4.0000000000D-03 1.7500D-03 2.2500D-03
 1.0000000000E+01 3.0000000000E-03 1.2500E-03 1.7500E-03
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            dos_file = (
                Path(tmpdir)
                / 'si.phdos'
            )
            dos_file.write_text(
                dos_text,
                encoding='utf-8',
            )

            result = parse_matdyn_dos_file(
                dos_file
            )

        self.assertEqual(
            result['npoints'],
            3,
        )
        self.assertEqual(
            result['natoms'],
            2,
        )
        self.assertEqual(
            result['frequencies_cm1'],
            [-10.0, 0.0, 10.0],
        )
        self.assertAlmostEqual(
            result['frequencies_thz'][2],
            10.0 * THZ_PER_CM_MINUS_ONE,
        )
        self.assertEqual(
            result['dos'],
            [0.0025, 0.004, 0.003],
        )
        self.assertEqual(
            result['atom_projected_dos'][0],
            [0.001, 0.00175, 0.00125],
        )

    def test_parse_matdyn_dos_file_supports_adjacent_values(self):
        dos_text = (
            "# Frequency[cm^-1] DOS PDOS\n"
            "-1.0000000000E+01-2.5000000000E-03"
            "-1.0000E-03-1.5000E-03\n"
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            dos_file = (
                Path(tmpdir)
                / 'adjacent.phdos'
            )
            dos_file.write_text(
                dos_text,
                encoding='utf-8',
            )

            result = parse_matdyn_dos_file(
                dos_file
            )

        self.assertEqual(
            result['frequencies_cm1'],
            [-10.0],
        )
        self.assertEqual(
            result['dos'],
            [-0.0025],
        )
        self.assertEqual(
            result['atom_projected_dos'],
            [[-0.001], [-0.0015]],
        )

    def test_parse_matdyn_dos_file_rejects_inconsistent_columns(self):
        dos_text = (
            "# Frequency[cm^-1] DOS PDOS\n"
            "0.0 1.0 0.4 0.6\n"
            "1.0 2.0 2.0\n"
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            dos_file = (
                Path(tmpdir)
                / 'bad.phdos'
            )
            dos_file.write_text(
                dos_text,
                encoding='utf-8',
            )

            with self.assertRaisesRegex(
                ValueError,
                'inconsistent column counts',
            ):
                parse_matdyn_dos_file(
                    dos_file
                )

    def test_parse_matdyn_dos_file_rejects_missing_data(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            dos_file = (
                Path(tmpdir)
                / 'empty.phdos'
            )
            dos_file.write_text(
                '# Frequency[cm^-1] DOS PDOS\n',
                encoding='utf-8',
            )

            with self.assertRaisesRegex(
                ValueError,
                'No QE matdyn DOS data',
            ):
                parse_matdyn_dos_file(
                    dos_file
                )

    def test_parse_matdyn_dos_file_requires_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaises(
                FileNotFoundError
            ):
                parse_matdyn_dos_file(
                    Path(tmpdir)
                    / 'missing.phdos'
                )

    def test_write_matdyn_band_data(self):
        band_path = {
            'distances': [0.0, 0.25],
        }
        frequencies = {
            'qpoints': [
                (0.0, 0.0, 0.0),
                (0.5, 0.0, 0.0),
            ],
            'frequencies_thz': [
                [-0.1, 1.0, 2.0],
                [0.2, 1.5, 2.5],
            ],
            'nqpoints': 2,
            'nmodes': 3,
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = (
                Path(tmpdir)
                / 'band-thz.dat'
            )
            result = write_matdyn_band_data(
                output_file,
                band_path,
                frequencies,
            )
            lines = output_file.read_text(
                encoding='utf-8'
            ).splitlines()

        self.assertEqual(
            result,
            output_file,
        )
        self.assertIn(
            'Frequency_3(THz)',
            lines[0],
        )
        self.assertEqual(
            lines[1],
            '0.0000000000 0.0000000000 0.0000000000 '
            '0.0000000000 -0.1000000000 1.0000000000 '
            '2.0000000000',
        )

    def test_write_matdyn_band_data_rejects_mismatched_counts(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaisesRegex(
                ValueError,
                'counts do not match',
            ):
                write_matdyn_band_data(
                    Path(tmpdir)
                    / 'bad-band.dat',
                    {'distances': [0.0, 1.0]},
                    {
                        'qpoints': [(0.0, 0.0, 0.0)],
                        'frequencies_thz': [[1.0, 2.0, 3.0]],
                        'nqpoints': 1,
                        'nmodes': 3,
                    },
                )

    def test_write_matdyn_dos_data_converts_density_units(self):
        dos_data = {
            'frequencies_thz': [0.0, 1.0],
            'dos': [
                2.0 * THZ_PER_CM_MINUS_ONE,
                3.0 * THZ_PER_CM_MINUS_ONE,
            ],
            'atom_projected_dos': [
                [
                    0.5 * THZ_PER_CM_MINUS_ONE,
                    1.0 * THZ_PER_CM_MINUS_ONE,
                ],
                [
                    1.5 * THZ_PER_CM_MINUS_ONE,
                    2.0 * THZ_PER_CM_MINUS_ONE,
                ],
            ],
            'npoints': 2,
            'natoms': 2,
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = (
                Path(tmpdir)
                / 'dos-thz.dat'
            )
            result = write_matdyn_dos_data(
                output_file,
                dos_data,
            )
            lines = output_file.read_text(
                encoding='utf-8'
            ).splitlines()

        self.assertEqual(
            result,
            output_file,
        )
        self.assertIn(
            'Atom_2_PDOS(1/THz)',
            lines[0],
        )
        self.assertEqual(
            lines[2],
            '1.0000000000 3.0000000000 1.0000000000 '
            '2.0000000000',
        )

    def test_write_matdyn_dos_data_rejects_bad_projection_shape(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaisesRegex(
                ValueError,
                'dimensions do not match',
            ):
                write_matdyn_dos_data(
                    Path(tmpdir)
                    / 'bad-dos.dat',
                    {
                        'frequencies_thz': [0.0, 1.0],
                        'dos': [1.0, 2.0],
                        'atom_projected_dos': [[0.5]],
                        'npoints': 2,
                        'natoms': 1,
                    },
                )

    def test_calculate_phonon_thermal_properties(self):
        dos_data = {
            'frequencies_thz': [1.0, 2.0, 3.0],
            'dos': [
                0.0,
                2.0 * THZ_PER_CM_MINUS_ONE,
                0.0,
            ],
        }

        result = calculate_phonon_thermal_properties(
            dos_data,
            t_min=0.0,
            t_max=300.0,
            t_step=300.0,
        )

        expected_zpe = (
            2.0
            * 0.5
            * 2.0
            * 4.135667696e-3
            * 96.48533212331002
        )

        self.assertEqual(
            result['temperatures_k'],
            [0.0, 300.0],
        )
        self.assertAlmostEqual(
            result['integrated_mode_weight'],
            2.0,
        )
        self.assertAlmostEqual(
            result['zero_point_energy_kj_mol'],
            expected_zpe,
        )
        self.assertAlmostEqual(
            result['free_energy_kj_mol'][0],
            expected_zpe,
        )
        self.assertEqual(
            result['entropy_j_k_mol'][0],
            0.0,
        )
        self.assertEqual(
            result['heat_capacity_j_k_mol'][0],
            0.0,
        )
        self.assertGreater(
            result['entropy_j_k_mol'][1],
            0.0,
        )
        self.assertGreater(
            result['heat_capacity_j_k_mol'][1],
            0.0,
        )

    def test_calculate_phonon_thermal_properties_rejects_negative_dos(self):
        with self.assertRaisesRegex(
            ValueError,
            'must not be negative',
        ):
            calculate_phonon_thermal_properties(
                {
                    'frequencies_thz': [1.0, 2.0],
                    'dos': [0.0, -1.0],
                }
            )

    def test_calculate_phonon_thermal_properties_rejects_bad_range(self):
        with self.assertRaisesRegex(
            ValueError,
            '0 <= t_min',
        ):
            calculate_phonon_thermal_properties(
                {
                    'frequencies_thz': [1.0, 2.0],
                    'dos': [0.0, 1.0],
                },
                t_min=300.0,
                t_max=100.0,
            )

    def test_write_phonon_thermal_properties(self):
        thermal_data = {
            'temperatures_k': [0.0, 300.0],
            'free_energy_kj_mol': [1.0, -2.0],
            'internal_energy_kj_mol': [1.0, 3.0],
            'entropy_j_k_mol': [0.0, 10.0],
            'heat_capacity_j_k_mol': [0.0, 8.0],
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = (
                Path(tmpdir)
                / 'thermal.csv'
            )
            result = write_phonon_thermal_properties(
                output_file,
                thermal_data,
            )
            lines = output_file.read_text(
                encoding='utf-8'
            ).splitlines()

        self.assertEqual(
            result,
            output_file,
        )
        self.assertEqual(
            lines[0],
            'T(K),Free_Energy(kJ/mol),Internal_Energy(kJ/mol),'
            'Entropy(J/K/mol),Cv(J/K/mol)',
        )
        self.assertEqual(
            lines[2],
            '300.0000000000,-2.0000000000,3.0000000000,'
            '10.0000000000,8.0000000000',
        )

    def test_run_ph_uses_existing_ground_state(self):
        output_text = """
        Program PHONON v.7.2 starts

        JOB DONE.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            state_dir = tmpdir / 'ground-state'
            save_dir = state_dir / 'nanoworks.save'
            save_dir.mkdir(parents=True)
            (
                save_dir
                / 'data-file-schema.xml'
            ).write_text(
                '<espresso/>',
                encoding='utf-8',
            )

            fildyn = tmpdir / 'phonon' / 'si.dyn'

            def fake_run_qe_program(**kwargs):
                Path(
                    kwargs['output_file']
                ).write_text(
                    output_text,
                    encoding='utf-8',
                )
                Path(
                    str(fildyn) + '0'
                ).write_text(
                    'q-grid metadata',
                    encoding='utf-8',
                )
                Path(
                    str(fildyn) + '1'
                ).write_text(
                    'dynamical matrix',
                    encoding='utf-8',
                )

                return {
                    'returncode': 0,
                }

            with patch(
                'nanoworks.engine.qe.run_qe_program',
                side_effect=fake_run_qe_program,
            ):
                workflow = run_ph(
                    input_file=tmpdir / 'ph.in',
                    output_file=tmpdir / 'ph.out',
                    state_dir=state_dir,
                    fildyn=fildyn,
                    qpoint_grid=(2, 2, 2),
                )

            input_text = workflow[
                'input_file'
            ].read_text(
                encoding='utf-8'
            )

        self.assertIn(
            "  prefix = 'nanoworks',",
            input_text,
        )
        self.assertIn(
            "  ldisp = .true.,",
            input_text,
        )
        self.assertEqual(
            workflow['result']['program'],
            'PHONON',
        )
        self.assertEqual(
            len(workflow['dynamical_matrix_files']),
            1,
        )

    def test_run_ph_requires_ground_state(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            with self.assertRaisesRegex(
                FileNotFoundError,
                'ground-state directory',
            ):
                run_ph(
                    input_file=tmpdir / 'ph.in',
                    output_file=tmpdir / 'ph.out',
                    state_dir=tmpdir / 'missing-state',
                    fildyn=tmpdir / 'si.dyn',
                    qpoint_grid=(2, 2, 2),
                )

    def test_run_ph_requires_dynamical_matrix_outputs(self):
        output_text = """
        Program PHONON v.7.2 starts

        JOB DONE.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            state_dir = tmpdir / 'ground-state'
            save_dir = state_dir / 'nanoworks.save'
            save_dir.mkdir(parents=True)
            (
                save_dir
                / 'data-file-schema.xml'
            ).write_text(
                '<espresso/>',
                encoding='utf-8',
            )

            def fake_run_qe_program(**kwargs):
                Path(
                    kwargs['output_file']
                ).write_text(
                    output_text,
                    encoding='utf-8',
                )

                return {
                    'returncode': 0,
                }

            with patch(
                'nanoworks.engine.qe.run_qe_program',
                side_effect=fake_run_qe_program,
            ):
                with self.assertRaisesRegex(
                    RuntimeError,
                    'q-grid metadata file',
                ):
                    run_ph(
                        input_file=tmpdir / 'ph.in',
                        output_file=tmpdir / 'ph.out',
                        state_dir=state_dir,
                        fildyn=tmpdir / 'si.dyn',
                        qpoint_grid=(2, 2, 2),
                    )

    def test_run_q2r_uses_phonon_grid_outputs(self):
        output_text = """
        Program Q2R v.7.2 starts

        JOB DONE.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            fildyn = tmpdir / 'si.dyn'
            flfrc = tmpdir / 'si.fc'
            Path(
                str(fildyn) + '0'
            ).write_text(
                'q-grid metadata',
                encoding='utf-8',
            )

            def fake_run_qe_program(**kwargs):
                Path(
                    kwargs['output_file']
                ).write_text(
                    output_text,
                    encoding='utf-8',
                )
                flfrc.write_text(
                    'real-space force constants',
                    encoding='utf-8',
                )

                return {
                    'returncode': 0,
                }

            with patch(
                'nanoworks.engine.qe.run_qe_program',
                side_effect=fake_run_qe_program,
            ):
                workflow = run_q2r(
                    input_file=tmpdir / 'q2r.in',
                    output_file=tmpdir / 'q2r.out',
                    fildyn=fildyn,
                    flfrc=flfrc,
                )

            input_text = workflow[
                'input_file'
            ].read_text(
                encoding='utf-8'
            )

        self.assertIn(
            "  fildyn = '",
            input_text,
        )
        self.assertIn(
            "  zasr = 'no',",
            input_text,
        )
        self.assertEqual(
            workflow['result']['program'],
            'Q2R',
        )
        self.assertEqual(
            workflow['flfrc'],
            flfrc,
        )

    def test_run_q2r_requires_grid_metadata(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            with self.assertRaisesRegex(
                FileNotFoundError,
                'q-grid metadata file',
            ):
                run_q2r(
                    input_file=tmpdir / 'q2r.in',
                    output_file=tmpdir / 'q2r.out',
                    fildyn=tmpdir / 'si.dyn',
                    flfrc=tmpdir / 'si.fc',
                )

    def test_run_q2r_requires_force_constants_output(self):
        output_text = """
        Program Q2R v.7.2 starts

        JOB DONE.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            fildyn = tmpdir / 'si.dyn'
            Path(
                str(fildyn) + '0'
            ).write_text(
                'q-grid metadata',
                encoding='utf-8',
            )

            def fake_run_qe_program(**kwargs):
                Path(
                    kwargs['output_file']
                ).write_text(
                    output_text,
                    encoding='utf-8',
                )

                return {
                    'returncode': 0,
                }

            with patch(
                'nanoworks.engine.qe.run_qe_program',
                side_effect=fake_run_qe_program,
            ):
                with self.assertRaisesRegex(
                    RuntimeError,
                    'force-constant file',
                ):
                    run_q2r(
                        input_file=tmpdir / 'q2r.in',
                        output_file=tmpdir / 'q2r.out',
                        fildyn=fildyn,
                        flfrc=tmpdir / 'si.fc',
                    )

    def test_run_matdyn_band_uses_force_constants(self):
        output_text = """
        Program MATDYN v.7.2 starts

        JOB DONE.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            flfrc = tmpdir / 'si.fc'
            flfrq = tmpdir / 'si.freq'
            flfrc.write_text(
                'real-space force constants',
                encoding='utf-8',
            )

            band_path = {
                'option': 'crystal',
                'kpoints': [
                    (0.0, 0.0, 0.0),
                    (0.0, 0.5, 0.0),
                ],
                'npoints': 2,
            }

            def fake_run_qe_program(**kwargs):
                Path(
                    kwargs['output_file']
                ).write_text(
                    output_text,
                    encoding='utf-8',
                )
                flfrq.write_text(
                    " &plot nbnd=   3, nks=   2 /\n"
                    " 0.000000 0.000000 0.000000\n"
                    " 1.0000 2.0000 3.0000\n"
                    " 0.000000 0.500000 0.000000\n"
                    " 4.0000 5.0000 6.0000\n",
                    encoding='utf-8',
                )

                return {
                    'returncode': 0,
                }

            with patch(
                'nanoworks.engine.qe.run_qe_program',
                side_effect=fake_run_qe_program,
            ):
                workflow = run_matdyn_band(
                    input_file=tmpdir / 'matdyn-band.in',
                    output_file=tmpdir / 'matdyn-band.out',
                    flfrc=flfrc,
                    flfrq=flfrq,
                    band_path=band_path,
                )

            input_text = workflow[
                'input_file'
            ].read_text(
                encoding='utf-8'
            )

        self.assertIn(
            '  dos = .false.,',
            input_text,
        )
        self.assertIn(
            '  q_in_cryst_coord = .true.,',
            input_text,
        )
        self.assertEqual(
            workflow['result']['program'],
            'MATDYN',
        )
        self.assertEqual(
            workflow['flfrq'],
            flfrq,
        )
        self.assertEqual(
            workflow['frequencies']['nqpoints'],
            2,
        )
        self.assertEqual(
            workflow['frequencies']['frequencies_cm1'][1],
            [4.0, 5.0, 6.0],
        )

    def test_run_matdyn_band_requires_force_constants(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            with self.assertRaisesRegex(
                FileNotFoundError,
                'force-constant file',
            ):
                run_matdyn_band(
                    input_file=tmpdir / 'matdyn-band.in',
                    output_file=tmpdir / 'matdyn-band.out',
                    flfrc=tmpdir / 'missing.fc',
                    flfrq=tmpdir / 'si.freq',
                    band_path={
                        'option': 'crystal',
                        'kpoints': [
                            (0.0, 0.0, 0.0),
                        ],
                        'npoints': 1,
                    },
                )

    def test_run_matdyn_band_requires_frequency_output(self):
        output_text = """
        Program MATDYN v.7.2 starts

        JOB DONE.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            flfrc = tmpdir / 'si.fc'
            flfrc.write_text(
                'real-space force constants',
                encoding='utf-8',
            )

            def fake_run_qe_program(**kwargs):
                Path(
                    kwargs['output_file']
                ).write_text(
                    output_text,
                    encoding='utf-8',
                )

                return {
                    'returncode': 0,
                }

            with patch(
                'nanoworks.engine.qe.run_qe_program',
                side_effect=fake_run_qe_program,
            ):
                with self.assertRaisesRegex(
                    RuntimeError,
                    'phonon frequency file',
                ):
                    run_matdyn_band(
                        input_file=tmpdir / 'matdyn-band.in',
                        output_file=tmpdir / 'matdyn-band.out',
                        flfrc=flfrc,
                        flfrq=tmpdir / 'si.freq',
                        band_path={
                            'option': 'crystal',
                            'kpoints': [
                                (0.0, 0.0, 0.0),
                            ],
                            'npoints': 1,
                        },
                    )

    def test_run_matdyn_dos_uses_force_constants(self):
        output_text = """
        Program MATDYN v.7.2 starts

        JOB DONE.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            flfrc = tmpdir / 'si.fc'
            fldos = tmpdir / 'si.dos'
            flfrc.write_text(
                'real-space force constants',
                encoding='utf-8',
            )

            def fake_run_qe_program(**kwargs):
                Path(
                    kwargs['output_file']
                ).write_text(
                    output_text,
                    encoding='utf-8',
                )
                fldos.write_text(
                    "# Frequency[cm^-1] DOS PDOS\n"
                    "0.0 1.0 0.4 0.6\n"
                    "1.0 2.0 0.8 1.2\n",
                    encoding='utf-8',
                )

                return {
                    'returncode': 0,
                }

            with patch(
                'nanoworks.engine.qe.run_qe_program',
                side_effect=fake_run_qe_program,
            ):
                workflow = run_matdyn_dos(
                    input_file=tmpdir / 'matdyn-dos.in',
                    output_file=tmpdir / 'matdyn-dos.out',
                    flfrc=flfrc,
                    fldos=fldos,
                    qpoint_grid=(20, 20, 20),
                )

            input_text = workflow[
                'input_file'
            ].read_text(
                encoding='utf-8'
            )

        self.assertIn(
            '  dos = .true.,',
            input_text,
        )
        self.assertIn(
            '  nk1 = 20,',
            input_text,
        )
        self.assertEqual(
            workflow['result']['program'],
            'MATDYN',
        )
        self.assertEqual(
            workflow['qpoint_grid'],
            (20, 20, 20),
        )
        self.assertEqual(
            workflow['dos']['frequencies_cm1'],
            [0.0, 1.0],
        )
        self.assertEqual(
            workflow['dos']['atom_projected_dos'],
            [[0.4, 0.8], [0.6, 1.2]],
        )

    def test_run_matdyn_dos_requires_force_constants(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            with self.assertRaisesRegex(
                FileNotFoundError,
                'force-constant file',
            ):
                run_matdyn_dos(
                    input_file=tmpdir / 'matdyn-dos.in',
                    output_file=tmpdir / 'matdyn-dos.out',
                    flfrc=tmpdir / 'missing.fc',
                    fldos=tmpdir / 'si.dos',
                    qpoint_grid=(20, 20, 20),
                )

    def test_run_matdyn_dos_requires_dos_output(self):
        output_text = """
        Program MATDYN v.7.2 starts

        JOB DONE.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            flfrc = tmpdir / 'si.fc'
            flfrc.write_text(
                'real-space force constants',
                encoding='utf-8',
            )

            def fake_run_qe_program(**kwargs):
                Path(
                    kwargs['output_file']
                ).write_text(
                    output_text,
                    encoding='utf-8',
                )

                return {
                    'returncode': 0,
                }

            with patch(
                'nanoworks.engine.qe.run_qe_program',
                side_effect=fake_run_qe_program,
            ):
                with self.assertRaisesRegex(
                    RuntimeError,
                    'phonon DOS file',
                ):
                    run_matdyn_dos(
                        input_file=tmpdir / 'matdyn-dos.in',
                        output_file=tmpdir / 'matdyn-dos.out',
                        flfrc=flfrc,
                        fldos=tmpdir / 'si.dos',
                        qpoint_grid=(20, 20, 20),
                    )

if __name__ == '__main__':
    unittest.main()
