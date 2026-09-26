"""Engine-neutral SCF settings and backend adapters."""

SCF_ACCURACY_LEVELS = (
    'loose',
    'normal',
    'tight',
    'very-tight',
)

ELECTRONIC_SOLVER_PROFILES = (
    'default',
    'fast',
    'robust',
)


def validate_scf_settings(
    accuracy='normal',
    max_steps=None,
    mixing=None,
    solver='default',
):
    """Validate the portable SCF controls used by every DFT backend."""
    accuracy = str(accuracy).strip().lower()
    solver = str(solver).strip().lower()

    if accuracy not in SCF_ACCURACY_LEVELS:
        raise ValueError(
            "SCF_accuracy must be one of: "
            + ', '.join(SCF_ACCURACY_LEVELS)
            + '.'
        )

    if solver not in ELECTRONIC_SOLVER_PROFILES:
        raise ValueError(
            "Electronic_solver must be one of: "
            + ', '.join(ELECTRONIC_SOLVER_PROFILES)
            + '.'
        )

    if max_steps is not None:
        max_steps = int(max_steps)

        if max_steps <= 0:
            raise ValueError(
                "SCF_max_steps must be a positive integer."
            )

    if mixing is not None:
        mixing = float(mixing)

        if not 0.0 < mixing <= 1.0:
            raise ValueError(
                "SCF_mixing must be greater than zero and at most one."
            )

    return {
        'accuracy': accuracy,
        'max_steps': max_steps,
        'mixing': mixing,
        'solver': solver,
    }


def resolve_qe_scf_settings(**settings):
    """Translate portable SCF intent to QE ``&ELECTRONS`` settings."""
    resolved = validate_scf_settings(**settings)
    convergence = {
        'loose': 1.0e-4,
        'normal': None,
        'tight': 1.0e-8,
        'very-tight': 1.0e-10,
    }
    diagonalization = {
        'default': None,
        'fast': 'david',
        'robust': 'cg',
    }

    return {
        'conv_thr': convergence[resolved['accuracy']],
        'mixing_beta': resolved['mixing'],
        'electron_maxstep': resolved['max_steps'],
        'diagonalization': diagonalization[resolved['solver']],
    }


def resolve_gpaw_scf_settings(**settings):
    """Translate portable SCF intent to GPAW calculator settings."""
    resolved = validate_scf_settings(**settings)
    convergence = {
        'loose': {
            'energy': 1.0e-3,
            'density': 1.0e-2,
            'eigenstates': 1.0e-4,
        },
        'normal': {},
        'tight': {
            'energy': 1.0e-6,
            'density': 1.0e-5,
            'eigenstates': 1.0e-8,
        },
        'very-tight': {
            'energy': 1.0e-8,
            'density': 1.0e-6,
            'eigenstates': 1.0e-10,
        },
    }
    eigensolver = {
        'default': None,
        'fast': 'rmm-diis',
        'robust': 'dav',
    }

    return {
        'convergence': convergence[resolved['accuracy']],
        'mixing': resolved['mixing'],
        'maxiter': resolved['max_steps'],
        'eigensolver': eigensolver[resolved['solver']],
    }


def qe_conv_thr_to_accuracy(conv_thr):
    """Map a native QE threshold to the nearest portable accuracy profile."""
    conv_thr = float(conv_thr)

    if conv_thr <= 1.0e-10:
        return 'very-tight'
    if conv_thr <= 1.0e-8:
        return 'tight'
    if conv_thr <= 1.0e-6:
        return 'normal'
    return 'loose'
