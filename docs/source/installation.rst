Installation
============

.. _installation:


Quick Installation of Nanoworks to Linux Systems 
------------------------------------------------

For Debian/Ubuntu systems, run the automated installation script:

.. code-block:: console

    $ curl -fsSL https://github.com/sblisesivdin/nanoworks/releases/latest/download/install-all-Debian-based.sh | bash

The release installer pins the Python package to the same Nanoworks version as
the GitHub release. For example, the exact-version URL for v26.8.0 is
``https://github.com/sblisesivdin/nanoworks/releases/download/v26.8.0/install-all-Debian-based.sh``.

Choose the DFT engine and optional components
---------------------------------------------

The development installer offers independent selections for the DFT engine
(QE, GPAW or both) and components (DFT only, DFT + MD, DFT + ML or all).
Defaults are both engines and all components. These selections will be
included in the next release installer; older release assets retain their
original menus. From a checkout of the current development branch:

.. code-block:: console

    $ NANOWORKS_ENGINE=qe NANOWORKS_COMPONENTS=dft bash install_scripts/install-all-Debian-based.sh --dry-run
    $ NANOWORKS_ENGINE=qe NANOWORKS_COMPONENTS=dft bash install_scripts/install-all-Debian-based.sh

``--dry-run`` prints the selected Python extras, system packages and resource
installation steps without creating a virtual environment, invoking sudo,
or downloading packages. Valid environment values are ``qe``, ``gpaw`` or
``both`` for ``NANOWORKS_ENGINE``, and ``dft``, ``md``, ``ml`` or ``all`` for
``NANOWORKS_COMPONENTS``. Invalid selections stop before installation.
``NANOWORKS_VERSION`` optionally pins the package version; the selected version
must provide the requested extras. Without a pin, pip installs the newest
available PyPI package, not the development checkout itself.

QE-only installation skips GPAW and its build configuration. KIM/ASAP3 and
LAMMPS dependencies are added only for MD; ML Python dependencies are added
only for ML. When QE is selected, the installer also installs the PseudoDojo
PBE scalar- and fully-relativistic sets under ``~/.nanoworks/pseudos/qe/``.
Existing GPAW build configuration is preserved. Selecting QE in an existing
virtual environment does not uninstall previously installed GPAW or other
components. Use a fresh environment when you need strict dependency isolation.

QE-only Python installation
---------------------------

For an externally managed QE installation, GPAW and its build libraries are
not required. Create a virtual environment and select only the QE extra:

.. code-block:: console

    $ python3 -m venv ~/.venv_nw_qe
    $ source ~/.venv_nw_qe/bin/activate
    (.venv_nw_qe) $ python -m pip install "nanoworks[qe]"
    (.venv_nw_qe) $ nanoworks --install-qe-pseudos
    (.venv_nw_qe) $ nanoworks --install-examples

To install the development checkout instead of the published package, run
``python -m pip install ".[qe]"`` from the repository directory. Set
``Engine = 'QE'`` in the input, or override it for one run:

.. code-block:: console

    (.venv_nw_qe) $ dftsolve -E QE -i input.py -g geometry.cif --check
    (.venv_nw_qe) $ dftsolve -E QE -i input.py -g geometry.cif

Engine selection during installation chooses dependencies; it does not change
the input's engine or the default calculation backend. QE executes its own
SCFs and postprocessors. Shared ASE/Phonopy utilities do not require GPAW.
QE 7.4.1 must be available separately; elasticity also needs thermo_pw 2.1.0.
Only the executables needed by the selected calculation stages are required.

Detailed Installation of Nanoworks to Linux Systems 
---------------------------------------------------

Installation of system libraries to Debian-based distributions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can also use the same commands on a pure Debian-based Linux system or Windows systems with WSL. If you do not know how to install Linux on Windows 11 with WSL, you can view `this video <https://www.youtube.com/watch?v=zZf4YH4WiZo>`_. On the WSL system, you can use either Debian or Ubuntu. We recommend Ubuntu due to the support provided by Microsoft. For a full GPAW + MD installation, install the following system files.
QE-only users can follow the minimal Python installation above:

.. code-block:: console

   $ sudo apt update && sudo apt upgrade -y
   $ sudo apt install -y python3-venv python3-pip unzip python-is-python3 \
                    python3-dev libopenblas-dev libxc-dev libscalapack-mpi-dev \
                    libfftw3-dev libkim-api-dev openkim-models libkim-api2 pkg-config \
                    task-spooler build-essential lammps

Installation of system libraries to Fedora-based distributions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For a full GPAW + MD installation, install the following system files.
QE-only users can follow the minimal Python installation above:

.. code-block:: console

   $ sudo dnf update
   $ sudo dnf install python3-devel openblas-devel libxc-devel scalapack-openmpi-devel fftw-devel pkgconf

You also must install `kim-api`, `kim-api-devel`, and `openkim-models`. At the time of writing these instructions, packages for Fedora 43 and 44 cannot be installed remotely. Therefore, we must download them, then install them with dnf locally. The order is important:

.. code-block:: console

   $ wget https://download.copr.fedorainfracloud.org/results/lecris/cmake-ninja/fedora-rawhide-x86_64/08840866-kim-api/kim-api-2.2.1-11.fc43.x86_64.rpm
   $ wget https://download.copr.fedorainfracloud.org/results/lecris/cmake-ninja/fedora-rawhide-x86_64/08840866-kim-api/kim-api-devel-2.2.1-11.fc43.x86_64.rpm
   $ wget https://download.copr.fedorainfracloud.org/results/lecris/cmake-ninja/fedora-rawhide-x86_64/08841484-openkim-models/openkim-models-2021.01.28-12.fc43.x86_64.rpm
   $ sudo dnf install kim-api-2.2.1-11.fc43.x86_64.rpm
   $ sudo dnf install kim-api-devel-2.2.1-11.fc43.x86_64.rpm
   $ sudo dnf install openkim-models-2021.01.28-12.fc43.x86_64.rpm
   
Also on Fedora-bades systems, we may need to specify g++ as the C and C++ compiler in the system. After GPAW 26.7.0, the GPAW C-extension will now be built as C++ code.

.. code-block:: console

   (.venv_nw) $ export CC=g++
   (.venv_nw) $ export CXX=g++

Creation of GPAW configuration file (GPAW installations only)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When installing the GPAW extra, GPAW needs configuration inputs for its compilation. QE-only installations skip this section. For this, a configuration file must be created before installing Nanoworks. Therefore, creating a config file called `siteconfig.py` file is important. You can use any text editor. Here, we are creating a file with the cat command, writing necessary information inside it, then closing it with the Ctrl-D command (^D).

.. code-block:: console

   $ mkdir -p ~/.gpaw
   $ cat > ~/.gpaw/siteconfig.py
   fftw = True
   scalapack = True
   libraries = ['xc', 'blas', 'fftw3', 'scalapack-openmpi']
   ^D

If you have problems with libraries fftw, scalapack, you can remove them from the `siteconfig.py` file. They are simply optional. However, for better performance, you need these configuration.

Python Virtual Environment Installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Then, if you do not have a Python environment, create one and activate it:

.. code-block:: console

   $ python -m venv ~/.venv_nw
   $ source ~/.venv_nw/bin/activate

OPTIONAL: DFT-D3 Dispersion Correction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the DFT-D3 dispersion correction, you need to install the dftd3 package on your system. Legacy dftd3 just works fine.

.. code-block:: console

   (.venv_nw) $ cd ~/.venv_nw/bin
   (.venv_nw) $ mkdir build_d3
   (.venv_nw) $ cd build_d3
   (.venv_nw) $ wget https://www.chemie.uni-bonn.de/grimme/de/software/dft-d3/dftd3.tgz
   (.venv_nw) $ tar -xzf dftd3.tgz
   (.venv_nw) $ make
   (.venv_nw) $ mv dftd3 ../
   (.venv_nw) $ cd ..
   (.venv_nw) $ rm -rf build_d3

OPTIONAL: Quantum ESPRESSO Backend
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Nanoworks can also use Quantum ESPRESSO as an alternative DFT engine.
Current native QE support covers PBE plane-wave ground-state,
geometry-optimization, DFT+U, spin-resolved DOS/PDOS, band-structure,
projected-band, pseudo-valence electron-density, DFPT phonon, and
``epsilon.x`` RPA optical workflows. Native QE
``HSE06``, ``HSE03``, and ``PBE0`` support ground-state, DOS/PDOS, band,
projected-band, and density calculations. QE hybrid DOS/PDOS uses a
dedicated SCF state, while hybrid band workflows use ``bands.x`` and, for
projected bands, ``projwfc.x``. QE hybrid geometry, elastic, phonon, and
optical workflows are not supported yet. QE DOS calculations currently use
tetrahedron occupations.

Semilocal QE elastic calculations are driven by ``thermo_pw.x``. This
Nanoworks version supports Quantum ESPRESSO 7.4.1 and thermo_pw 2.1.0. If a
different QE version is detected, Nanoworks warns that it has not been
validated but allows the calculation to continue.

Quantum ESPRESSO itself is not installed automatically by the Nanoworks
Python package or the Debian-based Nanoworks installer. This is intentional:
QE and thermo_pw remain externally managed calculation engines, commonly
provided through local builds or HPC environment modules. A working Quantum
ESPRESSO installation with ``pw.x``,
``dos.x``, ``projwfc.x``, ``bands.x``, ``pp.x``, ``ph.x``, ``q2r.x``,
``matdyn.x``, ``epsilon.x``, and, for elastic calculations,
``thermo_pw.x`` available in ``PATH`` is required for all
currently supported workflows.

Nanoworks can install the required PseudoDojo pseudopotential libraries
for the QE backend:

.. code-block:: console

   (.venv_nw) $ nanoworks --install-qe-pseudos

This installs the standard PBE scalar-relativistic and fully-relativistic
UPF pseudopotential sets under:

.. code-block:: text

   ~/.nanoworks/pseudos/qe/pseudodojo/pbe/

Native QE workflows use the scalar-relativistic PBE set by default. Set
``Pseudo_relativistic = 'full'`` in a ``dftsolve`` or ``dftconverge`` input
to select the fully-relativistic set consistently. This resource selection
does not itself enable spin-orbit coupling. Enable supported nonmagnetic
semilocal QE SOC calculations separately with ``SOC_calc=True``.

Quantum ESPRESSO 7.4.1 is the supported version. The supported thermo_pw
version is 2.1.0.

Installation of Nanoworks and Python Modules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The base ``nanoworks`` package installs common utilities. Choose engine extras
explicitly; GPAW is optional:

.. list-table:: Python dependency selections
   :header-rows: 1
   :widths: 35 65

   * - Selection
     - Install command
   * - QE only
     - ``python -m pip install "nanoworks[qe]"``
   * - GPAW only
     - ``python -m pip install "nanoworks[gpaw]"``
   * - Both DFT engines
     - ``python -m pip install "nanoworks[gpaw,qe]"``
   * - QE + MD
     - ``python -m pip install "nanoworks[qe,md]"``
   * - QE + ML
     - ``python -m pip install "nanoworks[qe,ml]"``
   * - QE + MD + ML
     - ``python -m pip install "nanoworks[qe,md,ml]"``
   * - All components including GPAW
     - ``python -m pip install "nanoworks[all]"``

``[dft]`` remains an alias for the GPAW dependency selection. ``[all]`` includes
GPAW and should not be used for a QE-only environment. Quote extras to prevent
shell glob expansion. Add ``--no-cache-dir`` if you need to bypass pip's cache.

The MD extra installs ASAP3 and kimpy. Its LAMMPS backend also requires the
LAMMPS executable; the Debian/Ubuntu installer includes it when MD is selected.
The ML extra installs PyTorch, MACE, CHGNet and SevenNet. Engine extras do not
install Quantum ESPRESSO or thermo_pw themselves.

Installation of Examples
^^^^^^^^^^^^^^^^^^^^^^^^
Nanoworks software has an examples directory containing numerous examples in the DFT, MD and ML topics. You can copy this directory to the easily accessible ~/.nanoworks/examples location with the following command.

.. code-block:: console

   (.venv_nw) $ nanoworks --install-examples

Each example contains `README.md` files. You can run the related example with a single command in each example.
