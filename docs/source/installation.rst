Installation
============

.. _installation:

Installation of Linux System Configuration 
------------------------------------------

Debian-based distributions (also Windows 11 with WSL)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can also use the same commands on a pure Debian-based Linux system or Windows systems with WSL. If you do not know how to install Linux on Windows 11 with WSL, you can view `this video <https://www.youtube.com/watch?v=zZf4YH4WiZo>`_. On the WSL system, you can use either Debian or Ubuntu. We recommend Ubuntu due to the support provided by Microsoft. First, install the required system files:

.. code-block:: console

   $ sudo apt update && sudo apt upgrade -y
   $ sudo apt install -y python3-venv python3-pip unzip python-is-python3 \
                    python3-dev libopenblas-dev libxc-dev libscalapack-mpi-dev \
                    libfftw3-dev libkim-api-dev openkim-models libkim-api2 pkg-config \
                    task-spooler

Fedora-based distributions
^^^^^^^^^^^^^^^^^^^^^^^^^^

First, install the required system files:

.. code-block:: console

   $ sudo dnf update
   $ sudo dnf install python3-devel openblas-devel libxc-devel scalapack-openmpi-devel fftw-devel pkgconf

You also must install `kim-api`, `kim-api-devel`, and `openkim-models`. At the time of writing these instructions, packages for Fedora 43 cannot be installed remotely. Therefore, we must download them, then install them with dnf locally. The order is important:

.. code-block:: console

   $ wget https://download.copr.fedorainfracloud.org/results/lecris/cmake-ninja/fedora-rawhide-x86_64/08840866-kim-api/kim-api-2.2.1-11.fc43.x86_64.rpm
   $ wget https://download.copr.fedorainfracloud.org/results/lecris/cmake-ninja/fedora-rawhide-x86_64/08840866-kim-api/kim-api-devel-2.2.1-11.fc43.x86_64.rpm
   $ wget https://download.copr.fedorainfracloud.org/results/lecris/cmake-ninja/fedora-rawhide-x86_64/08841484-openkim-models/openkim-models-2021.01.28-12.fc43.src.rpm
   $ sudo dnf install kim-api-2.2.1-11.fc43.x86_64.rpm
   $ sudo dnf install kim-api-devel-2.2.1-11.fc43.x86_64.rpm
   $ sudo dnf install openkim-models-2021.01.28-12.fc43.src.rpm

Python Virtual Environment Installation
---------------------------------------

Then, if you do not have a Python environment, create one and activate it:

.. code-block:: console

   $ python -m venv ~/.venv_nw
   $ source ~/.venv_nw/bin/activate

Python Libraries Installation
-----------------------------
For `dftsolve`, `mdsolve`, and `mlsolve` commands, we need to install many Python libraries. 

ASE and GPAW for dftsolve
^^^^^^^^^^^^^^^^^^^^^^^^^

After preparing your Linux system and environment, you must have the `ase` and `gpaw` codes on your computer. You can find more information about installing `ASE <https://wiki.fysik.dtu.dk/ase/install.html>`_ and `GPAW <https://wiki.fysik.dtu.dk/gpaw/install.html>`_ from their related sites.

.. code-block:: console

   (.venv_nw) $ pip3 install --upgrade ase
   (.venv_nw) $ pip3 install setuptools_scm
   (.venv_nw) $ pip3 install spglib docutils elastic requests phonopy pyrapl pymongo pandas

Creating a `siteconfig.py` file is important. You can use any text editor. Here, we are creating a file with the cat command, writing necessary information inside it, then closing it with the Ctrl-D command (^D).

.. code-block:: console

   (.venv_nw) $ mkdir -p ~/.gpaw
   (.venv_nw) $ cat > ~/.gpaw/siteconfig.py
   fftw = True
   scalapack = True
   libraries = ['xc', 'blas', 'fftw3', 'scalapack-openmpi']
   ^D

If you have problems with libraries fftw, scalapack, you can remove them from `siteconfig.py` file. They are simply optional. Then continue to install gpaw

.. code-block:: console

    (.venv_nw) $ pip3 install --upgrade gpaw

Use `gpaw info` to see installation information. However, PAW datasets are not installed yet. To install them, first create a directory under `~/.gpaw` and then install PAW datasets.

.. code-block:: console

    (.venv_nw) $ mkdir ~/.gpaw/gpaw-setups
    (.venv_nw) $ gpaw install-data --gpaw ~/.gpaw/gpaw-setups/

Installation of ASAP and KIM for mdsolve
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For quick optimization, we need simple interatomic modeling. For this, we need `ASAP3 <https://wiki.fysik.dtu.dk/asap/>`_ for ASE, then we must use `KIM <https://openkim.org/kim-api/>`_ with `OpenKIM <https://openkim.org/>`_ models and `kimpy <https://github.com/openkim/kimpy>`_ libraries.

.. code-block:: console

    (.venv_nw) $ pip3 install --upgrade asap3
    (.venv_nw) $ pip3 install --upgrade kimpy

Installation of Machine Learning libraries for mlsolve
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For machine learning features, we must install the libraries related with `PyTorch <https://pytorch.org/>`_, `MACE (Multi-Atomic Cluster Expansion) <https://github.com/ACEsuit/mace>`_, `CHGNet (Charge-Informed Graph Neural Network) <https://github.com/CederGroupHub/chgnet>`_, `SevenNet (Scalable Equivariance Enabled Neural Network) <https://github.com/MDIL-SNU/SevenNet>`_:

.. code-block:: console

   (.venv_nw) $ pip3 install torch mace-torch chgnet sevenn

Installation of Nanoworks
-------------------------
Then, lastly, install nanoworks:

.. code-block:: console
   
   (.venv_nw) $ pip3 install nanoworks
