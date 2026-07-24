# Example: Grimme-D3 correction on MoS2

Ground, DOS and Band calculations of Bulk MoS2 can be done with and without Grimme-D3 correction. You must find and install Grimme-D3 library before continue:

## 1. Go to your virtual environment's executable folder

    cd /home/YOUR_USERNAME/.venv/bin

Here, ".venv" is your virtual python environment's name. Change to yours if necessary.

## 2. Create a temporary build folder to avoid clutter

    mkdir build_d3
    cd build_d3

# 3. Download the original and pure D3 source code from the University of Bonn

    wget https://www.chemie.uni-bonn.de/grimme/de/software/dft-d3/dftd3.tgz

(This link is working at Apr 16th 2026. Find new link, if necessary)

# 4. Extract the archive and compile it (It automatically uses gfortran that we installed earlier)

    tar -xzf dftd3.tgz
    make

# 5. Move the successfully created real 'dftd3' program to the parent folder (bin)

    mv dftd3 ../

# 6. Delete the temporary folder and clean up

    cd ...
    rm -rf build_d3

To run the calculation with MPI on 4 cores for without vdW correction, please execute the following command in this folder.

    dftsolve -p 4 -i MoS2-wo-vdW.py -g Bulk_MoS2.cif

Then run with vdW correction:

    dftsolve -p 4 -i MoS-with-vdW.py -g Bulk_MoS2.cif

Because we use `Outdirname` variable, results are saved in different folders.
