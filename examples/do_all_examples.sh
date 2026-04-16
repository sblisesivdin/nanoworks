#!/usr/bin/env bash
echo "Nanoworks: "
echo "Adding all examples to tsp queue. Please use tsp after running this."
CORENUMBER=4
SCRIPT=`realpath $0`
SCRIPTPATH=`dirname $SCRIPT`
# Examples
# Bulk-Al-noCIF -------------------
echo "Adding: Bulk-GaAs-noCIF."
cd ./Bulk-GaAs-noCIF
tsp dftsolve -p $CORENUMBER -i bulk_gaas.py

# Cr2O-spin -------------------
echo "Adding: Cr2O-spin"
cd ../Cr2O-spin
tsp dftsolve -p $CORENUMBER -i Cr2O.py -g Cr2O_mp-1206821_primitive.cif

# Graphene-LCAO -------------
echo "Adding: Graphene-LCAO."
cd ../Graphene-LCAO
echo "Step 1: Pristine graphene."
tsp dftsolve -p $CORENUMBER -i graphene.py -g graphene4x4.cif
echo "Step 2: Graphene with defect."
tsp dftsolve -p $CORENUMBER -i graphene.py -g graphene4x4withdefect.cif

# Graphene-charged -------------
echo "Adding: Graphene-charged."
cd ../Graphene-charged
echo "Step 1: Neutral defective graphene."
tsp dftsolve -p $CORENUMBER -i graphene-neutral.py -g graphene4x4withdefect.cif
echo "Step 2: Charged defective Graphene."
tsp dftsolve -p $CORENUMBER -i graphene-charged.py -g graphene4x4withdefect.cif

# Si-2atoms-optical ----------------
echo "adding: Si-2atoms-optical"
cd ../Si-2atoms-optical
echo "Step 1: Ground, DOS, and Band."
tsp dftsolve -p $CORENUMBER -i Si-Step1-ground_dos_band.py -g Si_mp-149_primitive_Example.cif
echo "Step 2: Optical - RPA."
tsp dftsolve -i Si-Step2-optical-RPA.py -g Si_mp-149_primitive_Example.cif
echo "Step 3: Optical - BSE."
tsp dftsolve -i Si-Step3-optical-BSE.py -g Si_mp-149_primitive_Example.cif

# Wurtzite ZnO with DFT+U
echo "Adding: ZnO with DFT+U."
cd ../ZnO-with-Hubbard
echo "Step 1: Ground, DOS, and Band with DFT+U."
tsp dftsolve -p $CORENUMBER -i ZnO_withHubbard.py
echo "Step 2: Ground, DOS, and Band without DFT+U."
tsp dftsolve -p $CORENUMBER -i ZnO_woHubbard.py

# Rocksalt TiC with Elastic Calculations
echo "Adding: Rocksalt TiC."
cd ../TiC-elastic-electronic
tsp dftsolve -p $CORENUMBER -i TiC.py -g TiC_primitive_geooptimized.cif

# Phonon dispersion of Aluminum
echo "Adding: Phonon dispersion of bulk Aluminum."
cd ../Al-phonon
tsp dftsolve -p $CORENUMBER -i Al-phonon.py -g Al_mp-134_primitive.cif

# SOC WSe2
echo "Adding: SOC calculations of 2D WSe2"
cd ../SOC-WSe2-noCIF
tsp dftsolve -p $CORENUMBER -i WSe2-with-SOC.py
tsp dftsolve -p $CORENUMBER -i WSe2-wo-SOC.py

# vdW MoS2
echo "Adding: Bulk MoS2 with Grimme-D3 correction"
cd ../Bulk-MoS2-vdW
tsp dftsolve -p $CORENUMBER -i MoS2-wo-vdW.py -g Bulk_MoS2.cif
tsp dftsolve -p $CORENUMBER -i MoS2-with-vdW.py -g Bulk_MoS2.cif

# Finish
echo "All calculations except the HSE calculation are added. Due to consuming too much time, please run the HSE example separately."
