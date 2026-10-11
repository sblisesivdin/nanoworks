#!/usr/bin/env bash
# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

# Queue the DFT examples; MD/ML use their own example README commands.
set -euo pipefail
CORENUMBER=4
SCRIPTPATH=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
command -v tsp >/dev/null || { echo "task-spooler (tsp) is required." >&2; exit 1; }
RUNROOT=$(mktemp -d "$SCRIPTPATH/example-runs-XXXXXXXX")
printf 'JobID(1) Engine(1) Case(1)\n' > "$RUNROOT/jobs.dat"
printf 'Nanoworks DFT example run\nProcesses per calculation: %s\nSource: %s\n' \
    "$CORENUMBER" "$SCRIPTPATH" > "$RUNROOT/run-info.txt"
nanoworks --version >> "$RUNROOT/run-info.txt" 2>&1 || true
if command -v git >/dev/null && git -C "$SCRIPTPATH" rev-parse HEAD >/dev/null 2>&1; then
    printf 'Example checkout commit: ' >> "$RUNROOT/run-info.txt"
    git -C "$SCRIPTPATH" rev-parse HEAD >> "$RUNROOT/run-info.txt"
fi

# One queue job owns an entire case, including dependent optical stages.
# A failed stage stops that case; independent queued cases remain runnable.
prepare_case() {
    local engine=$1 folder=$2 name=$3 source
    CASEDIR="$RUNROOT/$engine/$name"
    mkdir -p -- "$CASEDIR"
    for source in "$SCRIPTPATH/$folder/"*.py "$SCRIPTPATH/$folder/"*.cif \
                  "$SCRIPTPATH/$folder/README.md"; do
        if [[ -f "$source" ]]; then cp -- "$source" "$CASEDIR/"; fi
    done
    printf 'QUEUED\n' > "$CASEDIR/status.txt"
    cat > "$CASEDIR/run-case.sh" <<'RUNNER'
#!/usr/bin/env bash
set -uo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")" || exit 1
printf 'RUNNING\n' > status.txt
run_stage() {
    local stage=$1 code
    shift
    printf '%q ' "$@" > "$stage-command.sh"
    printf '\n' >> "$stage-command.sh"
    "$@" > "$stage.log" 2>&1
    code=$?
    printf '%s\n' "$code" > "$stage-exit-code.txt"
    if (( code != 0 )); then
        printf 'FAILED %s %s\n' "$stage" "$code" > status.txt
        exit "$code"
    fi
}
RUNNER
}

add_stage() {
    local stage=$1
    shift
    printf 'run_stage %q ' "$stage" >> "$CASEDIR/run-case.sh"
    printf '%q ' "$@" >> "$CASEDIR/run-case.sh"
    printf '\n' >> "$CASEDIR/run-case.sh"
}

add_dft() {
    local stage=$1 engine=$2 input=$3 geometry=$4
    shift 4
    local args=(dftsolve -E "$engine" -p "$CORENUMBER" -i "$input")
    [[ -f "$CASEDIR/$input" ]] || { echo "Missing input: $input" >&2; exit 1; }
    if [[ -n "$geometry" ]]; then
        [[ -f "$CASEDIR/$geometry" ]] || { echo "Missing geometry: $geometry" >&2; exit 1; }
        args+=(-g "$geometry")
    fi
    args+=("$@")
    add_stage "$stage-check" "${args[@]}" --check
    add_stage "$stage" "${args[@]}"
}

submit_case() {
    local engine=$1 name=$2 job
    printf "printf 'SUCCESS\\n' > status.txt\n" >> "$CASEDIR/run-case.sh"
    job=$(tsp bash "$CASEDIR/run-case.sh")
    printf '%s %s %s\n' "$job" "$engine" "$name" >> "$RUNROOT/jobs.dat"
    printf 'Queued %s: %s/%s\n' "$job" "$engine" "$name"
}

queue_dft() {
    local engine=$1 folder=$2 name=$3 input=$4 geometry=$5
    shift 5
    prepare_case "$engine" "$folder" "$name"
    add_dft calculation "$engine" "$input" "$geometry" "$@"
    submit_case "$engine" "$name"
}

for engine in GPAW QE; do
    # QE variants use its supported tetrahedron DOS. LDA/LCAO examples use
    # PBE/PW with the managed PBE UPFs: these are workflow tests, not identical
    # physical setups for engine-to-engine comparisons.
    dos=()
    lda=()
    lcao=()
    if [[ "$engine" == QE ]]; then
        dos=(-s DOS_integration=tetrahedron)
        lda=(-s XC_calc=PBE)
        lcao=(-s Mode=PW)
    fi
    queue_dft "$engine" Bulk-GaAs-noCIF bulk-gaas bulk_gaas.py "" "${dos[@]}" "${lda[@]}"
    queue_dft "$engine" Cr2O-spin cr2o-spin Cr2O.py Cr2O_mp-1206821_primitive.cif "${dos[@]}"
    queue_dft "$engine" Graphene-LCAO graphene-pristine graphene.py graphene4x4.cif "${dos[@]}" "${lda[@]}" "${lcao[@]}"
    queue_dft "$engine" Graphene-LCAO graphene-defect graphene.py graphene4x4withdefect.cif "${dos[@]}" "${lda[@]}" "${lcao[@]}"
    queue_dft "$engine" Graphene-charged graphene-neutral graphene-neutral.py graphene4x4withdefect.cif "${dos[@]}" "${lda[@]}" "${lcao[@]}"
    queue_dft "$engine" Graphene-charged graphene-charged graphene-charged.py graphene4x4withdefect.cif "${dos[@]}" "${lda[@]}" "${lcao[@]}"

    prepare_case "$engine" Si-2atoms-optical silicon-optical-stages
    add_dft ground "$engine" Si-Step1-ground_dos_band.py Si_mp-149_primitive_Example.cif "${dos[@]}"
    add_dft rpa "$engine" Si-Step2-optical-RPA.py Si_mp-149_primitive_Example.cif
    if [[ "$engine" == GPAW ]]; then
        add_dft bse "$engine" Si-Step3-optical-BSE.py Si_mp-149_primitive_Example.cif
    fi
    submit_case "$engine" silicon-optical-stages

    queue_dft "$engine" ZnO-with-Hubbard zno-hubbard ZnO_withHubbard.py "" "${dos[@]}"
    queue_dft "$engine" ZnO-with-Hubbard zno-reference ZnO_woHubbard.py "" "${dos[@]}"
    queue_dft "$engine" TiC-elastic-electronic tic-elastic TiC.py TiC_primitive_geooptimized.cif "${dos[@]}"
    queue_dft "$engine" Al-phonon al-phonon Al-phonon.py Al_mp-134_primitive.cif
    queue_dft "$engine" Si-phonon si-phonon Si-phonon.py Si_mp-149_primitive.cif
    queue_dft "$engine" SOC-WSe2-noCIF wse2-no-soc WSe2-wo-SOC.py "" "${dos[@]}"
    if [[ "$engine" == QE ]]; then
        queue_dft QE SOC-WSe2-noCIF wse2-soc WSe2-QE-with-SOC.py ""
    else
        queue_dft GPAW SOC-WSe2-noCIF wse2-soc WSe2-with-SOC.py ""
    fi
    queue_dft "$engine" Bulk-MoS2-vdW mos2-reference MoS2-wo-vdW.py Bulk_MoS2.cif "${dos[@]}"
    queue_dft "$engine" Bulk-MoS2-vdW mos2-d3 MoS2-with-vdW.py Bulk_MoS2.cif "${dos[@]}"
    queue_dft "$engine" Si-with-HSE si-hse Si-with-HSE.py Si_mp-149_primitive.cif "${dos[@]}"

    prepare_case "$engine" Convergence silicon-convergence
    add_stage preflight dftconverge --check -p "$CORENUMBER" -g Si.cif -i "Si-$engine-convergence.py"
    add_stage calculation dftconverge -p "$CORENUMBER" -g Si.cif -i "Si-$engine-convergence.py"
    submit_case "$engine" silicon-convergence

    if [[ "$engine" == QE ]]; then
        queue_dft QE Si-2atoms-optical silicon-combined-rpa Si-QE-Combined-RPA-smoke.py Si_mp-149_primitive_Example.cif
    else
        queue_dft GPAW Si-2atoms-optical silicon-combined-rpa Si-Combined-RPA-smoke.py Si_mp-149_primitive_Example.cif
        queue_dft GPAW Si-2atoms-optical silicon-combined-bse Si-Combined-BSE-smoke.py Si_mp-149_primitive_Example.cif
    fi
done

queue_dft QE Fe-QE-spin-smoke fe-spin Fe-QE-spin-smoke.py ""
queue_dft QE Si-QE-elastic si-elastic Si-QE-elastic-smoke.py ""
queue_dft QE NiO-QE-Hubbard-phonon nio-hubbard-phonon NiO-QE-phonon.py ""
queue_dft QE ZnO-with-Hubbard zno-native-hubbard ZnO_QE_Hubbard.py ""

printf '\nResults and logs: %s\n' "$RUNROOT"
echo "Use tsp to monitor jobs. Upload the complete run folder after the queue finishes."
echo "Optical, hybrid, elastic and phonon cases can require substantial time and memory."
