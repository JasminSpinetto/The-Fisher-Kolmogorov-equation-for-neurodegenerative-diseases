#!/bin/bash

# This script runs the sensitivity analysis for Figure 8 and then
# automatically generates the final 2x2 composite plot.

echo "--- Starting Full Analysis for Figure 8 ---"

# --- STAGE 1: C++ SIMULATIONS ---
echo ""
echo "--- STAGE 1: Running C++ Simulations ---"
echo "Loading environment modules..."
module load gcc-glibc dealii

# --- Common Parameters and Paths ---
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
PROJECT_ROOT_DIR="$(dirname "$(dirname "$(dirname "$SCRIPT_DIR")")")"

EXECUTABLE="${PROJECT_ROOT_DIR}/build/neuro_disease_2D"
MESH_FILE="${PROJECT_ROOT_DIR}/meshes/new_slice_generated.msh"
VISUALIZATION_SCRIPT="${PROJECT_ROOT_DIR}/scripts/2D/figure_8/visualize_2d_figure_8.py"

DIMENSION=2; TOTAL_TIME=48.0; TIME_STEP=0.24; POLY_DEGREE=1
C_0_VALUE=0.2; GRAY_MATTER_THRESH=5.0
SEEDING_TYPE=2; FIBER_TYPE=2
BASELINE_DEXT=1.5; BASELINE_DAXN=3.0; BASELINE_ALPHA=0.6

OUTPUT_BASE_DIR="${PROJECT_ROOT_DIR}/results/2D_tests/figure_8"
rm -rf "$OUTPUT_BASE_DIR"; mkdir -p "$OUTPUT_BASE_DIR"

run_simulation() {
    local dext=$1; local daxn=$2; local alpha=$3; local output_sub_dir=$4
    local full_output_dir="${OUTPUT_BASE_DIR}/${output_sub_dir}/"
    mkdir -p "$full_output_dir"

    echo "---------------------------------------------------------"
    echo "RUNNING CASE: ${output_sub_dir}"
    echo "---------------------------------------------------------"

    n=$(nproc)
    mpirun -n $n "$EXECUTABLE" \
        -m "$MESH_FILE" -D $DIMENSION -T $TOTAL_TIME -t $TIME_STEP \
        -g $POLY_DEGREE -e "$dext" -x "$daxn" -a "$alpha" \
        -c $C_0_VALUE -H $GRAY_MATTER_THRESH -s $SEEDING_TYPE -f $FIBER_TYPE \
        -d "$full_output_dir" -o "solution"
}

if [ ! -f "$EXECUTABLE" ]; then
    echo "Error: Executable not found at $EXECUTABLE"; exit 1
fi

run_simulation $BASELINE_DEXT $BASELINE_DAXN $BASELINE_ALPHA "baseline"
run_simulation 6.0 $BASELINE_DAXN $BASELINE_ALPHA "4x_d_ext"
run_simulation $BASELINE_DEXT 24.0 $BASELINE_ALPHA "8x_d_axn"
run_simulation $BASELINE_DEXT $BASELINE_DAXN 1.2 "2x_alpha"

echo "--- C++ Simulations Finished ---"

# --- STAGE 2: PYTHON VISUALIZATION ---
echo ""
echo "--- STAGE 2: Generating Final Plot ---"
source /opt/miniforge3/bin/activate pv-env
python3 "$VISUALIZATION_SCRIPT"

if [ $? -ne 0 ]; then
    echo ""
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo "ERROR: Python visualization script failed."
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    exit 1
fi
conda deactivate

echo ""
echo "--- Full Analysis for Figure 8 Finished Successfully ---"