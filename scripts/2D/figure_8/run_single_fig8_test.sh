#!/bin/bash

# This script runs a single, specific test case from the Figure 8 sensitivity analysis
# and then automatically generates the visualization plot.

echo "--- Running Single Test from Figure 8 Analysis ---"

# --- Step 1: Load the C++ environment ---
echo "Loading environment modules..."
module load gcc-glibc dealii

# --- Step 2: CONFIGURE YOUR TEST HERE ---
# Choose one of the four case names: "baseline", "4x_d_ext", "8x_d_axn", "2x_alpha"
CASE_NAME="baseline"

# --- Step 3: Define paths and common parameters ---
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
PROJECT_ROOT_DIR="$(dirname "$(dirname "$(dirname "$SCRIPT_DIR")")")"

EXECUTABLE="${PROJECT_ROOT_DIR}/build/neuro_disease_2D"
MESH_FILE="${PROJECT_ROOT_DIR}/meshes/new_slice_generated.msh"
# --- THE CHANGE: Path to your new, generic visualization script ---
VISUALIZATION_SCRIPT="${PROJECT_ROOT_DIR}/scripts/2D/visualize_2d_single_case.py"

DIMENSION=2; TOTAL_TIME=48.0; TIME_STEP=0.24; POLY_DEGREE=1
C_0_VALUE=0.2; GRAY_MATTER_THRESH=4.0
SEEDING_TYPE=2; FIBER_TYPE=2

# --- Step 4: Automatically set parameters based on the chosen CASE_NAME ---
echo "Configuring parameters for case: $CASE_NAME"
case $CASE_NAME in
    "baseline")
        ALPHA_VALUE=0.6; D_EXT_VALUE=1.5; D_AXN_VALUE=3.0 ;;
    "4x_d_ext")
        ALPHA_VALUE=0.6; D_EXT_VALUE=6.0; D_AXN_VALUE=3.0 ;;
    "8x_d_axn")
        ALPHA_VALUE=0.6; D_EXT_VALUE=1.5; D_AXN_VALUE=24.0 ;;
    "2x_alpha")
        ALPHA_VALUE=1.2; D_EXT_VALUE=1.5; D_AXN_VALUE=3.0 ;;
    *) echo "Error: Invalid CASE_NAME '$CASE_NAME'"; exit 1 ;;
esac

# --- Step 5: Define the output location ---
OUTPUT_DIR="${PROJECT_ROOT_DIR}/results/2D_tests/single_test_fig8/${CASE_NAME}/"
FILENAME_PREFIX="solution"
rm -rf "$OUTPUT_DIR"; mkdir -p "$OUTPUT_DIR"

if [ ! -f "$EXECUTABLE" ]; then
    echo "Error: Executable not found at $EXECUTABLE"
    exit 1
fi

# --- Step 6: Run the simulation ---
echo "Executing simulation..."
n=$(nproc)
mpirun -n $n "$EXECUTABLE" \
  -m "$MESH_FILE" -D $DIMENSION -T $TOTAL_TIME -t $TIME_STEP \
  -g $POLY_DEGREE -a $ALPHA_VALUE -x $D_AXN_VALUE -e $D_EXT_VALUE \
  -c $C_0_VALUE -H $GRAY_MATTER_THRESH -s $SEEDING_TYPE -f $FIBER_TYPE \
  -d "$OUTPUT_DIR" -o "$FILENAME_PREFIX"

# --- Step 7: Automatic Visualization ---
echo ""
echo "--- Automatically generating visualization for ${CASE_NAME} ---"
source /opt/miniforge3/bin/activate pv-env
python3 "$VISUALIZATION_SCRIPT" "$OUTPUT_DIR"

# Check the exit code of the python script
if [ $? -ne 0 ]; then
    echo ""
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo "ERROR: Visualization script failed."
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    exit 1
fi
conda deactivate

echo ""
echo "--- Single Test and Visualization Finished Successfully ---"