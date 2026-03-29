
#!/bin/bash

# This script runs a single, specific test case from the Figure 9 sensitivity analysis
# and then automatically generates the visualization plot.

echo "--- Running Single Test from Figure 9 Analysis ---"

# --- Step 1: Load the C++ environment ---
echo "Loading environment modules..."
module load gcc-glibc dealii

# --- Step 2: CONFIGURE YOUR TEST HERE ---
# Choose one disease name and one fiber model name from the options below.

# Options for DISEASE_NAME: "alpha_synuclein", "amyloid_beta", "tau", "tdp43"
DISEASE_NAME="tdp43"

# Options for FIBER_MODEL_NAME: "isotropic", "radial", "circumferential", "axon_based"
FIBER_MODEL_NAME="axon_based"

# --- Step 3: Define paths and common parameters ---
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
PROJECT_ROOT_DIR="$(dirname "$(dirname "$(dirname "$SCRIPT_DIR")")")"

EXECUTABLE="${PROJECT_ROOT_DIR}/build/neuro_disease_2D"
MESH_FILE="${PROJECT_ROOT_DIR}/meshes/new_slice_generated.msh" # Using the specified mesh
VISUALIZATION_SCRIPT="${PROJECT_ROOT_DIR}/scripts/2D/visualize_2d_single_case.py"

DIMENSION=2; TOTAL_TIME=48.0; TIME_STEP=0.24; POLY_DEGREE=1
C_0_VALUE=0.2; GRAY_MATTER_THRESH=4.0
BASELINE_DEXT=1.5; BASELINE_DAXN=3.0; BASELINE_ALPHA=0.6

# --- Step 4: Translate names into numeric codes ---
case $DISEASE_NAME in
    "alpha_synuclein") seeding_type=0 ;; "amyloid_beta") seeding_type=1 ;;
    "tau") seeding_type=2 ;;             "tdp43") seeding_type=3 ;;
    *) echo "Error: Invalid DISEASE_NAME '$DISEASE_NAME'"; exit 1 ;;
esac

case $FIBER_MODEL_NAME in
    "radial") fiber_type=0 ;; "circumferential") fiber_type=1 ;;
    "axon_based" | "isotropic") fiber_type=2 ;;
    *) echo "Error: Invalid FIBER_MODEL_NAME '$FIBER_MODEL_NAME'"; exit 1 ;;
esac

if [ "$FIBER_MODEL_NAME" == "isotropic" ]; then
    d_axn_value=0.0
else
    d_axn_value=$BASELINE_DAXN
fi

# --- Step 5: Define the output location ---
OUTPUT_DIR="${PROJECT_ROOT_DIR}/results/2D_tests/single_test_fig9/${DISEASE_NAME}_${FIBER_MODEL_NAME}/"
rm -rf "$OUTPUT_DIR"; mkdir -p "$OUTPUT_DIR"

if [ ! -f "$EXECUTABLE" ]; then
    echo "Error: Executable not found at $EXECUTABLE"
    exit 1
fi

# --- Step 6: Run the simulation ---
echo "Executing simulation for: ${DISEASE_NAME} with ${FIBER_MODEL_NAME} fibers"
n=$(nproc)
mpirun -n $n "$EXECUTABLE" \
  -m "$MESH_FILE" -D $DIMENSION -T $TOTAL_TIME -t $TIME_STEP \
  -g $POLY_DEGREE -a $BASELINE_ALPHA -x $d_axn_value -e $BASELINE_DEXT \
  -c $C_0_VALUE -H $GRAY_MATTER_THRESH -s $seeding_type -f $fiber_type \
  -d "$OUTPUT_DIR" -o "solution"

# --- Step 7: Automatic Visualization ---
echo ""
echo "--- Automatically generating visualization for this case ---"
source /opt/miniforge3/bin/activate pv-env
python3 "$VISUALIZATION_SCRIPT" "$OUTPUT_DIR"

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