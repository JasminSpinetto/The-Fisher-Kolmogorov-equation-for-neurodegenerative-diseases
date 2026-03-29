#!/bin/bash

# This script replicates the 1D experiment from Weickenmeier et al. (2019)
# using parameters designed to show a traveling wave of activation.

echo "Starting 1D parameter sweep..."
echo "Loading environment modules..."
module load gcc-glibc dealii

# --- THE FIX: Define all paths relative to the project's root directory ---
# Find the directory where this script is located
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
# The project root is two levels up from the script's directory
PROJECT_ROOT_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"

EXECUTABLE="${PROJECT_ROOT_DIR}/build/neuro_disease_1D"
MESH_FILE="${PROJECT_ROOT_DIR}/meshes/mesh-1D-centered.msh"

# --- THE FIX: Set the new, organized output directory ---
OUTPUT_BASE_DIR="${PROJECT_ROOT_DIR}/results/1D_tests/figure_3"

# --- Simulation Parameters ---
TOTAL_TIME=20.0
TIME_STEP=0.1
POLY_DEGREE=1
BASE_ALPHA=1.0
BASE_D_AXN=0.0001 
C_0_VALUE=0.1

# Clean old results and create the new directory structure
rm -rf "$OUTPUT_BASE_DIR"
mkdir -p "$OUTPUT_BASE_DIR"

# Check if the executable exists before starting
if [ ! -f "$EXECUTABLE" ]; then
    echo "Error: Executable not found at $EXECUTABLE"
    exit 1
fi

# --- Main experiment loop ---
for alpha_multiplier in 1 2 4
do
  for d_multiplier in 1 2 4
  do
    CURRENT_ALPHA=$(awk -v base="$BASE_ALPHA" -v mult="$alpha_multiplier" 'BEGIN{print base * mult}')
    CURRENT_D_AXN=$(awk -v base="$BASE_D_AXN" -v mult="$d_multiplier" 'BEGIN{print base * mult}')
    
    # Subdirectories will be created inside the new OUTPUT_BASE_DIR
    OUTPUT_DIR="${OUTPUT_BASE_DIR}/alpha_${alpha_multiplier}a_d_${d_multiplier}d/"
    FILENAME_PREFIX="solution"
    mkdir -p "$OUTPUT_DIR"

    echo "---------------------------------------------------------"
    echo "RUNNING: alpha = $CURRENT_ALPHA, d_axn = $CURRENT_D_AXN"
    echo "Output will be in: $OUTPUT_DIR"
    echo "---------------------------------------------------------"

    # The command to run the C++ executable
    $EXECUTABLE \
      -m "$MESH_FILE" \
      -D 1 \
      -T $TOTAL_TIME \
      -t $TIME_STEP \
      -g $POLY_DEGREE \
      -a $CURRENT_ALPHA \
      -x $CURRENT_D_AXN \
      -e 0.0 \
      -c $C_0_VALUE \
      -d "$OUTPUT_DIR" \
      -o "$FILENAME_PREFIX"
  done
done

echo "Experiment finished. Results are in the '${OUTPUT_BASE_DIR}' directory."