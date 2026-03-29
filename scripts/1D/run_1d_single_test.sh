#!/bin/bash

# This script runs a single simulation with the original, known-working parameters.

echo "--- Running Single Initial Test ---"

# --- THE FIX: Load the necessary environment for running the program ---
# This sets up the LD_LIBRARY_PATH so the OS can find the shared libraries.
echo "Loading environment modules..."
module load gcc-glibc dealii

# --- Define paths and parameters ---
EXECUTABLE="./build/neuro_disease_1D"
# MESH_FILE="meshes/mesh-40.msh" # Using the original [0, 1] mesh
MESH_FILE="meshes/mesh-1D-centered.msh" # [-1,1]

# These are your original parameters from the static config
TOTAL_TIME=20.0
TIME_STEP=0.1
POLY_DEGREE=1
ALPHA=1.8
D_EXT=0.0
D_AXN=0.2
C_0_VALUE=0.4 # The initial peak concentration

# --- Define output location ---
OUTPUT_DIR="single_test_output/"
FILENAME_PREFIX="output"

# Create the output directory
mkdir -p $OUTPUT_DIR

# Check if the executable exists
if [ ! -f "$EXECUTABLE" ]; then
    echo "Error: Executable not found at $EXECUTABLE"
    echo "Please compile the project first by running ./build.sh"
    exit 1
fi

# --- Run the simulation ---
echo "Executing simulation..."
$EXECUTABLE \
  -m $MESH_FILE \
  -T $TOTAL_TIME \
  -t $TIME_STEP \
  -g $POLY_DEGREE \
  -a $ALPHA \
  -x $D_AXN \
  -e $D_EXT \
  -c $C_0_VALUE \
  -d $OUTPUT_DIR \
  -o $FILENAME_PREFIX

echo "--- Simulation Finished ---"
echo "Output files are in the '${OUTPUT_DIR}' directory."
echo "You can open '${OUTPUT_DIR}/${FILENAME_PREFIX}.pvd' in ParaView to see the results."