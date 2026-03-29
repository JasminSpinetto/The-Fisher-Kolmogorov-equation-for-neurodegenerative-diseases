#!/bin/bash

# This script runs the full 16-case sensitivity analysis for Figure 9.
# It uses the corrected mesh file and a robust command structure for output.

echo "--- Starting Full Sensitivity Analysis for Figure 9 ---"
echo "Loading environment modules..."
module load gcc-glibc dealii

# --- Common Parameters and Paths ---
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
PROJECT_ROOT_DIR="$(dirname "$(dirname "$(dirname "$SCRIPT_DIR")")")"

EXECUTABLE="${PROJECT_ROOT_DIR}/build/neuro_disease_2D"
# --- THE FIX: Use the correct mesh file you have specified ---
MESH_FILE="${PROJECT_ROOT_DIR}/meshes/new_slice_generated.msh" 

VISUALIZATION_SCRIPT="${PROJECT_ROOT_DIR}/scripts/2D/visualize_2d_single_case.py"

TOTAL_TIME=48.0; TIME_STEP=0.24; POLY_DEGREE=1; C_0_VALUE=0.2
GRAY_MATTER_THRESH=5.0; DIMENSION=2
BASELINE_DEXT=1.5; BASELINE_DAXN=3.0; BASELINE_ALPHA=0.6

OUTPUT_BASE_DIR="${PROJECT_ROOT_DIR}/results/2D_tests/figure_9"
rm -rf "$OUTPUT_BASE_DIR" # Clean old results
mkdir -p "$OUTPUT_BASE_DIR"

n=$(nproc)

# --- Main Experiment Loop ---
for seeding_type in 0 1 2 3; do
    case $seeding_type in
        0) disease_name="alpha_synuclein" ;; 1) disease_name="amyloid_beta" ;;
        2) disease_name="tau" ;;             3) disease_name="tdp43" ;;
    esac

    for fiber_model_name in "isotropic" "circumferential" "radial" "axon_based"; do
        
        case $fiber_model_name in
            "radial")          fiber_type=0 ;;
            "circumferential") fiber_type=1 ;;
            "axon_based" | "isotropic") fiber_type=2 ;;
        esac
        
        if [ "$fiber_model_name" == "isotropic" ]; then
            d_axn_value=0.0
        else
            d_axn_value=$BASELINE_DAXN
        fi

        output_dir="${OUTPUT_BASE_DIR}/${disease_name}/${fiber_model_name}/"
        mkdir -p "$output_dir"
        
        echo "========================================================="
        echo "RUNNING CASE: ${disease_name} with ${fiber_model_name}"
        echo "========================================================="
        
        # --- THE FIX: Use the proven, robust command with the correct, hardcoded filename ---
        mpirun -n $n "$EXECUTABLE" \
          -m "$MESH_FILE" -D $DIMENSION -T $TOTAL_TIME -t $TIME_STEP \
          -g $POLY_DEGREE -a $BASELINE_ALPHA -x $d_axn_value -e $BASELINE_DEXT \
          -c $C_0_VALUE -H $GRAY_MATTER_THRESH -s $seeding_type -f $fiber_type \
          -d "$output_dir" -o "solution"
          
        # --- Automatic Verification Step ---
        echo "--- Verifying output for ${disease_name}/${fiber_model_name} ---"
        source /opt/miniforge3/bin/activate pv-env
        python3 "$VISUALIZATION_SCRIPT" "$output_dir"
        
        if [ $? -ne 0 ]; then
            echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            echo "ERROR: Visualization failed for the previous case."
            echo "The output files may be corrupt. Aborting the experiment."
            echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            exit 1
        fi
        conda deactivate
    done
done

echo ""
echo "--- Figure 9 Sensitivity Analysis Finished Successfully ---"