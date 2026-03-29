

#!/bin/bash

# This script runs a comprehensive benchmark of different linear solver preconditioners.
# It captures detailed performance metrics including wall time, Newton iterations,
# and linear solver iterations, then generates a summary report both on-screen
# and to a text file.

echo "--- Running Preconditioner Benchmark ---"

# --- Step 1: Load Environment & Define Paths ---
echo "Loading environment modules..."
module load gcc-glibc dealii

# --- Ensure 'bc' is installed for timing calculations ---
if ! command -v bc &> /dev/null
then
    echo "ERROR: The 'bc' command-line calculator is not installed."
    echo "Please install it (e.g., 'sudo apt-get install bc') and try again."
    exit 1
fi

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
PROJECT_ROOT_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"

EXECUTABLE="${PROJECT_ROOT_DIR}/build/neuro_disease_2D"
MESH_FILE="${PROJECT_ROOT_DIR}/meshes/new_slice_generated.msh"
RESULTS_BASE_DIR="${PROJECT_ROOT_DIR}/results/2D_tests/preconditioner_benchmark"
SUMMARY_FILE="${RESULTS_BASE_DIR}/benchmark_summary.txt"

# --- Step 2: Define Simulation Parameters ---
DIMENSION=2; TOTAL_TIME=48.0; TIME_STEP=0.24; POLY_DEGREE=1
ALPHA_VALUE=0.6; D_EXT_VALUE=1.5; D_AXN_VALUE=3.0
C_0_VALUE=0.2; GRAY_MATTER_THRESH=4.0
SEEDING_TYPE=2; FIBER_TYPE=2
PRECONDITIONERS=("none" "jacobi" "ilu" "amg")

# --- Step 3: Verify Executable and Prepare Output Directory ---
if [ ! -f "$EXECUTABLE" ]; then
    echo "ERROR: Executable not found at $EXECUTABLE"
    exit 1
fi
rm -rf "$RESULTS_BASE_DIR"; mkdir -p "$RESULTS_BASE_DIR"

# --- Step 4: Create the run_simulation function ---
run_simulation() {
    local precon_name=$1
    local output_dir=$2
    local log_file=$3
    local n_procs=$(nproc)

    echo "Running simulation for preconditioner: '${precon_name}' with ${n_procs} processes..."
    echo "Log file: ${log_file}"

    # Direct execution, redirecting both stdout and stderr (2>&1) to the log file.
    mpirun -n "$n_procs" "$EXECUTABLE" \
        -m "$MESH_FILE" -D "$DIMENSION" -T "$TOTAL_TIME" -t "$TIME_STEP" \
        -g "$POLY_DEGREE" -a "$ALPHA_VALUE" -x "$D_AXN_VALUE" -e "$D_EXT_VALUE" \
        -c "$C_0_VALUE" -H "$GRAY_MATTER_THRESH" -s "$SEEDING_TYPE" -f "$FIBER_TYPE" \
        -d "$output_dir" -o "solution" \
        -p "$precon_name" > "$log_file" 2>&1
}

# --- Step 5: Prepare Storage for Results ---
declare -A status_messages

# --- Step 6: Main Benchmark Loop ---
for PRECON in "${PRECONDITIONERS[@]}"; do
    echo ""
    echo "========================================================"
    echo "--- TESTING PRECONDITIONER: $PRECON ---"
    echo "========================================================"

    OUTPUT_DIR="${RESULTS_BASE_DIR}/${PRECON}/"
    LOG_FILE="${OUTPUT_DIR}/simulation.log"
    mkdir -p "$OUTPUT_DIR"

    # Time the execution of the function
    start_time=$(date +%s.%N)
    run_simulation "$PRECON" "$OUTPUT_DIR" "$LOG_FILE"
    end_time=$(date +%s.%N)
    wall_time=$(echo "$end_time - $start_time" | bc -l)

    # --- Step 7: Analyze the log file ---
    if grep -q "Simulation finished." "$LOG_FILE"; then
        TOT_NEWT_IT=$(grep "Total Newton iterations:" "$LOG_FILE" | awk '{print $4}')
        TOT_LIN_IT=$(grep "Total linear iterations:" "$LOG_FILE" | awk '{print $4}')
        TIME_STEPS=$(grep "Total time steps:" "$LOG_FILE" | awk '{print $4}')
        
        AVG_NEWT_PER_STEP=$(echo "scale=2; $TOT_NEWT_IT / $TIME_STEPS" | bc -l)
        AVG_LIN_PER_STEP=$(echo "scale=2; $TOT_LIN_IT / $TIME_STEPS" | bc -l)

        status_messages[$PRECON]=$(printf "Time: %6.2f | Avg Newton/Step: %s | Avg Linear/Step: %s | Total Linear It: %s" \
            "$wall_time" "$AVG_NEWT_PER_STEP" "$AVG_LIN_PER_STEP" "$TOT_LIN_IT")
        echo "RESULT: SUCCESS"
    else
        echo "!!! SIMULATION FAILED. See last 10 lines of log:"
        tail -n 10 "$LOG_FILE"
        status_messages[$PRECON]=$(printf "FAILED after %6.2f seconds. Check log for errors." "$wall_time")
        echo "RESULT: FAILED"
    fi
done

# --- Step 8: Generate the final summary report (both to screen and file) ---
(
# This parenthesis creates a subshell, so we can use 'tee' to send output
# to both the console and the summary file at the same time.
echo ""
echo ""
echo "======================================================================================================================"
echo "---                                           BENCHMARK SUMMARY                                                  ---"
echo "======================================================================================================================"
printf "%-15s | %-16s | %-22s | %-24s | %-20s\n" "Preconditioner" "Wall Time (s)" "Avg Newton Its/Step" "Avg Linear Its/Step" "Total Linear Its"
echo "----------------------------------------------------------------------------------------------------------------------"

for PRECON in "${PRECONDITIONERS[@]}"; do
    if [[ ${status_messages[$PRECON]} == *"Time:"* ]]; then
        wall_time=$(echo "${status_messages[$PRECON]}" | awk -F'|' '{print $1}' | awk '{print $2}')
        # CORRECTED: Use field $3 instead of $4
        avg_newt=$(echo "${status_messages[$PRECON]}" | awk -F'|' '{print $2}' | awk '{print $3}')
        # CORRECTED: Use field $3 instead of $4
        avg_lin=$(echo "${status_messages[$PRECON]}" | awk -F'|' '{print $3}' | awk '{print $3}')
        tot_lin=$(echo "${status_messages[$PRECON]}" | awk -F'|' '{print $4}' | awk '{print $4}')
        printf "%-15s | %-16s | %-22s | %-24s | %-20s\n" "$PRECON" "$wall_time" "$avg_newt" "$avg_lin" "$tot_lin"
    else
        printf "%-15s | %s\n" "$PRECON" "${status_messages[$PRECON]:-DID NOT RUN}"
    fi
done
echo "======================================================================================================================"
) | tee "${SUMMARY_FILE}"

echo ""
echo "Benchmark finished."
echo "Summary report saved to: ${SUMMARY_FILE}"
echo "Detailed logs are in subdirectories of: ${RESULTS_BASE_DIR}"