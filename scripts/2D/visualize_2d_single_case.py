import os
import subprocess
import sys

# This script orchestrates the visualization for a SINGLE case from the Figure 9 experiment.

if __name__ == "__main__":
    # --- Input Validation ---
    if len(sys.argv) < 2:
        print("Usage: python3 ./scripts/2D/visualize_2d_single_case.py <path_to_result_directory>")
        print("Example: python3 ./scripts/2D/visualize_2d_single_case.py results/2D_tests/single_test_fig9/alpha_synuclein_axon_based")
        sys.exit(1)
        
    # --- THE FIX: The command-line argument is the full path to the results ---
    full_path = sys.argv[1]
    case_name = os.path.basename(os.path.normpath(full_path)) # Get the last part of the path for titles

    # --- Configuration ---
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))
    
    # Paths to the generic pvpython helper scripts
    COMPUTE_SCRIPT_PATH = os.path.join(PROJECT_ROOT, "scripts/2D/compute_activation.py")
    PLOT_SCRIPT_PATH = os.path.join(PROJECT_ROOT, "scripts/2D/generate_single_plot.py")
    
    if not os.path.isdir(full_path):
        print(f"Error: Directory not found at path: {full_path}")
        sys.exit(1)
        
    # Define the input/output filenames for the scripts
    input_vtu_for_plotting = os.path.join(full_path, "activation_time.vtu")
    output_png = os.path.join(full_path, "activation_plot.png")

    # --- Stage 1: Pre-process data to create activation_time.vtu ---
    print(f"\n--- STAGE 1: Pre-processing data for case: {case_name} ---")
    compute_command = f'bash -c "source /opt/miniforge3/bin/activate pv-env; pvpython {COMPUTE_SCRIPT_PATH} {full_path}"'
    print(f"Executing: {compute_command}")
    
    compute_result = subprocess.run(compute_command, shell=True, capture_output=True, text=True)
    if compute_result.returncode != 0:
        print(f"--- PVPYTHON (compute) FAILED for {case_name} ---")
        print("STDERR:", compute_result.stderr)
        sys.exit(1)
    else:
        print(compute_result.stdout)

    # --- Stage 2: Generate the final plot from the pre-processed file ---
    print(f"\n--- STAGE 2: Generating final plot for case: {case_name} ---")
    plot_command = f'bash -c "source /opt/miniforge3/bin/activate pv-env; pvpython {PLOT_SCRIPT_PATH} {input_vtu_for_plotting} {output_png}"'
    print(f"Executing: {plot_command}")
    
    plot_result = subprocess.run(plot_command, shell=True, capture_output=True, text=True)
    if plot_result.returncode != 0:
        print(f"--- PVPYTHON (plot) FAILED for {case_name} ---")
        print("STDERR:", plot_result.stderr)
    else:
        print(plot_result.stdout)
        print(f"\n--- SUCCESS! ---")
        print(f"Final plot saved to: {output_png}")