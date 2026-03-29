import os
import subprocess
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.colors as mcolors
import numpy as np
import sys

# This script orchestrates the visualization process and assembles the final plot.

if __name__ == "__main__":
    # --- Configuration ---
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..', '..'))

    RESULTS_BASE_DIR = os.path.join(PROJECT_ROOT, "results/2D_tests/figure_8")
    COMPUTE_SCRIPT_PATH = os.path.join(PROJECT_ROOT, "scripts/2D/compute_activation.py")
    PLOT_SCRIPT_PATH = os.path.join(PROJECT_ROOT, "scripts/2D/generate_single_plot.py")
    
    cases = ["baseline", "4x_d_ext", "8x_d_axn", "2x_alpha"]
    
    # --- THE FIX: Define the maximum visualization time for the colorbar ---
    VIZ_MAX_TIME = 200.0

    # --- Stage 1: Generate a single plot for each case using pvpython ---
    print("--- STAGE 1: Generating individual plots with pvpython ---")
    for case_name in cases:
        full_path = os.path.join(RESULTS_BASE_DIR, case_name)
        input_vtu_for_plotting = os.path.join(full_path, "activation_time.vtu")
        output_png = os.path.join(full_path, "activation_plot.png")
        
        if not os.path.isdir(full_path):
            print(f"Warning: Directory not found for case '{case_name}'. Skipping.")
            continue

        # (The pre-processing and plotting execution calls remain unchanged)
        compute_command = f'bash -c "source /opt/miniforge3/bin/activate pv-env; pvpython {COMPUTE_SCRIPT_PATH} {full_path}"'
        subprocess.run(compute_command, shell=True, capture_output=True, text=True)

        plot_command = f'bash -c "source /opt/miniforge3/bin/activate pv-env; pvpython {PLOT_SCRIPT_PATH} {input_vtu_for_plotting} {output_png}"'
        subprocess.run(plot_command, shell=True, capture_output=True, text=True)


    # --- Stage 2: Assemble the final figure from the generated plots ---
    print("\n--- STAGE 2: Assembling final 2x2 figure ---")
    
    fig, axes = plt.subplots(2, 2, figsize=(10, 11))
    fig.suptitle("2D Sensitivity Analysis (Replication of Figure 8)", fontsize=16)
    
    case_map = { "baseline": (0, 0), "4x_d_ext": (0, 1), "8x_d_axn": (1, 0), "2x_alpha": (1, 1) }
    titles = { "baseline": "baseline simulation", "4x_d_ext": "4-fold diffusion d_ext", "8x_d_axn": "8-fold axonal transport d_axn", "2x_alpha": "2-fold growth rate α" }

    for case_name, (r, c) in case_map.items():
        ax = axes[r, c]
        image_path = os.path.join(RESULTS_BASE_DIR, case_name, "activation_plot.png")
        
        if os.path.exists(image_path):
            img = mpimg.imread(image_path)
            ax.imshow(img)
        else:
            ax.text(0.5, 0.5, 'Image not found', ha='center', va='center')
        
        ax.set_title(titles[case_name])
        ax.axis('off')

    # --- THE FIX: Add back the shared colorbar with your customizations ---
    fig.subplots_adjust(right=0.85, wspace=0.05, hspace=0.05)
    
    # 1. Use the 'turbo' colormap
    cmap = plt.get_cmap('turbo')
    # 2. Normalize the color range from 0 to 200
    norm = mcolors.Normalize(vmin=0, vmax=VIZ_MAX_TIME)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    # 3. Create the colorbar axes and the colorbar itself
    c_ax = fig.add_axes([0.88, 0.15, 0.04, 0.7])
    cbar = fig.colorbar(sm, cax=c_ax)
    cbar.set_label("activation time (t)")
    
    # 4. Set the custom t0 and tmax labels
    cbar.set_ticks([0, VIZ_MAX_TIME])
    cbar.set_ticklabels([r'$t_0$', r'$t_{max}$'])

    output_filename = os.path.join(PROJECT_ROOT, "results/2D_tests/figure_8/figure8_2D_replication.png")
    plt.savefig(output_filename, dpi=300)
    print(f"\nFinal composite plot saved as {output_filename}")
    
    # plt.show()