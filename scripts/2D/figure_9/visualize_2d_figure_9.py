import os
import subprocess
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.colors as mcolors
import numpy as np
import sys

# This script orchestrates the visualization process and assembles the final plot for Figure 9.

if __name__ == "__main__":
    # --- Configuration ---
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..', '..'))

    RESULTS_BASE_DIR = os.path.join(PROJECT_ROOT, "results/2D_tests/figure_9")
    COMPUTE_SCRIPT_PATH = os.path.join(PROJECT_ROOT, "scripts/2D/compute_activation.py")
    PLOT_SCRIPT_PATH = os.path.join(PROJECT_ROOT, "scripts/2D/generate_single_plot.py")
    
    diseases = ["alpha_synuclein", "amyloid_beta", "tau", "tdp43"]
    fiber_models = ["isotropic", "circumferential", "radial", "axon_based"]
    
    # Define the maximum visualization time for the colorbar
    VIZ_MAX_TIME = 200.0

    # --- Stage 1 & 2: Pre-process data and generate a single plot for each case ---
    print("--- STAGE 1 & 2: Generating individual plots for all 16 cases ---")
    # (The logic for this stage is correct and remains unchanged)
    for disease_name in diseases:
        for fiber_model_name in fiber_models:
            case_name = f"{disease_name} with {fiber_model_name}"
            full_path = os.path.join(RESULTS_BASE_DIR, disease_name, fiber_model_name)
            input_vtu_for_plotting = os.path.join(full_path, "activation_time.vtu")
            output_png = os.path.join(full_path, "activation_plot.png")
            
            if not os.path.isdir(full_path):
                print(f"Warning: Directory not found for case '{case_name}'. Skipping.")
                continue
            if not os.path.exists(input_vtu_for_plotting):
                compute_command = f'bash -c "source /opt/miniforge3/bin/activate pv-env; pvpython {COMPUTE_SCRIPT_PATH} ."'
                subprocess.run(compute_command, shell=True, capture_output=True, text=True, cwd=full_path)
            if not os.path.exists(output_png):
                plot_command = f'bash -c "source /opt/miniforge3/bin/activate pv-env; pvpython {PLOT_SCRIPT_PATH} activation_time.vtu activation_plot.png"'
                subprocess.run(plot_command, shell=True, capture_output=True, text=True, cwd=full_path)

    # --- Stage 3: Assemble the final 4x4 figure ---
    print("\n--- STAGE 3: Assembling final 4x4 figure ---")
    
    fig, axes = plt.subplots(4, 4, figsize=(14, 14))
    fig.suptitle("2D Fiber & Seeding Sensitivity Analysis (Replication of Figure 9)", fontsize=16)
    
    for i, disease_name in enumerate(diseases):
        for j, fiber_model_name in enumerate(fiber_models):
            ax = axes[i, j]
            image_path = os.path.join(RESULTS_BASE_DIR, disease_name, fiber_model_name, "activation_plot.png")
            
            if os.path.exists(image_path):
                img = mpimg.imread(image_path)
                ax.imshow(img)
            else:
                ax.text(0.5, 0.5, 'Image\nnot found', ha='center', va='center', fontsize=9)
            
            ax.axis('off')
            
            if i == 0:
                ax.set_title(fiber_model_name.replace('_', ' ').capitalize(), fontsize=12)
            if j == 0:
                fig.text(0.06, 0.77 - i*0.175, disease_name.replace('_', ' ').capitalize(), 
                         ha='center', va='center', rotation='vertical', fontsize=14)

    # --- THE FIX: Add back the shared colorbar with your customizations ---
    fig.subplots_adjust(left=0.1, right=0.85, wspace=0.05, hspace=0.05)
    
    # 1. Use the 'turbo' colormap
    cmap = plt.get_cmap('turbo')
    # 2. Normalize the color range from 0 to 200
    norm = mcolors.Normalize(vmin=0, vmax=VIZ_MAX_TIME)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    # 3. Create the colorbar axes and the colorbar itself
    c_ax = fig.add_axes([0.88, 0.15, 0.03, 0.7])
    cbar = fig.colorbar(sm, cax=c_ax)
    cbar.set_label("activation time (t)")
    
    # 4. Set the custom t0 and tmax labels
    cbar.set_ticks([0, VIZ_MAX_TIME])
    cbar.set_ticklabels([r'$t_0$', r'$t_{max}$'])

    output_filename = os.path.join(RESULTS_BASE_DIR, "figure9_2D_replication.png")
    plt.savefig(output_filename, dpi=300)
    print(f"\nFinal composite plot saved as {output_filename}")
    
    # plt.show()