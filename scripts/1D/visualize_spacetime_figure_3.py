
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import meshio
import sys

# --- Configuration ---
# --- THE FIX: Define all paths relative to the project root ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))

RESULTS_BASE_DIR = os.path.join(PROJECT_ROOT, "results/1D_tests/figure_3")
TOTAL_TIME = 20.0

# --- Main function to process one simulation run ---
def process_one_simulation(directory_path):
    """
    Finds and reads the sequence of raw .0.vtu files and returns
    the spatial coordinates and the stacked space-time data.
    """
    vtu_files = sorted(glob.glob(os.path.join(directory_path, "solution_[0-9][0-9][0-9].0.vtu")))
    
    if not vtu_files:
        print(f"Warning: No solution files found in {directory_path}")
        return None, None
        
    mesh = meshio.read(vtu_files[0])
    x_coords = mesh.points[:, 0]
    
    concentration_data = []

    for file_path in vtu_files:
        mesh = meshio.read(file_path)
        if "u" not in mesh.point_data:
            print(f"Warning: Data 'u' not found in a timestep. Skipping.")
            continue
        
        sort_indices = np.argsort(mesh.points[:, 0])
        sorted_concentration = mesh.point_data["u"][sort_indices]
        concentration_data.append(sorted_concentration)

    x_coords_sorted = np.sort(x_coords)
        
    if not concentration_data:
        print(f"Error: No valid data found in {directory_path}")
        return None, None
    
    spacetime_data = np.vstack(concentration_data)

    return x_coords_sorted, spacetime_data

# --- Main plotting logic ---
if __name__ == "__main__":
    fig, axes = plt.subplots(3, 3, figsize=(12, 10), sharex=True, sharey=True)
    fig.suptitle("Space-Time Evolution of 1D Spreading (Concentration)", fontsize=16)

    alpha_multipliers = [1, 2, 4]
    d_multipliers = [1, 2, 4]
    
    im = None
    for i, d_mult in enumerate(d_multipliers):
        for j, alpha_mult in enumerate(alpha_multipliers):
            ax = axes[i, j]
            dir_name = f"alpha_{alpha_mult}a_d_{d_mult}d"
            full_path = os.path.join(RESULTS_BASE_DIR, dir_name)
            
            x, spacetime_data = process_one_simulation(full_path)
            
            if x is None: continue

            im = ax.imshow(spacetime_data, extent=[-1, 1, 0, TOTAL_TIME], origin='lower', 
                           aspect='auto', cmap='jet', vmin=0, vmax=1.0)
            
            if i == 0:
                ax.set_title(f"growth {alpha_mult}xα")
            if j == 2:
                ax.text(1.1, 0.5, f"spreading {d_mult}xd", rotation=-90,
                        verticalalignment='center', transform=ax.transAxes)

    for ax in axes[-1, :]:
        ax.set_xlabel("space (x)")
    for ax in axes[:, 0]:
        ax.set_ylabel("time (t)")
        
    plt.setp(axes, xticks=[-1, 0, 1])

    fig.subplots_adjust(right=0.85, wspace=0.15, hspace=0.15)
    if im:
        cbar_ax = fig.add_axes([0.88, 0.15, 0.04, 0.7])
        cbar = fig.colorbar(im, cax=cbar_ax)
        cbar.set_label("concentration (u)")
    
    # --- THE FIX: Save the plot in the correct output directory ---
    output_filename = os.path.join(RESULTS_BASE_DIR, "figure3_spacetime_replication.png")
    plt.savefig(output_filename, dpi=300)
    print(f"\nPlot saved as {output_filename}")
    
    # plt.show() # Comment out for non-GUI environments