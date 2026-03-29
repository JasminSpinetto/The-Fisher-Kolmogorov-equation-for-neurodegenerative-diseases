import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import meshio
import sys

# --- Configuration ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))
RESULTS_BASE_DIR = os.path.join(PROJECT_ROOT, "results/1D_tests/figure_3")
TOTAL_TIME = 20.0
TIME_STEP = 0.1
ACTIVATION_THRESHOLD = 0.2

# --- Main function to process one simulation run ---
def process_one_simulation(directory_path):
    vtu_files = sorted(glob.glob(os.path.join(directory_path, "solution_[0-9][0-9][0-9].0.vtu")))
    if not vtu_files:
        print(f"Warning: No solution files found in {directory_path}")
        return None, None

    mesh = meshio.read(vtu_files[0])
    x_coords = mesh.points[:, 0]
    activation_times = np.full_like(x_coords, np.inf)

    for i, file_path in enumerate(vtu_files):
        current_time = i * TIME_STEP
        mesh = meshio.read(file_path)
        if "u" in mesh.point_data:
            concentration = mesh.point_data["u"]
            newly_activated_mask = (concentration > ACTIVATION_THRESHOLD) & (np.isinf(activation_times))
            activation_times[newly_activated_mask] = current_time

    # Set any points that never activated to the maximum time for consistent coloring
    activation_times[np.isinf(activation_times)] = TOTAL_TIME
            
    return x_coords, activation_times

# --- Main plotting logic ---
if __name__ == "__main__":
    fig, axes = plt.subplots(3, 3, figsize=(12, 10), sharex=True, sharey=True)
    fig.suptitle("Replication of Fisher-Kolmogorov 1D Spreading", fontsize=16)

    alpha_multipliers = [1, 2, 4]
    d_multipliers = [1, 2, 4]
    im = None

    for i, d_mult in enumerate(d_multipliers):
        for j, alpha_mult in enumerate(alpha_multipliers):
            ax = axes[i, j]
            dir_name = f"alpha_{alpha_mult}a_d_{d_mult}d"
            full_path = os.path.join(RESULTS_BASE_DIR, dir_name)
            
            x, t_act = process_one_simulation(full_path)
            
            if x is None: continue

            sort_indices = np.argsort(x)
            x_sorted, t_act_sorted = x[sort_indices], t_act[sort_indices]
            
            # --- THE DEFINITIVE FIX: Create a simple rectangular image ---
            # We tile the 1D activation time data vertically. No masking is used.
            display_grid = np.tile(t_act_sorted, (100, 1))
            
            cmap = plt.get_cmap('jet_r', 12) # Use a discrete colormap
            
            im = ax.imshow(display_grid, extent=[-1, 1, 0, 1], origin='lower', 
                           aspect='auto', cmap=cmap, vmin=0, vmax=TOTAL_TIME)
            
            if i == 0: ax.set_title(f"growth {alpha_mult}xα")
            if j == 2: ax.text(1.1, 0.5, f"spreading {d_mult}xd", rotation=-90,
                               verticalalignment='center', transform=ax.transAxes)

    for ax in axes[-1, :]: ax.set_xlabel("x")
    for ax in axes[:, 0]: ax.set_ylabel("concentration")
    plt.setp(axes, yticks=[0, 1], xticks=[-1, 0, 1])

    fig.subplots_adjust(right=0.85, wspace=0.1, hspace=0.1)
    if im:
        c_ax = fig.add_axes([0.88, 0.15, 0.03, 0.7])
        cbar = fig.colorbar(im, cax=c_ax)
        cbar.set_label("activation time (t)")
    
    output_filename = os.path.join(RESULTS_BASE_DIR, "figure3_replication.png")
    plt.savefig(output_filename, dpi=300)
    print(f"\nPlot saved as {output_filename}")
    
    plt.show()