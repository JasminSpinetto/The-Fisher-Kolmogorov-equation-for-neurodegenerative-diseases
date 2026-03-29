Of course. You are absolutely right to want a final, clean `README.md` that reflects the complete, robust, and automated workflow you have built.

Here is the definitive, updated `README.md`. It is clean, accurate, and provides clear instructions for the entire project as it currently stands.

---


# Neurodegenerative Disease Simulation

This project provides a C++ implementation of the reaction-diffusion model for neurodegenerative diseases presented in Weickenmeier et al. (2019), "A physics-based model explains the prion-like features of neurodegeneration...".

The code is designed to be compiled and run in a Linux-based environment (e.g., a Docker container) with the deal.II finite element library. Visualization is handled by a set of Python 3 scripts, which require a Conda environment with ParaView.

### 1. Project Structure

-   `/src`: Contains all C++ source code.
-   `/include`: Contains all C++ header files.
-   `/meshes`: Contains all input mesh files.
-   `/scripts`: Contains all bash and Python scripts for running experiments and visualizations.
-   `/results`: The default location where simulation and visualization outputs are saved.

### 2. Compiling the Code

All compilation should be done inside your configured Linux/Docker environment.

1.  **Load Dependencies:** Before compiling, ensure your environment is set up.
    ```bash
    module load gcc-glibc dealii
    ```
2.  **Run the Build Script:** From the main project directory (`nmpde-neurodegenerative-diseases`), run the provided build script. It will create a `build` directory, configure the project with CMake, and compile all C++ executables.
    ```bash
    ./build.sh
    ```
    This creates the executables `neuro_disease_1D`, `neuro_disease_2D`, etc., inside the `build` folder.

### 3. Running the 1D Experiment (Figure 3 Replication)

This experiment runs a 3x3 parameter sweep to test the effects of growth (`alpha`) and spreading (`d`) on a 1D domain.

#### **Step 1: Run the C++ Simulations**

From the main project directory, execute the experiment script. This will run 9 simulations and save the raw data.

```bash
./scripts/1D/run_1d_experiment.sh
```
-   **Output Location:** `/results/1D_tests/figure_3/`

#### **Step 2: Visualize the 1D Results**

This project includes two different Python scripts to visualize the 1D results.

*   **To create the "Activation Time" plot (like the paper's Figure 3):**
    ```bash
    python3 ./scripts/1D/visualize_1d_experiment.py
    ```
    -   **Output:** `figure3_replication.png` saved in the results directory.

*   **To create the "Space-Time" plot (raw concentration):**
    ```bash
    python3 ./scripts/1D/visualize_spacetime_figure_3.py
    ```
    -   **Output:** `figure3_spacetime_replication.png` saved in the results directory.

### 4. Running the 2D Experiments (Figure 8 & 9)

These experiments run sensitivity analyses on a 2D brain slice mesh. The visualization process requires ParaView.

#### **Step 1: Install ParaView via Conda (One-Time Setup)**

The visualization scripts rely on `pvpython`. The recommended way to install this in the container is with Conda.

1.  **Install Conda/Miniforge:**
    ```bash
    wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
    bash Miniforge3-Linux-x86_64.sh -b -p /opt/miniforge3
    source /opt/miniforge3/bin/activate
    ```
2.  **Create and Activate a ParaView Environment:**
    ```bash
    conda create -n pv-env -c conda-forge paraview -y
    source /opt/miniforge3/bin/activate pv-env
    ```
3.  **Install Python Libraries:** Inside the active environment, install the necessary libraries for plotting.
    ```bash
    pip install matplotlib meshio numpy
    ```
*Note: You must activate the Conda environment (`source /opt/miniforge3/bin/activate pv-env`) in every new terminal session before running the visualization scripts.*

---

#### **4.1 Figure 8: Parameter Sensitivity Analysis**

This runs a 2x2 analysis on the `baseline`, `4x d_ext`, `8x d_axn`, and `2x alpha` cases and automatically generates the final plot.

1.  **Run the Full Workflow:** From the project root, activate your Conda environment and run the main script:
    ```bash
    source /opt/miniforge3/bin/activate pv-env
    ./scripts/2D/figure_8/ext_axn_alpha_sensitivity_analysis_figure_8.sh
    ```
    -   **What it Does:** This single script runs the 4 C++ simulations and then immediately calls the Python visualization pipeline to produce the final composite image.
    -   **Output Location:** Raw data is saved in `/results/2D_tests/figure_8/`. The final plot, `figure8_2D_replication.png`, is saved in the main project directory.

---

#### **4.2 Figure 9: Disease and Fiber Model Analysis**

This runs a full 4x4 analysis on 4 disease types and 4 fiber models and automatically verifies each run.

1.  **Run the Full Workflow:** From the project root, activate your Conda environment and run the main script:
    ```bash
    source /opt/miniforge3/bin/activate pv-env
    ./scripts/2D/figure_9/fiber_seed_sensitivity_analysis_figure_9.sh
    ```
    -   **What it Does:** This script runs all 16 C++ simulations. After each one, it automatically runs a verification/visualization step.
    -   **Output Location:** Raw data and individual plots (`activation_plot.png`) are saved in subfolders inside `/results/2D_tests/figure_9/`.

2.  **Assemble the Final Figure 9 Plot:** After the main script finishes, you can create the final 4x4 composite image. From the project root, run:
    ```bash
    python3 ./scripts/2D/figure_9/visualize_2d_figure9.py
    ```
    -   **Output:** `figure9_2D_replication.png` saved in the `/results/2D_tests/figure_9/` directory.

---

#### **4.3 Running and Visualizing a Single 2D Case**

To debug or analyze a specific experiment, use the single-case scripts. They run both the simulation and the visualization in one command.

*   **For a Figure 8 Case:**
    1.  Open `./scripts/2D/figure_8/run_single_fig8_test.sh`.
    2.  In **Step 2**, set the `CASE_NAME` variable (e.g., to `"8x_d_axn"`).
    3.  Save and run the script from the project root. It will run the C++ code and then immediately generate the plot.
        ```bash
        ./scripts/2D/figure_8/run_single_fig8_test.sh
        ```

*   **For a Figure 9 Case:**
    1.  Open `./scripts/2D/figure_9/run_single_fig9_test.sh`.
    2.  In **Step 2**, set the `DISEASE_NAME` and `FIBER_MODEL_NAME` variables.
    3.  Save and run the script from the project root.
        ```bash
        ./scripts/2D/figure_9/run_single_fig9_test.sh
        ```

### Environment 
run the following command to be able to see the execution time in the preconditioner test 
apt-get update && apt-get install -y bc

## Credits
This project has been developed together with:
- Daniel Yezid Guarnizo Orjuela
- Anna Paola Izzo
- Veronica Di Gennaro