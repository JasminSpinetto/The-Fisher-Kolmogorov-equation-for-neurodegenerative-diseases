
# This script must be run with pvpython.
from paraview.simple import *
import numpy as np
import vtk
import os
import sys
import glob
from vtk.util import numpy_support

# --- Configuration ---
TIME_STEP = 0.24
ACTIVATION_THRESHOLD = 0.25 

if len(sys.argv) < 2:
    sys.exit(1)
result_dir = sys.argv[1]

# --- File Loading ---
pvtu_files = sorted(glob.glob(os.path.join(result_dir, "solution_[0-9][0-9][0-9].pvtu")))
if not pvtu_files:
    print(f"Error: No solution*.pvtu files found in {result_dir}")
    sys.exit(1)

reader = OpenDataFile(pvtu_files)
scene = GetAnimationScene()
# --- THE FIX: Use the modern, correct syntax for getting timesteps ---
times = scene.TimeKeeper.TimestepValues

# --- Core Activation Time Calculation ---
print(f"Processing {len(times)} timesteps for {result_dir}...")
reader.UpdatePipeline(times[0])
first_data = servermanager.Fetch(reader)
num_points = first_data.GetNumberOfPoints()
activation = np.full(num_points, np.inf)

for i, t in enumerate(times):
    reader.UpdatePipeline(t)
    data = servermanager.Fetch(reader)
    u_array = data.GetPointData().GetArray("u")
    u_numpy = numpy_support.vtk_to_numpy(u_array)
    newly_activated_mask = (u_numpy >= ACTIVATION_THRESHOLD) & (np.isinf(activation))
    activation[newly_activated_mask] = t

activation[np.isinf(activation)] = times[-1]
print("Activation calculation complete.")

# --- Create and Save the New Summary File ---
final_data_object = servermanager.Fetch(reader)
act_vtk_array = numpy_support.numpy_to_vtk(activation)
act_vtk_array.SetName("activationTime")
final_data_object.GetPointData().AddArray(act_vtk_array)

output_filename = os.path.join(result_dir, "activation_time.vtu")
writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName(output_filename)
writer.SetInputData(final_data_object)
writer.Write()
print(f"Success! Wrote pre-processed data to {output_filename}")