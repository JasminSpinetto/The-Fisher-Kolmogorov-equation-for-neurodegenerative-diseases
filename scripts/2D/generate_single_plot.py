from paraview.simple import *
import os
import sys

# --- Configuration ---
IMAGE_RESOLUTION = [1200, 1000]
# The TOTAL_TIME of the simulation is needed to set a CONSISTENT color scale
TOTAL_TIME = 200

# --- Input Validation ---
if len(sys.argv) < 3:
    print("Usage: pvpython generate_single_plot.py <path_to_vtu> <output_png>")
    sys.exit(1)
vtu_file_path = sys.argv[1]
output_png_path = sys.argv[2]

if not os.path.exists(vtu_file_path):
    print(f"Error: Input file not found at {vtu_file_path}")
    sys.exit(1)

# --- 1. Load the pre-computed data ---
reader = OpenDataFile(vtu_file_path)

# --- 2. Visualization Pipeline ---
view = CreateView('RenderView')
view.ViewSize = IMAGE_RESOLUTION
view.Background = [0.32, 0.34, 0.43] # ParaView's default gray background

# --- THE DEFINITIVE FIX: Hide all extra annotations ---
view.OrientationAxesVisibility = 0 # Hide the little XYZ axes

display = Show(reader, view)
ColorBy(display, ('POINTS', 'activationTime'))
lut = GetColorTransferFunction('activationTime')
lut.ApplyPreset('Turbo', True)

# --- THE DEFINITIVE FIX: Set a CONSISTENT data range for all plots ---
# This ensures that a specific color means the same time in every image.
lut.RescaleTransferFunction(0, TOTAL_TIME)

# This command gets a reference to the color bar
scalar_bar = GetScalarBar(lut, view)
# --- THE DEFINITIVE FIX: Make the color bar invisible ---
scalar_bar.Visibility = 0

# Set final representation properties for a clean look
display.Representation = 'Surface'
display.ColorArrayName = ['POINTS', 'activationTime']
display.LookupTable = lut
display.Interpolation = 'Gouraud' 

# --- 3. Frame the Shot and Save ---
ResetCamera()
view.CameraParallelScale *= 0.9 # Zoom in slightly

SaveScreenshot(output_png_path, view, ImageResolution=IMAGE_RESOLUTION)
print(f"Successfully saved final plot to {output_png_path}")