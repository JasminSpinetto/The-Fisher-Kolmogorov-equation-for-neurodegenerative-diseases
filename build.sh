#!/bin/bash

# Delete build folder if it exists
if [ -d build ]; then
    rm -rf build
fi
# Load neccessary library
module load gcc-glibc dealii

# Create build directory and compile
mkdir build && cd build
cmake ..
# cmake --build .
cmake --build . -j 1


# The C++ code is set to write the mesh to ../meshes/ from the build directory.
# So we check for the file at that location.
# Note the use of "!" for "not", correct straight quotes, and correct path.
if [ ! -f "../meshes/mesh-1D-centered.msh" ]; then
    echo "Mesh file not found. Generating new mesh..."
    ./mesh_1D_generator
else
    echo "Mesh file already exists. Skipping generation."
fi


