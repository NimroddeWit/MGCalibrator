#!/bin/bash

# Set the script to exit immediately on error
set -e

# Define the build directory
BUILD_DIR="./build"

# Step 1: Create the build directory
echo "Creating the build directory..."
mkdir -p "$BUILD_DIR"

# Step 2: Run CMake and build the project
echo "Configuring and building the project..."
cd "$BUILD_DIR"
cmake ..                   # Configure the project
make -s                    # Build the project

# Step 3: Optionally install the project
echo "Installation Step: You might need administrative privileges for this step."
sudo make install

# Step 4: Clean up - remove the build directory
echo "Cleaning up the build directory..."
cd ..
rm -rf "$BUILD_DIR"

echo "Build and clean-up completed."

