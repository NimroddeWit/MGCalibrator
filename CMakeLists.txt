# Set the minimum required version of cmake for the project.
cmake_minimum_required(VERSION 3.10)

# Define the project name and optionally include the version and description here
project(Absolutifier VERSION 1.0 DESCRIPTION "Absolute quantification of shotgun sequencing data")

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set binary output directory to bin within the project folder
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

# Add executable target with source files
add_executable(Absolutifier src/main.cpp src/count_fastq_bases.cpp src/read_metadata.cpp)

# Include directories, modify this if there are any header files in a separate include directory
include_directories(include)

# If there are any libraries that need to be linked, find them and link them here
# find_package(SomeLibrary REQUIRED)
# target_link_libraries(Absolutifier PRIVATE ${SomeLibrary_LIBRARIES})

# Set up installation rules. The install command needs to specify the target and destination.
install(TARGETS Absolutifier DESTINATION bin)
