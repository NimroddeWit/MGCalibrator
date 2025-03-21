# Set the minimum required version of CMake
cmake_minimum_required(VERSION 3.10)

# Define the project
project(absolutifier VERSION 1.0 DESCRIPTION "Absolute quantification of shotgun sequencing data")

# Require C++20 (SeqAn3 requires it)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Disable warnings
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    add_compile_options(-w)
endif()

# Set binary output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/..)

# Include directories
include_directories(external/seqan3/include include)

# Add SeqAn3 as a subdirectory (since we manually cloned it)
add_subdirectory(external/seqan3)

# Source files
set(SOURCES
    src/main.cpp
    src/count_bases.cpp
    src/read_metadata.cpp
    src/io_matrix.cpp
)

# Create executable
add_executable(absolutifier ${SOURCES})

# Link with SeqAn3
target_link_libraries(absolutifier PRIVATE seqan3::seqan3)

# Print where the executable will be created
message(STATUS "Executable will be created in: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
