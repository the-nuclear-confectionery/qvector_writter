cmake_minimum_required(VERSION 3.12)
project(qvector_writter)

set(CMAKE_CXX_STANDARD 20)

# Use CONDA_PREFIX as the base prefix
if(DEFINED ENV{CONDA_PREFIX})
    set(PREFIX $ENV{CONDA_PREFIX})
else()
    message(FATAL_ERROR "CONDA_PREFIX environment variable is not set. Please activate a Conda environment.")
endif()

message(STATUS "Using CONDA_PREFIX: ${PREFIX}")

# Locate ROOT
find_package(ROOT 6 REQUIRED)

# Locate HepMC3
find_package(HepMC3 REQUIRED)

# Smash include and lib
set(SMASH_INCLUDE_DIR "${PREFIX}/include/smash-3.1")
set(SMASH_LIB_DIR "${PREFIX}/lib/smash-3.1")

# Cuba (for Smash)
set(CUBA_INCLUDE_DIR "${PREFIX}/include/smash-3.1/cuba")

# YAML
find_package(yaml-cpp REQUIRED)

# Include directories
include_directories(SYSTEM
    ${ROOT_INCLUDE_DIRS}
    ${SMASH_INCLUDE_DIR}
    ${CUBA_INCLUDE_DIR}
    ${YAML_CPP_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Source files
set(SOURCES
    src/main.cpp
    src/qvector_writter.cpp
)

# Create executable
add_executable(qvector_writter.exe ${SOURCES})

# Link libraries
target_link_libraries(qvector_writter.exe
    ${ROOT_LIBRARIES}
    ${HEPMC3_ROOTIO_LIB}
    yaml-cpp
    -Wl,-rpath,${SMASH_LIB_DIR}
    -L${SMASH_LIB_DIR} -lsmash
    -lHepMC3 -lHepMC3rootIO
)

# Endianness detection
include(TestBigEndian)
TEST_BIG_ENDIAN(IS_BIG_ENDIAN)
if(IS_BIG_ENDIAN)
    message(STATUS "Big endian architecture detected.")
    add_definitions("-DBIG_ENDIAN_ARCHITECTURE")
else()
    message(STATUS "Little endian architecture detected.")
    add_definitions("-DLITTLE_ENDIAN_ARCHITECTURE")
endif()
