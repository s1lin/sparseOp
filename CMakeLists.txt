cmake_minimum_required(VERSION 3.5)
project(SparseOp2 CXX)
set(CMAKE_AUTOMOC ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

endif ()

INCLUDE_DIRECTORIES(/usr/local/include/eigen3/)
# Add sub directories
add_subdirectory(helper)
add_subdirectory(dataStructure)
add_subdirectory(solver)
add_subdirectory(bin)