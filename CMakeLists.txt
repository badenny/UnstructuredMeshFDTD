cmake_minimum_required(VERSION 3.22)
project(UnstructuredMeshFDTD)

set(CMAKE_CXX_STANDARD 20)
#set(CMAKE_CXX_COMPILER "clang++")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3 -v -stdlib=libc++ -std=c++20")
set(CMAKE_CXX_COMPILER "/opt/homebrew/Cellar/gcc/12.2.0/bin/aarch64-apple-darwin21-g++-12")
set(CMAKE_C_COMPILER "/opt/homebrew/Cellar/gcc/12.2.0/bin/aarch64-apple-darwin21-gcc-12")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -v -stdlib=libstdc++ -std=c++20 -Wl,-ld_classic -fopenmp")
set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES "/opt/homebrew/Cellar/gcc/12.2.0/include/c++/12")

add_executable(UnstructuredMeshFDTD main.cpp EdgeHash.h)

#Set-up Gmsh
#set(GMSH_DIR /Users/dennyba/third-party/builds/gmsh/install)
set(GMSH_DIR /Users/dennyba/third-party/builds/gmsh-gcc/install)
include_directories(${GMSH_DIR}/include)
link_directories(${GMSH_DIR}/lib)
target_link_libraries(UnstructuredMeshFDTD -L${GMSH_DIR}/lib)
set(GMSH_LINK_FLAGS "-lgmsh")
set(CMAKE_EXE_LINKER_FLAGS  ${CMAKE_EXE_LINKER_FLAGS} ${GMSH_LINK_FLAGS})

#Setup Eigen
set(EIGEN_DIR /Users/dennyba/third-party/source/eigen)
include_directories(${EIGEN_DIR})

