cmake_minimum_required(VERSION 3.19)
project(MolecularDynamics)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

find_package(Boost REQUIRED COMPONENTS)

include_directories(${Boost_INCLUDE_DIRS})

add_executable(MolecularDynamics main.cpp)



