cmake_minimum_required(VERSION 3.0)

project(gp)
set(PYBIND11_CPP_STANDARD -std=c++14)
set(CMAKE_CXX_FLAGS)
#set(CMAKE_CXX_FLAGS "-O3 -fopenmp")
add_subdirectory(pybind11)

pybind11_add_module(gp SHARED gp.cpp)  
