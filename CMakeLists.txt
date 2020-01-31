cmake_minimum_required(VERSION 3.13)
project(prog1)

set(CMAKE_CXX_STANDARD 17)

include_directories(.)
include(Warnings.cmake)

add_executable(prog1
        fourier_motzkin.cpp
        fourier_motzkin.hpp
        ineq.cpp
        ineq.hpp
        main.cpp)
