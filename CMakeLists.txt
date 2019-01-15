cmake_minimum_required(VERSION 3.8.2)
project(Matrix)

set(CMAKE_CXX_STANDARD 17)

add_executable(Matrix timeChecker.cpp Matrix.hpp Complex.cpp Complex.h)