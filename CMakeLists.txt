cmake_minimum_required(VERSION 3.16)
project(Modeling)

set(CMAKE_CXX_STANDARD 17)

add_executable(Modeling main.cpp
        Id.cpp Id.h
        utils.h)

target_link_libraries(Modeling armadillo)