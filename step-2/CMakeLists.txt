cmake_minimum_required(VERSION 3.12)

project(step-2 LANGUAGES CXX)

include(${CMAKE_CURRENT_SOURCE_DIR}/../1_common/cmake-common.cmake)

add_executable(run src/step_2.cpp)

deal_ii_setup_target(run)
