cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(atomic-strain VERSION 1.0.0 LANGUAGES CXX)

if(NOT TARGET coretoolkit::coretoolkit)
    find_package(coretoolkit REQUIRED)
endif()

volt_add_plugin(atomic-strain
    SOURCES src/atomic_strain_engine.cpp
    MAIN_SOURCE src/plugin.cpp
)
