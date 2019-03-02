cmake_minimum_required(VERSION 2.8.9)

project(imageProcessing)

include_directories(include)

set(SOURCES
        src/main.cpp
        src/Task2.cpp
        src/Task3
        src/Fourier include/Histogram.h include/ImageUtils.h)

add_executable(imageProcessing ${SOURCES})

target_link_libraries(imageProcessing -lX11 -pthread)