cmake_minimum_required (VERSION 3.8)

add_library(projector_lib 
include/projector_lib/Geometries.h 
include/projector_lib/Geometry.h 
include/projector_lib/Line.h 
include/projector_lib/Point.h 
include/projector_lib/Projector.hpp 
include/projector_lib/SumAlgorithm.h
include/projector_lib/MyImg.h
src/Line.cpp
src/MyImg.cpp
src/Projector.cpp
)

target_include_directories(projector_lib PUBLIC include)
