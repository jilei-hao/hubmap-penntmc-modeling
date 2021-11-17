cmake_minimum_required(VERSION 2.8)

PROJECT(OvaryModelGenerator)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_executable(OvaryModelGenerator MACOSX_BUNDLE OvaryModelGenerator.cxx)

if(VTK_LIBRARIES)
  target_link_libraries(OvaryModelGenerator ${VTK_LIBRARIES})
else()
  target_link_libraries(OvaryModelGenerator vtkHybrid vtkWidgets)
endif()
