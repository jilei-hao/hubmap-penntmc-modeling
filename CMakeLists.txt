cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)

PROJECT(hubmap-model-generator)

find_package(VTK COMPONENTS 
  vtkCommonColor
  vtkCommonCore
  vtkCommonDataModel
  vtkCommonTransforms
  vtkFiltersCore
  vtkFiltersGeneral
  vtkFiltersSources
  vtkImagingStencil
  vtkIOLegacy
  vtkIOXML
  QUIET
)

if (NOT VTK_FOUND)
  message("Stopping Configuration: ${VTK_NOT_FOUND_MESSAGE}")
  return()
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")

# vtk_module_autoinit is needed
vtk_module_autoinit(
  TARGETS ${PROJECT_NAME}
  MODULES ${VTK_LIBRARIES}
  )

# set binary output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory(Source)




