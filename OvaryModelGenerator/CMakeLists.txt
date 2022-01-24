cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

PROJECT(OvaryModelGenerator)

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
  QUIET
)

if (NOT VTK_FOUND)
  message("Stopping Configuration: ${VTK_NOT_FOUND_MESSAGE}")
  return()
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")

add_executable(OvaryModelGenerator OvaryModelGenerator.cxx)
target_link_libraries(OvaryModelGenerator PRIVATE ${VTK_LIBRARIES})
# vtk_module_autoinit is needed
vtk_module_autoinit(
  TARGETS OvaryModelGenerator
  MODULES ${VTK_LIBRARIES}
  )
