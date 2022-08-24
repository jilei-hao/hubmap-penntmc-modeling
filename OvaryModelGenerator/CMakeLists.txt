cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

PROJECT(OvaryModelGenerator)
message(${VTK_DIR})
find_package(VTK COMPONENTS 
  CommonColor
  CommonCore
  CommonDataModel
  CommonTransforms
  FiltersCore
  FiltersGeneral
  FiltersGeometry
  FiltersSources
  FiltersPoints
  ImagingStencil
  IOGeometry
  IOLegacy
  IOImage
  IOXML
  IOExport
  InteractionStyle
  RenderingContextOpenGL2
  RenderingCore
  RenderingFreeType
  RenderingGL2PSOpenGL2
  RenderingOpenGL2
  IOExport
  QUIET
)

# include(FindLibArchive)
# if (NOT LibArchive_FOUND)
#   message(FATAL_ERROR "Stopping Configuration: LibArchive not found")
# endif()


if (NOT VTK_FOUND)
  message(FATAL_ERROR "Stopping Configuration: ${VTK_NOT_FOUND_MESSAGE}")
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/VH_F_Pelvis.glb ${CMAKE_CURRENT_BINARY_DIR}/data/VH_F_Pelvis.glb COPYONLY)

add_library(TestingHelper TestingHelper.cxx)
target_link_libraries(TestingHelper PRIVATE ${VTK_LIBRARIES})


# target_link_libraries(Main PRIVATE ${VTK_LIBRARIES})
add_library(OvaryModelGenerator OvaryModelGenerator.cxx)
add_executable(ModelGenerator Main.cxx PelvisModelGenerator.cxx)
target_link_libraries(ModelGenerator PRIVATE ${VTK_LIBRARIES} OvaryModelGenerator)
target_link_libraries(OvaryModelGenerator PRIVATE ${VTK_LIBRARIES} TestingHelper)
target_include_directories(OvaryModelGenerator PRIVATE)


# vtk_module_autoinit is needed
vtk_module_autoinit(
  TARGETS OvaryModelGenerator
  MODULES ${VTK_LIBRARIES}
  )
