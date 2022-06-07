add_executable(${PROJECT_NAME} main.cxx)

target_link_libraries(${PROJECT_NAME} PRIVATE ${VTK_LIBRARIES})