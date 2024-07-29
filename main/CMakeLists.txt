set (CMAKE_CXX_STANDARD 20)
set (MYAPP hypercutter)


add_executable(${MYAPP}
${CMAKE_SOURCE_DIR}/src/bsp/pmp/surface_mesh.cpp
${CMAKE_SOURCE_DIR}/src/bsp/pmp/io/read_obj.cpp
${CMAKE_SOURCE_DIR}/src/bsp/pmp/io/write_obj.cpp
${CMAKE_SOURCE_DIR}/src/bsp/pmp/algorithms/normals.cpp
${CMAKE_SOURCE_DIR}/src/bsp/mesh.cpp 
${CMAKE_SOURCE_DIR}/src/bsp/mesh.h
${CMAKE_SOURCE_DIR}/src/bsp/bsp.h 
${CMAKE_SOURCE_DIR}/src/bsp/bsp.cpp 
${CMAKE_SOURCE_DIR}/src/bsp/bvh.cpp
${CMAKE_SOURCE_DIR}/src/bsp/util.h 
${CMAKE_SOURCE_DIR}/src/bsp/util.cpp 
${CMAKE_SOURCE_DIR}/src/bsp/oracle.cpp
${CMAKE_SOURCE_DIR}/src/bsp/bsp_util.cpp
${CMAKE_SOURCE_DIR}/src/bsp/tribox3.cpp
main.cpp


)

# pmp use double
add_compile_definitions(PMP_SCALAR_TYPE_64)

# find_package(Catch2 REQUIRED)
#target_link_libraries(${MYAPP} PRIVATE Catch2::Catch2WithMain)




# linux windows flag
# https://stackoverflow.com/questions/9160335/os-specific-instructions-in-cmake-how-to


# message(STATUS "Looking for sub directory - found" ${CMAKE_SOURCE_DIR})