cmake_minimum_required(VERSION 3.16)
project(Modeling)

set(CMAKE_CXX_STANDARD 17)

add_executable(Modeling main.cpp
        Id.h Id.cpp
        utils.h utils.cpp
        ovf_wrapper.hpp ovf_wrapper.cpp)

#target_link_libraries(Modeling armadillo)

message(STATUS "The release flags are: ${CMAKE_CXX_FLAGS_RELEASE}")

#add_subdirectory(libraries/ovf)
add_subdirectory(${PROJECT_SOURCE_DIR}/libraries/ovf)
set(OVF_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/libraries/ovf/include)
target_include_directories(Modeling PRIVATE ${OVF_INCLUDE_DIRS})
target_link_libraries(Modeling PUBLIC ovf_static armadillo)

#add_subdirectory(libraries/ovf)
#target_include_directories(Modeling PRIVATE libraries/ovf/include)
##target_link_libraries(Modeling PUBLIC )
#
#target_link_libraries(Modeling
#        armadillo
#        ${CMAKE_BINARY_DIR}/libraries/ovf/CMakeFiles/libovf_static.a
#        )

#add_subdirectory(libraries/ovf)
#target_include_directories(ovf PUBLIC libraries/ovf/include)
#target_link_libraries(ovf ${OVF_LIBRARIES_STATIC})