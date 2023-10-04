cmake_minimum_required (VERSION 3.20)

set (project_name "linear_algebra")

set (cmake_cxx_standard 20)

project (${project_name})

add_library (${project_name} linear_algebra.cpp)

install (TARGETS ${project_name} DESTINATION ${CMAKE_INSTALL_PREFIX})

set_target_properties(${project_name}
    PROPERTIES
    DEBUG_POSTFIX "d"
    RELEASE_POSTFIX ""
)