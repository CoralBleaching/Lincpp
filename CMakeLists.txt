cmake_minimum_required (VERSION 3.12)

set (project_name "linear_algebra")

set (cmake_cxx_standard 20)

project (${project_name})

add_library (${project_name} linear_algebra.cpp)

set_target_properties(${project_name}
    PROPERTIES
    DEBUG_POSTFIX "d"
    RELEASE_POSTFIX ""
)