#[[ Custom find script for JsonCpp as
    CMake doesn't provide one yet.

    Creates an imported interface target called JsonCpp that
    can be added to the package dependencies of an MITK module or plugin. ]]

find_path(JsonCpp_INCLUDE_DIR
  NAMES jsoncpp/jsoncpp
  PATHS ${MITK_EXTERNAL_PROJECT_PREFIX}
  PATH_SUFFIXES include
)

if(NOT TARGET JsonCpp)
  add_library(JsonCpp INTERFACE IMPORTED GLOBAL)
  target_include_directories(JsonCpp INTERFACE ${JsonCpp_INCLUDE_DIR})
  # if(CMAKE_CXX_COMPILER_ID MATCHES Clang OR CMAKE_CXX_COMPILER_ID STREQUAL GNU)
  #   target_compile_options(JsonCpp INTERFACE -fno-strict-aliasing)
  # endif()
endif()

find_package_handle_standard_args(JsonCpp
  FOUND_VAR JsonCpp_FOUND
  REQUIRED_VARS JsonCpp_INCLUDE_DIR
)
