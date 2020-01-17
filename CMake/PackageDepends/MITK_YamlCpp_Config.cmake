#[[ This file is empty as everything is already set up in
    CMake/YamlCpp.cmake. However,
    MITK relies on the existence of this file to
    determine if the package was found. ]]

list(APPEND ALL_LIBRARIES ${YamlCpp_LIBRARIES})
if(YamlCpp_INCLUDE_DIRS)
    list(APPEND ALL_INCLUDE_DIRECTORIES ${YamlCpp_INCLUDE_DIRS})
endif()