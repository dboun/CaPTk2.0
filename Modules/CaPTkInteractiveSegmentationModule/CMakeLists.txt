mitk_create_module(
  INCLUDE_DIRS PUBLIC include third_party/jsoncpp/include #third_party/yaml-cpp
  PACKAGE_DEPENDS ITK Qt5|Core+Widgets
  DEPENDS PUBLIC MitkCore MitkMultilabel #MitkCaPTkCommon
  PACKAGE_DEPENDS PUBLIC YamlCpp
)

add_subdirectory(cmdapps)