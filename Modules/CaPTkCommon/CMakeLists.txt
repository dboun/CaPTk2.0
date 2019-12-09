mitk_create_module(
  INCLUDE_DIRS
    PUBLIC include src/yaml-cpp/include
    PRIVATE src
  #DEPENDS PUBLIC MitkCore
  #PACKAGE_DEPENDS PRIVATE GuidelinesSupportLibrary
)

#add_subdirectory(cmdapps)