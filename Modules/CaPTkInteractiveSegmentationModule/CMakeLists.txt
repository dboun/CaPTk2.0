mitk_create_module( #CaPTkInteractiveSegmentationModule
  INCLUDE_DIRS
    PUBLIC include
    PRIVATE src
  DEPENDS PUBLIC MitkSegmentationUI MitkCore
  #PACKAGE_DEPENDS PRIVATE GuidelinesSupportLibrary
)

add_subdirectory(cmdapps)