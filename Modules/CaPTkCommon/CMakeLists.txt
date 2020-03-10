include_directories(
	src
)

mitk_create_module(
  INCLUDE_DIRS PUBLIC include src
  PACKAGE_DEPENDS ITK Qt5|Core
  DEPENDS MitkMultilabel MitkCore
)

if(MODULE_IS_ENABLED)
  add_subdirectory(test)
endif()

add_subdirectory(cmdapps)