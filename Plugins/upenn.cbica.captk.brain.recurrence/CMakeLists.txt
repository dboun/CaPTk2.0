project(upenn_cbica_captk_brain_recurrence)

#include_directories(${CTK_INCLUDE_DIRS})

mitk_create_plugin(
  EXPORT_DIRECTIVE MITK_CAPTK_RECURRENCE
  EXPORTED_INCLUDE_SUFFIXES src
  PACKAGE_DEPENDS ITK OpenCV
  #MODULE_DEPENDS MitkCaPTkRecurrence
)
