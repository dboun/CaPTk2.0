option(BUILD_CaPTkPseudoProgressionCmdApp "Build command-line app for CaPTkPseudoProgression" ON)

if(BUILD_CaPTkPseudoProgressionCmdApp)
  mitkFunctionCreateCommandLineApp(
    NAME CaPTkPseudoProgressionCmdApp
    CPP_FILES CaPTkPseudoProgressionCmdApp.cpp
    PACKAGE_DEPENDS ITK Qt5|Core
    DEPENDS MitkCaPTkPseudoProgression
  )
endif()
