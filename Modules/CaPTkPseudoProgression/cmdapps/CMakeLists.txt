option(BUILD_CaPTkPseudoProgressionPredictionCmdApp "Build command-line app for CaPTkPseudoProgressionPrediction" ON)

if(BUILD_CaPTkPseudoProgressionPredictionCmdApp)
  mitkFunctionCreateCommandLineApp(
    NAME CaPTkPseudoProgressionPredictionCmdApp
    CPP_FILES CaPTkPseudoProgressionPredictionCmdApp.cpp
    PACKAGE_DEPENDS Qt5|Core+Widgets
    DEPENDS MitkCaPTkPseudoProgression
  )
endif()
