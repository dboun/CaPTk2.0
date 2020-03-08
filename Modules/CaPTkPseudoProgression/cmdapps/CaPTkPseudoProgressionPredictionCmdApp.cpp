/** \file CaPTkPseudoProgressionPredictionCmdApp.cpp
* \brief CLI program for invoking the Glioblastoma PseudoProgression Prediction Index algorithms
* 
* These algorithms can be used to train a new pseudoprogression prediction model
* or to infer a pseudoprogression prediction for subjects.
*/

#include <mitkCommandLineParser.h>
#include <mitkLogMacros.h>

#include "CaPTkPseudoProgressionPredictionAlgorithm.h"


#include <QString>
#include <QStringList>
#include <QCoreApplication>
#include <string>
#include <map>





int main(int argc, char* argv[])
{
    QCoreApplication cmdapp(argc, argv); // needed for QCoreApplication::applicationDirPath() to work anywhere
    // consts placed here for now
    const QString CBICA_MODEL_DIR = QCoreApplication::applicationDirPath() + QString("/models/pseudoprogression_model");

    mitkCommandLineParser parser;

    /**** Set general information about the command line app ****/

    parser.setCategory("CaPTk Cmd App Category");
    parser.setTitle("CaPTk PseudoProgression Prediction Cmd App");
    parser.setContributor("CBICA");
    parser.setDescription(
                "This command line app accepts a directory of subjects and produces either a"
                " newly-trained model or a pseudoprogression prediction result "
                " depending on the provided parameters.");

    parser.setArgumentPrefix("--", "-");

    /**** Add arguments. Unless specified otherwise, each argument is optional.
            See mitkCommandLineParser::addArgument() for more information. ****/


    parser.addArgument(
                "usage",
                "u",
                mitkCommandLineParser::Bool,
                "Usage",
                "Show the usage menu and ignore all other input.",
                us::Any(),
                true); // optional

    parser.addArgument(
                "train",
                "t",
                mitkCommandLineParser::Bool,
                "Enable training",
                "Pass this parameter to train a new model using the input data. By default, testing will be performed.",
                us::Any(),
                true); // optional

    parser.addArgument(
                "input",
                "i",
                mitkCommandLineParser::String,
                "Input subjects",
                "Path to the directory containing all subjects as subdirectories",
                us::Any(),
                false);

    parser.addArgument(
                "output",
                "o",
                mitkCommandLineParser::String,
                "Output Directory",
                "Path to the desired output directory where the model or results will be created",
                us::Any(),
                false);

    parser.addArgument(
                "model",
                "m",
                mitkCommandLineParser::String,
                "Custom Model",
                "Path to the directory containing the model to use. By default, the CBICA CaPTk Model is used.",
                us::Any(),
                true);
    bool parseSuccess;
    std::map<std::string, us::Any> parsedArgs = parser.parseArguments(argc, argv, &parseSuccess);

    if (parser.argumentParsed("u"))
    {
        MITK_INFO << parser.helpText();
        return EXIT_SUCCESS;
    }
    if (!parseSuccess)
    {
        MITK_INFO << "Failed to parse command-line parameters.";
        MITK_INFO << parser.helpText();
        return EXIT_FAILURE;

    }

    bool trainNewModel = false;
    bool useCustomModel = false;

    if (parsedArgs.empty())
    {
        MITK_INFO << "No arguments provided.";
        return EXIT_FAILURE; // Exit, usage information is already displayed
    }

    if (parsedArgs["input"].Empty() || parsedArgs["output"].Empty())
    {
        MITK_INFO << "Missing a required parameter.";
        MITK_INFO << parser.helpText();
        return EXIT_FAILURE;
    }


    // Handle train/predict mode switching
    if (parser.argumentParsed("train"))
    {
        trainNewModel = true;
        MITK_INFO << "operation type: train";
    }
    else
    {
        trainNewModel = false;
        MITK_INFO << "operation type: test";
    }

    MITK_INFO << "subjects: " << parsedArgs["input"].ToString();
    MITK_INFO << "output: " << parsedArgs["output"].ToString();

    if (parser.argumentParsed("model"))
    {
        useCustomModel = true;
        MITK_INFO << "model: " << parsedArgs["model"].ToString();
    }

    else if (parsedArgs["model"].Empty() && !trainNewModel)
    {
        useCustomModel = false;
        MITK_INFO << "Defaulting to the CBICA CaPTk pseudoprogression prediction index model.\n"
                     "This is a model trained on de novo glioblastoma cases.\n"
                     "Please note that this model was created following certain assumptions\n"
                     "(described in the publication below)\n"
                     "It can be used for research purposes only.\n"
                     "H.Akbari, et al. Quantitative radiomics and machine learning to distinguish\n"
                     "true progression from pseudoprogression in patients with GBM, \n"
                     "ASNR 56th Annual Meeting, 2018.";
    }

    // get strings from mandatory arguments
    QString outputDir = QString::fromStdString(
                parsedArgs["output"].ToString()
            );

    QString subjectDir = QString::fromStdString(
                parsedArgs["input"].ToString()
            );


    QString modelDir;
    if (useCustomModel)
    {
        modelDir = QString::fromStdString(parsedArgs["model"].ToString()
                );
    }


    captk::PseudoProgressionPredictionModuleAlgorithm algorithm = captk::PseudoProgressionPredictionModuleAlgorithm();
    auto resAlgorithm = algorithm.Run(
                modelDir,
                subjectDir,
                outputDir,
                trainNewModel,
                useCustomModel,
                CBICA_MODEL_DIR);

    bool ok = std::get<0>(resAlgorithm);
    std::string errorMessage = std::get<1>(resAlgorithm);

    if (ok)
    {
        MITK_INFO << "The CaPTk pseudoprogression algorithm ran successfully.";
    }
    else
    {
        MITK_ERROR << "The CaPTk pseudoprogression algorithm failed while running with the following error:";
        MITK_ERROR << errorMessage;
        return EXIT_FAILURE;
    }
    if (trainNewModel)
    {
        MITK_INFO << "The trained model has been saved in " << outputDir.toStdString();
    }
    else
    {
        MITK_INFO << "PseudoProgression prediction results have been saved in " << outputDir.toStdString();
    }

    return 0;

}
