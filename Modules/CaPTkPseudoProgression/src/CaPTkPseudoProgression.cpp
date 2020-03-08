#include "CaPTkPseudoProgression.h"

#include <QtConcurrent/QtConcurrent>
#include <QMessageBox>
#include <QFileInfo>

#include <iostream>
#include <fstream>
#include "mitkLogMacros.h"

#include "CaPTkPseudoProgressionPredictionAlgorithm.h"
namespace captk {
CaPTkPseudoProgression::CaPTkPseudoProgression(
	QObject *parent)
	: QObject(parent)
{
	connect(&m_Watcher, SIGNAL(finished()), this, SLOT(OnAlgorithmFinished()));
    m_CbicaModelDir = QCoreApplication::applicationDirPath() + QString("/models/pseudoprogression_model");
}

void CaPTkPseudoProgression::Run(
        QString modelDir,
        QString subjectDir,
        QString outputDir,
        bool trainNewModel,
        bool useCustomModel
	)
{
    MITK_INFO << "[CaPTkPseudoProgression::Run]";

	/* ---- Check if it's already running ---- */

	if (m_IsRunning)
	{
		QMessageBox msgError;
		msgError.setText(
			"The algorithm is already running!\nPlease wait for it to finish."
		);
		msgError.setIcon(QMessageBox::Critical);
		msgError.setWindowTitle("Please wait");
		msgError.exec();
		return;
	}
	m_IsRunning = true;

	/* ---- Check requirements ---- */

	bool ok = true;              // Becomes false if there is an issue
	std::string problemStr = ""; // Populated if there is an issue



	// TODO: Check requirements (now everything is assumed ok)
    if (subjectDir.isEmpty()) {
        ok = false;
        problemStr += "Please specify a directory containing subject data.\n";
    }
    if (useCustomModel && modelDir.isEmpty()) {
        ok = false;
        problemStr += "Please specify a directory containing a model to use.\n";
    }
    if (outputDir.isEmpty()) {
        ok = false;
        problemStr += "Please specify a valid output directory.\n";
    }

	// Return if there is an issue
	if (!ok)
	{
		QMessageBox msgError;
		msgError.setText(problemStr.c_str());
		msgError.setIcon(QMessageBox::Critical);
		msgError.setWindowTitle("Incorrect state.");
		msgError.exec();
		m_IsRunning = false;
		return;
	}

	/* ---- Run ---- */
    // This should only be reached if all requirements are met.
    // std::bind is used because normally QtConcurrent::run accepts max=5 function arguments
    m_FutureResult = QtConcurrent::run(std::bind(
        &CaPTkPseudoProgression::RunThread, this,
        modelDir,
        subjectDir,
        outputDir,
        trainNewModel,
        useCustomModel,
        m_CbicaModelDir
    ));
    m_Watcher.setFuture(m_FutureResult);
}


void CaPTkPseudoProgression::OnAlgorithmFinished()
{
    MITK_INFO << "[CaPTkPseudoProgression::OnAlgorithmFinished]";

	if (m_FutureResult.result().ok)
	{
		// Execution finished successfully
        QMessageBox msgSuccess;
        QString msg = "A PseudoProgression Prediction Index (SPI) has been calculated for "
                      "the given subjects by applying the specified model. \n\n";
        msg += "SPI index saved in 'results.csv' file in the output directory. \n\n";
        msgSuccess.setText(msg);
        msgSuccess.setIcon(QMessageBox::Information);
        msgSuccess.setWindowTitle("CaPTk PseudoProgression Module Success!");
        msgSuccess.exec();
	}
	else
	{
		// Something went wrong
		QMessageBox msgError;

		msgError.setText(m_FutureResult.result().errorMessage.c_str());
		msgError.setIcon(QMessageBox::Critical);
		msgError.setWindowTitle("CaPTk PseudoProgression Module Error!");
		msgError.exec();
	}

	m_IsRunning = false;
    emit done(); // notify that the module is done all work
}

CaPTkPseudoProgression::Result
CaPTkPseudoProgression::RunThread(
        QString modelDir,
        QString subjectDir,
        QString outputDir,
        bool trainNewModel,
        bool useCustomModel,
        QString cbicaModelDir
	)
{
    MITK_INFO << "[CaPTkPseudoProgression::RunThread]";

    CaPTkPseudoProgression::Result runResult;

    try
    {
        PseudoProgressionPredictionModuleAlgorithm algorithm = PseudoProgressionPredictionModuleAlgorithm();
        // auto resAlgorithm = algorithm.Run(
        //            modelDir,
        //            subjectDir,
        //            outputDir,
        //            trainNewModel,
        //            useCustomModel,
        //            cbicaModelDir
        //             );

        // runResult.ok = std::get<0>(resAlgorithm);
        // runResult.errorMessage = std::get<1>(resAlgorithm);
    }
    catch (const std::exception& e)
    {
        runResult.ok = false;
        runResult.errorMessage = e.what();
    }
    catch(...)
    {
        runResult.ok = false;
        runResult.errorMessage = "Unexpected error!";        
    }

	return runResult;
}


}
