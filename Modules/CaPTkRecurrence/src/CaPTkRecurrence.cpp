#include "CaPTkRecurrence.h"

#include <QtConcurrent/QtConcurrent>
#include <QMessageBox>
#include <QFileInfo>

#include <iostream>
#include <fstream>
#include "mitkLogMacros.h"

namespace captk 
{

CaPTkRecurrence::CaPTkRecurrence(
	QObject *parent)
	: QObject(parent)
{
	connect(&m_Watcher, SIGNAL(finished()), this, SLOT(OnAlgorithmFinished()));
}

void CaPTkRecurrence::Train(QString subjectDirPath, QString outputModelDirPath)
{
    MITK_INFO << "[CaPTkRecurrence::Train]";

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

	/* ---- Run ---- */

    // std::bind is used because normally QtConcurrent::run accepts max=5 function arguments
    m_FutureResult = QtConcurrent::run(std::bind(
        &CaPTkRecurrence::RunTrainThread, this,
        subjectDirPath,
        outputModelDirPath
    ));
    m_Watcher.setFuture(m_FutureResult);
}

void CaPTkRecurrence::Inference(QString subjectDirPath, QString modelDirPath, QString outputDirPath)
{
    MITK_INFO << "[CaPTkRecurrence::Inference]";

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

	/* ---- Run ---- */

    // std::bind is used because normally QtConcurrent::run accepts max=5 function arguments
    m_FutureResult = QtConcurrent::run(std::bind(
        &CaPTkRecurrence::RunTrainThread, this,
        subjectDirPath,
        modelDirPath,
        outputDirPath
    ));
    m_Watcher.setFuture(m_FutureResult);
}

void CaPTkRecurrence::OnAlgorithmFinished()
{
    MITK_INFO << "[CaPTkRecurrence::OnAlgorithmFinished]";

	if (m_FutureResult.result().ok)
	{
		// Execution finished successfully
        QMessageBox msgSuccess;
        QString msg = "Glioblastoma Infiltration was predicted for "
                      "the given subjects by applying the specified model.";
        msgSuccess.setText(msg);
        msgSuccess.setIcon(QMessageBox::Information);
        msgSuccess.setWindowTitle("CaPTk Recurrence Module Success!");
        msgSuccess.exec();
	}
	else
	{
		// Something went wrong
		QMessageBox msgError;

		msgError.setText(m_FutureResult.result().errorMessage.c_str());
		msgError.setIcon(QMessageBox::Critical);
		msgError.setWindowTitle("CaPTk Recurrence Module Error!");
		msgError.exec();
	}

	m_IsRunning = false;
    emit done(); // notify that the module is done all work
}

CaPTkRecurrence::Result
CaPTkRecurrence::RunTrainThread(QString subjectDirPath, QString outputModelDirPath)
{
    MITK_INFO << "[CaPTkRecurrence::RunTrainThread]";

    CaPTkRecurrence::Result runResult;

    std::string errorMessage;
    bool ok = captk::RecurrenceTrain(
        subjectDirPath.toStdString(), 
        outputModelDirPath.toStdString(),
        errorMessage
    );

	runResult.ok = ok;
	runResult.errorMessage = QString(errorMessage.c_str());

	return runResult;
}

CaPTkRecurrence::Result
CaPTkRecurrence::RunInferenceThread(QString subjectDirPath, QString modelDirPath, QString outputDirPath)
{
    MITK_INFO << "[CaPTkRecurrence::RunInferenceThread]";

    CaPTkRecurrence::Result runResult;

    std::string errorMessage;
    bool ok = captk::RecurrenceInference(
        subjectDirPath.toStdString(), 
        modelDirPath.toStdString(),
        outputDirPath.toStdString(),
        errorMessage
    );

	runResult.ok = ok;
	runResult.errorMessage = QString(errorMessage.c_str());

	return runResult;
}

}
