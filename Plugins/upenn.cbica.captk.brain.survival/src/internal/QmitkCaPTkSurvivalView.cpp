/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center,
Division of Medical and Biological Informatics.
All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#include "QmitkCaPTkSurvivalView.h"

// blueberry
#include <berryConstants.h>
#include <berryIWorkbenchPage.h>

// mitk
#include "mitkApplicationCursor.h"
#include "mitkStatusBar.h"
//#include "mitkPlanePositionManager.h"
#include "mitkPluginActivator.h"

// Qmitk
//#include "QmitkRenderWindow.h"

// us
#include <usGetModuleContext.h>
#include <usModule.h>
#include <usModuleContext.h>
#include <usModuleResource.h>
#include <usModuleResourceStream.h>

// Qt
#include <QDateTime>
#include <QFileDialog>
#include <QInputDialog>
#include <QMessageBox>
#include <QFileDialog>

// CaPTk
#include "CaPTkSurvival.h"

//#include "tinyxml.h"

//#include <itksys/SystemTools.hxx>

//#include <regex>

const std::string QmitkCaPTkSurvivalView::VIEW_ID = "upenn.cbica.captk.views.survival";

QmitkCaPTkSurvivalView::QmitkCaPTkSurvivalView()
  : m_Parent(nullptr)
{
  m_CaPTkSurvival = new captk::CaPTkSurvival(this);
}

QmitkCaPTkSurvivalView::~QmitkCaPTkSurvivalView()
{

}

void QmitkCaPTkSurvivalView::CreateQtPartControl(QWidget *parent)
{
  // setup the basic GUI of this view
  m_Parent = parent;

  m_Controls.setupUi(parent);

  //m_CaPTkSurvival->SetProgressBar(m_Controls.progressBar);

  /**** Connect signals & slots ****/
  // Note: done ahead of time so that initialization can trigger the right events

  connect(m_Controls.m_cbUsageSelector, SIGNAL(currentTextChanged(QString)),
          this, SLOT(OnUsageComboBoxCurrentTextChanged(QString)));

  connect(m_Controls.m_cbModelSourceSelector, SIGNAL(currentTextChanged(QString)),
          this, SLOT(OnModelSourceComboBoxCurrentTextChanged(QString)));

  connect(m_Controls.pushButton_CustomModelDirBrowse, SIGNAL(clicked()),
          this, SLOT(OnCustomModelDirectoryButtonClicked()));

  connect(m_Controls.pushButton_OutputDirBrowse, SIGNAL(clicked()),
          this, SLOT(OnOutputDirectoryButtonClicked()));

  connect(m_Controls.pushButton_SubjectDirBrowse, SIGNAL(clicked()),
          this, SLOT(OnSubjectDirectoryButtonClicked()));

  connect(m_Controls.pushButtonRun, SIGNAL(clicked()),
          this, SLOT(OnRunButtonPressed()));

  /**** Initialize widgets ****/

  // Initialize the usage type combo box
  m_Controls.m_cbUsageSelector->addItems(QStringList() << "Use Existing Model" << "Train New Model");
  // Set combo box to the last user selected value
  m_Controls.m_cbUsageSelector->setCurrentText(
    this->GetPreferences()->Get("SurvivalUsageComboBox", "Use Existing Model")
  );

  // Initialize model source combo box
  m_Controls.m_cbModelSourceSelector->addItems(
    QStringList() << "CBICA CaPTk Model" << "Custom"
  );
  // Set combo box to the last user selected value
  m_Controls.m_cbModelSourceSelector->setCurrentText(
              this->GetPreferences()->Get("SurvivalModelSourceComboBox", "CBICA CaPTk Model")
  );





}

void QmitkCaPTkSurvivalView::Activated()
{

}

void QmitkCaPTkSurvivalView::Deactivated()
{
  // Not yet implemented
}

void QmitkCaPTkSurvivalView::Visible()
{
  // Not yet implemented
}

void QmitkCaPTkSurvivalView::Hidden()
{
  // Not yet implemented
}

int QmitkCaPTkSurvivalView::GetSizeFlags(bool width)
{
  if (!width)
  {
    return berry::Constants::MIN | berry::Constants::MAX | berry::Constants::FILL;
  }
  else
  {
    return 0;
  }
}

int QmitkCaPTkSurvivalView::ComputePreferredSize(bool width,
                                                          int /*availableParallel*/,
                                                          int /*availablePerpendicular*/,
                                                          int preferredResult)
{
  if (width == false)
  {
    return 100;
  }
  else
  {
    return preferredResult;
  }
}

/************************************************************************/
/* protected slots                                                      */
/************************************************************************/
void QmitkCaPTkSurvivalView::OnUsageComboBoxCurrentTextChanged(const QString& text)
{
  // Remember it
  this->GetPreferences()->Put("SurvivalUsageComboBox", text);
  this->GetPreferences()->Flush();

  if (text == "Use Existing Model") {
      m_Controls.m_cbModelSourceSelector->setVisible(true);
      m_Controls.label_Source->setVisible(true);
      if (m_Controls.m_cbModelSourceSelector->currentText() == "Custom") {
          m_Controls.lineEdit_CustomModelDir->setVisible(true);
          m_Controls.label_CustomModelDir->setVisible(true);
          m_Controls.pushButton_CustomModelDirBrowse->setVisible(true);
      }
      else if (m_Controls.m_cbModelSourceSelector->currentText() == "CBICA CaPTk Model") {
          m_Controls.groupBox_PaperInfo->setVisible(true);
      }

  }
  else if (text == "Train New Model") {
      m_Controls.m_cbModelSourceSelector->setVisible(false);
      m_Controls.label_Source->setVisible(false);
      m_Controls.lineEdit_CustomModelDir->setVisible(false);
      m_Controls.label_CustomModelDir->setVisible(false);
      m_Controls.pushButton_CustomModelDirBrowse->setVisible(false);
      m_Controls.groupBox_PaperInfo->setVisible(false);

  }
}

void QmitkCaPTkSurvivalView::OnModelSourceComboBoxCurrentTextChanged(const QString& text)
{
  // Remember it
  this->GetPreferences()->Put("SurvivalSourceComboBox", text);
  this->GetPreferences()->Flush();

  // Show/Hide views below it
  if (text == "CBICA CaPTk Model")
  {
    m_Controls.lineEdit_CustomModelDir->setVisible(false);
    m_Controls.label_CustomModelDir->setVisible(false);
    m_Controls.pushButton_CustomModelDirBrowse->setVisible(false);
    m_Controls.groupBox_PaperInfo->setVisible(true);
    m_Controls.label_PaperInformation->setVisible(true);
  }
  else if (text == "Custom")
  {
    m_Controls.lineEdit_CustomModelDir->setVisible(true);
    m_Controls.pushButton_CustomModelDirBrowse->setVisible(true);
    m_Controls.label_CustomModelDir->setVisible(true);
    m_Controls.groupBox_PaperInfo->setVisible(false);
    m_Controls.label_PaperInformation->setVisible(false);
  }

}


void QmitkCaPTkSurvivalView::OnSubjectDirectoryButtonClicked()
{
  auto dirName = QFileDialog::getExistingDirectory(m_Parent,
    tr("Select subject directory"), this->GetLastFileOpenPath());

  if(dirName.isEmpty() || dirName.isNull()) { return; }

  // Find and save file information
  QFileInfo file(dirName);
  if (!file.isDir()) { return; }
  this->SetLastFileOpenPath(dirName);



  // Set the path to the QLineEdit
  m_Controls.lineEdit_SubjectDir->setText(dirName);
}

void QmitkCaPTkSurvivalView::OnCustomModelDirectoryButtonClicked()
{
    auto dirName = QFileDialog::getExistingDirectory(m_Parent,
                                                     tr("Select model directory"),
                                                     this->GetLastFileOpenPath());

    if(dirName.isEmpty() || dirName.isNull()) { return; }

    // Find and save file information
    QFileInfo file(dirName);
    if (!file.isDir()) { return; }
    this->SetLastFileOpenPath(dirName);


    // Set the path to the QLineEdit
    m_Controls.lineEdit_CustomModelDir->setText(dirName);
}

void QmitkCaPTkSurvivalView::OnOutputDirectoryButtonClicked()
{
  auto dirName = QFileDialog::getExistingDirectory(m_Parent, 
    tr("Select output directory"), this->GetLastFileOpenPath());

  if(dirName.isEmpty() || dirName.isNull()) { return; }

  // Find and save file information
  QFileInfo file(dirName);
  if (!file.isDir()) { return; }
  this->SetLastFileOpenPath(dirName);

  // Set the path to the QLineEdit
  m_Controls.lineEdit_OutputDir->setText(dirName);
}

void QmitkCaPTkSurvivalView::OnRunButtonPressed()
{
  QString modelCsvPath = m_Controls.lineEdit_CustomModelDir->text();
  QString subjectDirPath = m_Controls.lineEdit_SubjectDir->text();
  QString outputDirPath = m_Controls.lineEdit_OutputDir->text();

//  m_CaPTkSurvival->Run(
//    featuresCsvPath,
//    responsesCsvPath,
//    classificationKernelStr,
//    configurationStr,
//    folds,
//    samples,
//    modelDirPath,
//    outputDirPath
//  );
}

/************************************************************************/
/* protected                                                            */
/************************************************************************/
void QmitkCaPTkSurvivalView::OnSelectionChanged(berry::IWorkbenchPart::Pointer, const QList<mitk::DataNode::Pointer> &)
{

}

void QmitkCaPTkSurvivalView::OnPreferencesChanged(const berry::IBerryPreferences*)
{

}

void QmitkCaPTkSurvivalView::NodeAdded(const mitk::DataNode *)
{

}

void QmitkCaPTkSurvivalView::NodeRemoved(const mitk::DataNode *)
{

}

void QmitkCaPTkSurvivalView::SetFocus()
{
}

void QmitkCaPTkSurvivalView::UpdateControls()
{
  // Hide views that are not useful
  // m_Controls.label_PatientImage->setVisible(false);

  this->RequestRenderWindowUpdate(mitk::RenderingManager::REQUEST_UPDATE_ALL);
}

void QmitkCaPTkSurvivalView::InitializeListeners()
{

}

QString QmitkCaPTkSurvivalView::GetLastFileOpenPath()
{
  return this->GetPreferences()->Get("LastFileOpenPath", "");
}

void QmitkCaPTkSurvivalView::SetLastFileOpenPath(const QString &path)
{
  this->GetPreferences()->Put("LastFileOpenPath", path);
  this->GetPreferences()->Flush();
}