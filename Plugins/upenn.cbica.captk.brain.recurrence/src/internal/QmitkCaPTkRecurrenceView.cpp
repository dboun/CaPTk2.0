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

#include "QmitkCaPTkRecurrenceView.h"

// blueberry
#include <berryConstants.h>
#include <berryIWorkbenchPage.h>

// mitk
#include "mitkApplicationCursor.h"
#include "mitkStatusBar.h"
//#include "mitkPlanePositionManager.h"
#include "mitkPluginActivator.h"
#include "mitkLogMacros.h"

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
#include <QCoreApplication>

// CaPTk
// #include "CaPTkRecurrence.h"

const std::string QmitkCaPTkRecurrenceView::VIEW_ID = "upenn.cbica.captk.views.recurrence";

QmitkCaPTkRecurrenceView::QmitkCaPTkRecurrenceView()
  : m_Parent(nullptr)
{
  // m_CaPTkRecurrence = new captk::CaPTkRecurrence(this);
}

QmitkCaPTkRecurrenceView::~QmitkCaPTkRecurrenceView()
{

}

void QmitkCaPTkRecurrenceView::CreateQtPartControl(QWidget *parent)
{
  // setup the basic GUI of this view
  m_Parent = parent;

  m_Controls.setupUi(parent);


  /**** Connect signals & slots ****/

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

  connect(m_Controls.lineEdit_CustomModelDir, SIGNAL(textChanged(QString)),
          this, SLOT(OnModelPathLineEditTextChanged(QString)));

  // connect(m_CaPTkRecurrence, SIGNAL(done()),
  //         this, SLOT(OnModuleDone()));

  /**** Initialize widgets ****/

  // Initialize the usage type combo box
  m_Controls.m_cbUsageSelector->addItems(QStringList() << "Use Existing Model" << "Train New Model");
  // Set combo box to the last user selected value
  m_Controls.m_cbUsageSelector->setCurrentText(
    this->GetPreferences()->Get("RecurrenceUsageComboBox", "Use Existing Model")
  );

  // Initialize model source combo box
  m_Controls.m_cbModelSourceSelector->addItems(
    QStringList() << "CBICA CaPTk Model" << "Custom"
  );
  // Set combo box to the last user selected value
  m_Controls.m_cbModelSourceSelector->setCurrentText(
              this->GetPreferences()->Get("RecurrenceModelSourceComboBox", "CBICA CaPTk Model")
  );

  // Set model information
  m_Controls.label_PaperInformation->setWordWrap(true);
  QString modelInfo = ("This is a model trained on de novo glioblastoma cases.<br>"
                       "Please note that this model was created following certain assumptions<br>"
                       "(described in the paper below)<br>"
                       "It can be used for research purposes only.<br><br>"
                       "S. Rathore, et al. Radiomic signature of infiltration in peritumoral edema predicts subsequent recurrence in glioblastoma: implications for personalized radiotherapy planning.<br>"
                       "Journal of medical imaging (Bellingham, Wash.)<br>"
                       "doi: 10.1117/1.JMI.5.2.021219.<br><br>"
                       "<a href=\"https://www.ncbi.nlm.nih.gov/pubmed/29531967\" style=\"color:red\">[LINK]</a>");
  m_Controls.label_PaperInformation->setText(modelInfo);
  m_Controls.label_PaperInformation->setTextFormat(Qt::RichText);
  m_Controls.label_PaperInformation->setTextInteractionFlags(Qt::TextBrowserInteraction);
  m_Controls.label_PaperInformation->setOpenExternalLinks(true);

  this->SelectedModelChanged();
}

void QmitkCaPTkRecurrenceView::Activated()
{

}

void QmitkCaPTkRecurrenceView::Deactivated()
{
  // Not yet implemented
}

void QmitkCaPTkRecurrenceView::Visible()
{
  // Not yet implemented
}

void QmitkCaPTkRecurrenceView::Hidden()
{
  // Not yet implemented
}

int QmitkCaPTkRecurrenceView::GetSizeFlags(bool width)
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

int QmitkCaPTkRecurrenceView::ComputePreferredSize(bool width,
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
void QmitkCaPTkRecurrenceView::OnUsageComboBoxCurrentTextChanged(const QString& text)
{
  // set preference to remember usage choice
  this->GetPreferences()->Put("RecurrenceUsageComboBox", text);
  this->GetPreferences()->Flush();

  if (text == "Use Existing Model") {
      m_Controls.m_cbModelSourceSelector->setVisible(true);
      m_Controls.label_Source->setVisible(true);
      m_Controls.groupBox_ModelRequirements->setVisible(true);
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
      m_Controls.groupBox_ModelRequirements->setVisible(false);

  }
}

void QmitkCaPTkRecurrenceView::OnModelSourceComboBoxCurrentTextChanged(const QString& text)
{
  // Set preference to remember model source choice
  this->GetPreferences()->Put("RecurrenceSourceComboBox", text);
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

  this->SelectedModelChanged();
}


void QmitkCaPTkRecurrenceView::OnSubjectDirectoryButtonClicked()
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

void QmitkCaPTkRecurrenceView::OnCustomModelDirectoryButtonClicked()
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

void QmitkCaPTkRecurrenceView::OnModelPathLineEditTextChanged(QString)
{
  if (m_Controls.m_cbModelSourceSelector->currentText() == "Custom")
  {
    this->SelectedModelChanged();
  }
}

void QmitkCaPTkRecurrenceView::OnOutputDirectoryButtonClicked()
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

void QmitkCaPTkRecurrenceView::OnRunButtonPressed()
{
  QString modelDirPath = m_Controls.lineEdit_CustomModelDir->text();
  QString subjectDirPath = m_Controls.lineEdit_SubjectDir->text();
  QString outputDirPath = m_Controls.lineEdit_OutputDir->text();
  bool trainNewModel = false; // true if training, false if using an existing model
  bool useCustomModel = false; // true if using a custom model, false if using CBICA's CaPTk model

  if (m_Controls.m_cbUsageSelector->currentText() == "Train New Model") {
      trainNewModel = true;
      useCustomModel = false;
  }
  else if (m_Controls.m_cbUsageSelector->currentText() == "Use Existing Model") {
      trainNewModel = false;

      if  (m_Controls.m_cbModelSourceSelector->currentText() == "CBICA CaPTk Model") {
          useCustomModel = false;
      }
      else if (m_Controls.m_cbModelSourceSelector->currentText() == "Custom") {
          useCustomModel = true;
      }

  }
  m_Controls.pushButtonRun->setDisabled(true);
  m_Controls.pushButtonRun->setText("Running Recurrence Prediction...");
  
  modelDirPath = modelDirPath;
  subjectDirPath = subjectDirPath;
  outputDirPath = outputDirPath;
  trainNewModel = trainNewModel;
  useCustomModel = useCustomModel;
  // m_CaPTkRecurrence->Run(
  //  modelDirPath,
  //  subjectDirPath,
  //  outputDirPath,
  //  trainNewModel,
  //  useCustomModel
  //  );
}

void QmitkCaPTkRecurrenceView::OnModuleDone()
{
    m_Controls.pushButtonRun->setEnabled(true);
    m_Controls.pushButtonRun->setText("Run Recurrence Prediction");
}

/************************************************************************/
/* protected                                                            */
/************************************************************************/
void QmitkCaPTkRecurrenceView::OnSelectionChanged(berry::IWorkbenchPart::Pointer, const QList<mitk::DataNode::Pointer> &)
{

}

void QmitkCaPTkRecurrenceView::OnPreferencesChanged(const berry::IBerryPreferences*)
{

}

void QmitkCaPTkRecurrenceView::NodeAdded(const mitk::DataNode *)
{

}

void QmitkCaPTkRecurrenceView::NodeRemoved(const mitk::DataNode *)
{

}

void QmitkCaPTkRecurrenceView::SetFocus()
{
}

void QmitkCaPTkRecurrenceView::UpdateControls()
{
  this->RequestRenderWindowUpdate(mitk::RenderingManager::REQUEST_UPDATE_ALL);
}

void QmitkCaPTkRecurrenceView::InitializeListeners()
{

}

QString QmitkCaPTkRecurrenceView::GetLastFileOpenPath()
{
  return this->GetPreferences()->Get("LastFileOpenPath", "");
}

void QmitkCaPTkRecurrenceView::SetLastFileOpenPath(const QString &path)
{
  this->GetPreferences()->Put("LastFileOpenPath", path);
  this->GetPreferences()->Flush();
}

QString QmitkCaPTkRecurrenceView::GetModelPath()
{
  QString modelPath;

  if (m_Controls.m_cbModelSourceSelector->currentText() == "CBICA CaPTk Model")
  {
    modelPath = QCoreApplication::applicationDirPath() + QString("/models/recurrence_model");
  }
  else
  {
    modelPath = m_Controls.lineEdit_CustomModelDir->text();
  }

  return modelPath;
}

void QmitkCaPTkRecurrenceView::SelectedModelChanged()
{
  MITK_INFO << "QmitkCaPTkRecurrenceView::SelectedModelChanged" << this->GetModelPath();

  QFile file(this->GetModelPath() + "/requirements.txt");
  if(!file.open(QIODevice::ReadOnly)) {
      // QMessageBox::information(0, 
      //       "This model does not have the correct requirements.txt file", 
      //       file.errorString()
      // );
      m_Controls.label_ModelRequirements->setText("requirements.txt file does not exist for this model.");
      return;
  }

  QTextStream in(&file);
  QString all;

  while(!in.atEnd()) {
      all += "\n • " + in.readLine();    
  }

  m_Controls.label_ModelRequirements->setText(all);

  file.close();
}
