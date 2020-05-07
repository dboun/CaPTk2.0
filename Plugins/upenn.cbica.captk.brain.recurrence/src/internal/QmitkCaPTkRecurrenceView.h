#ifndef QmitkCaPTkRecurrenceView_h
#define QmitkCaPTkRecurrenceView_h

#include <QmitkAbstractView.h>

#include <mitkILifecycleAwarePart.h>

#include "ui_QmitkCaPTkRecurrenceControls.h"
#include "CaPTkRecurrence.h"

class CaPTkRecurrence;

// berry
#include <berryIBerryPreferences.h>

class QmitkRenderWindow;

/**
 * \ingroup ToolManagerEtAl
 * \ingroup upenn_cbica_captk_brain_recurrence_internal
 */
class QmitkCaPTkRecurrenceView : public QmitkAbstractView, public mitk::ILifecycleAwarePart
{
  Q_OBJECT

public:
  static const std::string VIEW_ID;

  QmitkCaPTkRecurrenceView();
  virtual ~QmitkCaPTkRecurrenceView();

  // GUI setup
  void CreateQtPartControl(QWidget *parent) override;

  // ILifecycleAwarePart interface
public:
  void Activated() override;
  void Deactivated() override;
  void Visible() override;
  void Hidden() override;

  virtual int GetSizeFlags(bool width);
  virtual int ComputePreferredSize(bool width,
                                   int /*availableParallel*/,
                                   int /*availablePerpendicular*/,
                                   int preferredResult);  
protected slots:

  void OnUsageComboBoxCurrentTextChanged(const QString& text);

  void OnModelSourceComboBoxCurrentTextChanged(const QString& text);

  void OnSubjectDirectoryButtonClicked();

  void OnCustomModelDirectoryButtonClicked();

  void OnModelPathLineEditTextChanged(QString);

  void OnOutputDirectoryButtonClicked();

  void OnModuleDone();

  /** \brief CaPTk Recurrence Plugin Run Button clicked slot */
  void OnRunButtonPressed();


protected:

  // reimplemented from QmitkAbstractView
  void OnSelectionChanged(berry::IWorkbenchPart::Pointer part, const QList<mitk::DataNode::Pointer> &nodes) override;

  // reimplemented from QmitkAbstractView
  void OnPreferencesChanged(const berry::IBerryPreferences* prefs) override;

  // reimplemented from QmitkAbstractView
  void NodeAdded(const mitk::DataNode* node) override;

  // reimplemented from QmitkAbstractView
  void NodeRemoved(const mitk::DataNode* node) override;

  void SetFocus() override;

  void UpdateControls();

  void InitializeListeners();

  QString GetLastFileOpenPath();

  void SetLastFileOpenPath(const QString &path);

  /// \brief Returns the selected model path. It can distinguish between CBICA and custom
  QString GetModelPath();

  /// \brief This is called when the selected model is changed in the UI, to update the requirements
  void SelectedModelChanged();

  /// \brief the Qt parent of our GUI (NOT of this object)
  QWidget *m_Parent;

  /// \brief Qt GUI file
  Ui::QmitkCaPTkRecurrenceControls m_Controls;

  captk::CaPTkRecurrence* m_CaPTkRecurrence;
};

#endif // QmitkCaPTkRecurrenceView_h
