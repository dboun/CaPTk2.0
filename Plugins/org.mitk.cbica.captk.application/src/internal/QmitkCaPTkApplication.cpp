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

#include "QmitkCaPTkApplication.h"

#include <berryPlatformUI.h>

#include <berryPlatform.h>
#include "berryPlatform.h"
#include "berryPlatformUI.h"
#include "berryIWorkbench.h"
#include "berryIConfigurationElement.h"
#include "berryIExtensionRegistry.h"
#include "berryIExtension.h"
#include <berryIBerryPreferencesService.h>
#include <berryIQtPreferencePage.h>

#include "QmitkCaPTkAppWorkbenchAdvisor.h"

QmitkCaPTkApplication::QmitkCaPTkApplication()
{

}

QVariant QmitkCaPTkApplication::Start(berry::IApplicationContext* /*context*/)
{
  berry::Platform::GetPreferencesService()->GetSystemPreferences()->Node(
    "org.mitk.editors.stdmultiwidget"
  )->Put("DepartmentLogo", ":/org.mitk.cbica.captk.application/cbica-logo.jpg");

  QScopedPointer<berry::Display> display(berry::PlatformUI::CreateDisplay());

  QScopedPointer<QmitkCaPTkAppWorkbenchAdvisor> wbAdvisor(new QmitkCaPTkAppWorkbenchAdvisor());
  int code = berry::PlatformUI::CreateAndRunWorkbench(display.data(), wbAdvisor.data());

  // exit the application with an appropriate return code
  return code == berry::PlatformUI::RETURN_RESTART
              ? EXIT_RESTART : EXIT_OK;
}

void QmitkCaPTkApplication::Stop()
{

}
