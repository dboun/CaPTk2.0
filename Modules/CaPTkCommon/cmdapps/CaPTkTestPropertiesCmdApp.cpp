/** \file CaPTkTestPropertiesCmdApp.cpp
 * \brief CLI program for testing properties
 */

#include <mitkCommandLineParser.h>

#include "CaPTkCohortOperations.h"

#include <mitkException.h>
#include <mitkLogMacros.h>
#include <mitkIOUtil.h>
#include <mitkImage.h>

#include <QString>
#include <QStringList>
#include <QFile>
#include <QJsonDocument>
#include <QSharedPointer>

#include <algorithm>
#include <string>
#include <fstream>
#include <iostream>
#include <memory>

/** \brief command-line app for testing properties
 *
 * This command-line app takes a list of directories and merges them to create a cohort json file.
 */
int main(int /*argc*/, char *argv[])
{
    QString image_path = QString(argv[1]);
    std::cout << image_path.toStdString() << "\n";

    auto image = mitk::IOUtil::Load<mitk::Image>(image_path.toStdString());

    std::cout << "\n\n\nPROPERTIES:\n";

    for (auto &prop : image->GetPropertyKeys())
    {
        std::cout << prop << ": " << image->GetProperty(prop.c_str())->GetValueAsString() << "\n";
    }

    std::cout << "\n\n";
    std::cout << image->GetProperty("name")->GetValueAsString() << "\n";
    std::cout << image->GetProperty("modality")->GetValueAsString() << "\n";
    std::cout << image->GetProperty("dicom.series.SeriesDescription")->GetValueAsString() << "\n";
}