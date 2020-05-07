#include "CaPTkRecurrenceTrainCheckInputStructure.h"

#include "mitkLogMacros.h"

#include <QDirIterator>

#include <set>

namespace captk 
{

MITKCAPTKRECURRENCE_EXPORT bool RecurrenceTrainCheckInputStructure(
    QSharedPointer<captk::Cohort> cohort,
    std::string& errorMessage)
{
    // The required modality + description combinations for training
    std::vector<std::tuple<QString, QString>> requiredImageTypes = {
        std::make_tuple<QString, QString>("seg", "reccurence"),
        std::make_tuple<QString, QString>("seg", "edema"),
        std::make_tuple<QString, QString>("seg", "tc")
    };

    return internal::RecurrenceCheckInputStructure(cohort, errorMessage, requiredImageTypes);
}

MITKCAPTKRECURRENCE_EXPORT bool RecurrenceInferenceCheckInputStructure(
    QSharedPointer<captk::Cohort> cohort,
    std::string& errorMessage,
    const std::vector<std::tuple<QString, QString>> requiredImageTypes)
{
    return internal::RecurrenceCheckInputStructure(cohort, errorMessage, requiredImageTypes);
}

bool internal::RecurrenceCheckInputStructure(
    QSharedPointer<captk::Cohort> cohort, 
    std::string& errorMessage, 
    const std::vector<std::tuple<QString, QString>> requiredImageTypes)
{
    auto subjects = cohort->GetSubjects();

    if (subjects.size() == 0)
    {
        errorMessage = "No subjects";
        return false;
    }

    // Check that only one study per subject
    {
        for (auto subj : subjects)
        {
            if (subj->GetStudies().size() != 1)
            {
                errorMessage = "Only one study per subject is permitted";
                return false;
            }
        }
    }

    // Parse first subject
    std::set<std::tuple<QString, QString>> imageTypesOfFirstSubject;
    {
        firstStudy = subjects[0]->GetStudies()[0];

        for (auto series : firstStudy->GetSeries())
        {
            auto modality = series->GetModality();
            auto desc = ((modality != "seg") ? 
                series->GetSeriesDescription() :
                series->GetSegmentLabel() 
            );

            imageTypesOfFirstSubject.insert(
                std::make_tuple<QString, QString>(modality, desc)
            );
        }
    }

    // Check that the mandatory image types exist for the first subject
    {
        for (auto req : requiredImageTypes)
        {
            bool found = false;

            std::set<std::tuple<QString, QString>>::iterator it;
            for (it = imageTypesOfFirstSubject.begin(); 
                 it != imageTypesOfFirstSubject.end(); 
                 ++it) 
            {
                modality = std::get<0>(*it);
                desc = std::get<1>(*it);

                if (std::get<0>(*it) == std::get<0>(req) && 
                    std::get<1>(*it) == std::get<1>(req))
                {
                    found = true;
                    break;
                }
            }

            if (!found)
            {
                errorMessage  = "Not all required image types found ";
                return false;
            }
        }
    }

    // Check that at least 1 actual feature exists for the first subject
    {
        if (subjects[0]->GetStudies()[0]->GetSeries().size() == 3)
        {
            errorMessage = "Please add at least one actual feature";
            return false;
        }
    }

    // Check that the rest of the subjects match the first
    {
        for (auto subj : subjects)
        {
            auto study = subj->GetStudies()[0];

            for (auto series : study->GetSeries())
            {
                auto modality = series->GetModality();
                auto desc = ((modality != "seg") ? 
                    series->GetSeriesDescription() :
                    series->GetSegmentLabel() 
                );

                auto t = std::make_tuple<QString, QString>(modality, desc);

                if (imageTypesOfFirstSubject.find(t) == imageTypesOfFirstSubject.end())
                {
                    errorMessage = "Image types don\'t match for all subjects";
                    return false;    
                }
            }
        }
    }
}

}