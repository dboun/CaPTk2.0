/**
\file CaPTkPseudoProgressionEstimator.h
This file holds the declaration of the class CaPTkPseudoProgressionEstimator.
https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu
Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html
*/

#ifndef CaPTkPseudoProgressionEstimator_h
#define CaPTkPseudoProgressionEstimator_h

#include "NiftiDataManager.h"
#include "OutputWritingManager.h"
#include "FeatureReductionClass.h"
#include "FeatureScalingClass.h"
#include "FeatureExtractionClass.h"
#include "fProgressDialog.h"
#include "itkMeanImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkShapeLabelObject.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkScalarImageToCooccurrenceMatrixFilter.h"
#include "itkHistogramToTextureFeaturesFilter.h"
#include "itkScalarImageToRunLengthFeaturesFilter.h"
#include "itkCSVArray2DFileReader.h"
#include "itkConnectedComponentImageFilter.h"
#include "CaPTkEnums.h"
#include "CaPTkClassifierUtils.h"
#include "cbicaLogging.h"
#include "itkEnhancedScalarImageToRunLengthFeaturesFilter.h"
#include "itkRoundImageFilter.h"

#define RECURRENCE_MODEL_G 0.5
#define RECURRENCE_MODEL_RHO 0.0896

#define NO_OF_PCS 5 // total number of principal components used

typedef itk::Image<float, 3> ImageType;
typedef itk::Image<float, 4> PerfusionImageType;
typedef itk::CSVArray2DFileReader<double> CSVFileReaderType;
typedef vnl_matrix<double> MatrixType;

typedef std::tuple<std::vector<ImageType::IndexType>, VariableSizeMatrixType> PerfusionTupleType;
typedef std::map<int, PerfusionTupleType> PerfusionMapType;

namespace captk
{
/**
\class PseudoProgressionEstimator
\brief PseudoProgression Estimation for glioblastomas

References:
@inproceedings{,
title={Imaging Surrogates of Infiltration Obtained Via Multiparametric Imaging Pattern Analysis Predict Subsequent Location of Recurrence of Glioblastoma},
author={Akbari, Hamed MD, PhD; Macyszyn, Luke MD, MA; Da, Xiao MSc; Bilello, Michel MD, PhD; Wolf, Ronald L. MD, PhD; Martinez-Lage, Maria MD; Biros, George PhD; Alonso-Basanta, Michelle MD, PhD; O'Rourke, Donald M. MD; Davatzikos, Christos PhD},
pages={572�580},
year={2016},
organization={Neurosurgery}
}
@inproceedings{,
title={Pattern Analysis of Dynamic Susceptibility Contrast-enhanced MR Imaging Demonstrates Peritumoral Tissue Heterogeneity},
author={Hamed Akbari, MD, PhD Luke Macyszyn, MD, MA Xiao Da, MS Ronald L. Wolf, MD, PhD Michel Bilello, MD, PhD Ragini Verma, PhD Donald M. O�Rourke, MD Christos Davatzikos, PhD},
pages={502-510},
year={2014},
organization={RSNA Radiology}
}
*/
class PseudoProgressionEstimator
{
public:
    PseudoProgressionEstimator();

    ~PseudoProgressionEstimator();

    std::string mLastEncounteredError; /// Gets last error
    std::string mPseudoTrainedFile;
    std::string mRecurrenceTrainedFile;
    cbica::Logging logger;

    void WriteCSVFiles(VariableSizeMatrixType inputdata, std::string filepath);
    void WriteCSVFiles(vtkSmartPointer<vtkTable> inputdata, std::string filepath);
    void WriteCSVFiles(VectorVectorDouble inputdata, std::string filepath);
    void WriteCSVFiles(VariableLengthVectorType inputdata, std::string filepath);
    void WriteCSVFiles(std::vector<double> inputdata, std::string filepath);

    VariableSizeMatrixType GetModelSelectedFeatures(
        VariableSizeMatrixType & ScaledFeatureSetAfterAddingLabel, 
        VariableLengthVectorType & psuSelectedFeatures
    );
    
    void WritePCAOutputs(
        std::string suffix, 
        std::string outputdirectory, 
        const VariableLengthVectorType mean, 
        const VariableSizeMatrixType coefs
    );

    VectorVectorDouble CombineAllThePerfusionFeaures(
        VectorVectorDouble T1IntensityHistogram,
        VectorVectorDouble TCIntensityHistogram,
        VectorVectorDouble T1TCIntensityHistogram,
        VectorVectorDouble T2IntensityHistogram,
        VectorVectorDouble FLIntensityHistogram,
        VectorVectorDouble T2FLIntensityHistogram,
        VectorVectorDouble AXIntensityHistogram,
        VectorVectorDouble FAIntensityHistogram,
        VectorVectorDouble RDIntensityHistogram,
        VectorVectorDouble TRIntensityHistogram,
        VectorVectorDouble PHIntensityHistogram,
        VectorVectorDouble PSIntensityHistogram,
        VectorVectorDouble RCIntensityHistogram,
        VectorVectorDouble PCA1IntensityHistogram,
        VectorVectorDouble PCA2IntensityHistogram,
        VectorVectorDouble PCA3IntensityHistogram,
        VectorVectorDouble PCA4IntensityHistogram,
        VectorVectorDouble PCA5IntensityHistogram,
        VectorVectorDouble PCA6IntensityHistogram,
        VectorVectorDouble PCA7IntensityHistogram,
        VectorVectorDouble PCA8IntensityHistogram,
        VectorVectorDouble PCA9IntensityHistogram,
        VectorVectorDouble PCA10IntensityHistogram,
        std::string output
    );

    vtkSmartPointer<vtkTable> MakePCAMatrix(int NumberOfFeatures, int NumberOfSamples);

    VariableSizeMatrixType ColumnWiseScaling(VariableSizeMatrixType PerfusionData);

    // ---- PERFUSION & PCA ----

    PerfusionMapType CombineAndCalculatePerfusionPCA(
        PerfusionMapType PerfusionDataMap, 
        VariableSizeMatrixType & TransformationMatrix, 
        VariableLengthVectorType & MeanVector
    );
    
    PerfusionMapType CombineAndCalculatePerfusionPCAForTestData(
        PerfusionMapType PerfusionDataMap, 
        VariableSizeMatrixType & TransformationMatrix,
        VariableLengthVectorType & MeanVector
    );

    template <class PerfusionImageType, class ImageType>
    VariableSizeMatrixType LoadPerfusionData(
        typename ImageType::Pointer maskImagePointerNifti, 
        typename PerfusionImageType::Pointer perfImagePointerNifti, 
        std::vector<typename ImageType::IndexType> & indices
    );

    // ---- ... ----

    VectorDouble CombineEstimates(
        const VariableLengthVectorType &estimates1, 
        const VariableLengthVectorType &estimates2
    );

    VariableLengthVectorType DistanceFunction(
        const VariableSizeMatrixType &testData, 
        const std::string &filename, 
        const double &rho, const double &bestg
    );

    int GetFeatureVectorSize(
        bool &useConvetionalData, 
        bool &useDTIData, 
        bool &usePerfData, 
        bool &useDistData
    );

    /** \brief Image Intensity using the itk::RescaleIntensityImageFilter to the 0-255 range */
    ImageTypeFloat3D::Pointer RescaleImageIntensity(ImageTypeFloat3D::Pointer image);

    template <class ImageType>
    std::tuple<VectorDouble, VectorDouble, VectorDouble, VectorDouble, VectorDouble> 
    GetAllFeaturesPerImagePerROI(
        typename ImageType::Pointer image, 
        typename ImageType::Pointer mask, 
        std::string modality
    );

    PerfusionMapType CombinePerfusionDataAndApplyExistingPerfusionModel(
        PerfusionMapType PerfusionDataMap, 
        VariableSizeMatrixType TransformationMatrix, 
        VariableLengthVectorType MeanVector
    );

    void ReadAllTheModelParameters(
        std::string modeldirectory,
        VariableSizeMatrixType & PCA_PERF,
        VariableSizeMatrixType & PCA_T1,
        VariableSizeMatrixType & PCA_T1CE,
        VariableSizeMatrixType & PCA_T2,
        VariableSizeMatrixType & PCA_FL,
        VariableSizeMatrixType & PCA_T1T1CE,
        VariableSizeMatrixType & PCA_T2FL,
        VariableSizeMatrixType & PCA_AX,
        VariableSizeMatrixType & PCA_FA,
        VariableSizeMatrixType & PCA_RAD,
        VariableSizeMatrixType & PCA_TR,
        VariableSizeMatrixType & PCA_PH,
        VariableSizeMatrixType & PCA_PSR,
        VariableSizeMatrixType & PCA_RCBV,
        VariableSizeMatrixType & PCA_PC1,
        VariableSizeMatrixType & PCA_PC2,
        VariableSizeMatrixType & PCA_PC3,
        VariableSizeMatrixType & PCA_PC4,
        VariableSizeMatrixType & PCA_PC5,
        VariableSizeMatrixType & PCA_PC6,
        VariableSizeMatrixType & PCA_PC7,
        VariableSizeMatrixType & PCA_PC8,
        VariableSizeMatrixType & PCA_PC9,
        VariableSizeMatrixType & PCA_PC10,
        VariableLengthVectorType & Mean_PERF,
        VariableLengthVectorType & Mean_T1,
        VariableLengthVectorType & Mean_T1CE,
        VariableLengthVectorType & Mean_T2,
        VariableLengthVectorType & Mean_FL,
        VariableLengthVectorType & Mean_T1T1CE,
        VariableLengthVectorType & Mean_T2FL,
        VariableLengthVectorType & Mean_AX,
        VariableLengthVectorType & Mean_FA,
        VariableLengthVectorType & Mean_RAD,
        VariableLengthVectorType & Mean_TR,
        VariableLengthVectorType & Mean_PH,
        VariableLengthVectorType & Mean_PSR,
        VariableLengthVectorType & Mean_RCBV,
        VariableLengthVectorType & Mean_PC1,
        VariableLengthVectorType & Mean_PC2,
        VariableLengthVectorType & Mean_PC3,
        VariableLengthVectorType & Mean_PC4,
        VariableLengthVectorType & Mean_PC5,
        VariableLengthVectorType & Mean_PC6,
        VariableLengthVectorType & Mean_PC7,
        VariableLengthVectorType & Mean_PC8,
        VariableLengthVectorType & Mean_PC9,
        VariableLengthVectorType & Mean_PC10
    );

    /** \brief Train new model and save it to a directory */
    bool TrainNewModelOnGivenData(
        const std::vector<std::map<CAPTK::ImageModalityType, std::string>> qualifiedsubjects, 
        const std::string &outputdirectory, 
        bool useConventioalData, bool useDTIData, bool usePerfData, bool useDistData
    );

    VariableSizeMatrixType LoadPseudoProgressionTrainingData(
        const std::vector<std::map<CAPTK::ImageModalityType, std::string>> &trainingsubjects, 
        std::vector<double> &traininglabels, 
        std::string outputdirectory
    );

    VariableSizeMatrixType LoadPseudoProgressionTestingData(
        const std::vector<std::map<CAPTK::ImageModalityType, std::string>> &trainingsubjects, 
        std::vector<double> &traininglabels, 
        std::string outputdirectory, 
        std::string modeldirectory
    );

    /** \brief Recurrence Estimation using existing model */
    bool PseudoProgressionEstimateOnExistingModel(
        std::vector<std::map<CAPTK::ImageModalityType, std::string>> qualifiedsubjects,
        const std::string &modeldirectory,
        const std::string &inputdirectory,
        const std::string &outputdirectory,
        bool useConventionalData, bool useDTIData, bool usePerfData, bool useDistData
    );

    template <typename T>
    std::vector<size_t> sort_indexes(const std::vector<T> &v);

    VariableSizeMatrixType SelectModelFeatures(const VariableSizeMatrixType &ModelFeatures);
    
    VectorDouble EffectSizeFeatureSelection(
        const VariableSizeMatrixType training_features, 
        std::vector<double> target
    );

    template <class ImageType>
    typename ImageType::Pointer ReadNiftiImage(const std::string &filename);

    template <class ImageType>
    typename ImageType::Pointer RescaleImageIntensity(const typename ImageType::Pointer &image);

    VectorDouble GetHistogramFeatures(std::vector<float> intensities, int m_Bins);
    VectorDouble GetHistogramFeaturesWhole(std::vector<float> intensities);

    template <class TImageTypeShape = TImageType>
    VectorDouble GetShapeFeatures(typename TImageTypeShape::Pointer mask);

    VectorDouble GetIntensityFeatures(std::vector<float> m_intensities);

    template <class ImageType>
    VectorDouble GetGLCMFeatures(
        typename ImageType::Pointer image, 
        typename ImageType::Pointer mask, 
        double minvalue, 
        double maxvalue
    );

    template <class ImageType>
    typename ImageType::Pointer MakeAdditionalModality(
        typename ImageType::Pointer image1, 
        typename ImageType::Pointer image2
    );

    template <class ImageType>
    VectorDouble GetRunLengthFeatures(
        typename ImageType::Pointer image, 
        typename ImageType::Pointer mask, 
        double minvalue, double maxvalue
    );

    template <class PerfusionImageType, class ImageType>
    VectorVectorDouble GetPerfusionFeatures(
        typename PerfusionImageType::Pointer image, 
        typename ImageType::Pointer mask
    );

    /** \brief Postprocessing of a recurrence map */
    template <class ImageType>
    VectorDouble RecurrenceMapPostprocessing(
        VectorDouble result, 
        std::vector<typename ImageType::IndexType> indices, 
        typename ImageType::Pointer RecurrenceProbabilityMap, 
        typename ImageType::Pointer edemaMap
    );

    /** \brief Recurrence Estimation using existing model on given subject */
    template <class ImageType>
    void PseudoProgressionEstimateOnGivenSubject(
        typename ImageType::Pointer edemaMask,
        typename ImageType::Pointer tumorMask,
        typename ImageType::Pointer FinalT1CEImagePointer,
        typename ImageType::Pointer FinalT2FlairImagePointer,
        typename ImageType::Pointer FinalT1ImagePointer,
        typename ImageType::Pointer FinalT2ImagePointer,
        std::vector<typename ImageType::Pointer> FinalPerfusionImagePointer,
        std::vector<typename ImageType::Pointer> FinalDTIImagePointer,
        int imagetype, 
        bool conDataPresent, 
        bool dtiDataPresent, 
        bool perfusionDataPresent, 
        bool distanceDataPresent,
        bool useOtherModalities, 
        std::string t1cebasefilename,
        const VectorVectorDouble &nearRegionIndices,
        const VectorVectorDouble &farRegionIndices
    );

    template <class ImageType>
    void PseudoProgressionEstimateOnGivenSubjectUsingExistingModel(
        typename ImageType::Pointer edemaMask,
        typename ImageType::Pointer tumorMask,
        typename ImageType::Pointer FinalT1CEImagePointer,
        typename ImageType::Pointer FinalT2FlairImagePointer,
        typename ImageType::Pointer FinalT1ImagePointer,
        typename ImageType::Pointer FinalT2ImagePointer,
        std::vector<typename ImageType::Pointer> FinalPerfusionImagePointer,
        std::vector<typename ImageType::Pointer> FinalDTIImagePointer,
        int imagetype, bool conDataPresent, bool dtiDataPresent, bool perfusionDataPresent, bool distanceDataPresent,
        bool useOtherModalities, std::string t1cebasefilename,
        const VectorVectorDouble &nearRegionIndices,
        const VectorVectorDouble &farRegionIndices, std::string modeldirectory
    );

    NiftiDataManager mNiftiLocalPtr;
    OutputWritingManager mOutputLocalPtr;
    FeatureReductionClass mFeatureReductionLocalPtr;
    FeatureScalingClass mFeatureScalingLocalPtr;
    FeatureExtractionClass mFeatureExtractionLocalPtr;

    template <class ImageType>
    void Run(
        typename ImageType::Pointer edemaMask,
        typename ImageType::Pointer tumorMask,
        typename ImageType::Pointer FinalT1CEImagePointer,
        typename ImageType::Pointer FinalT2FlairImagePointer,
        typename ImageType::Pointer FinalT1ImagePointer,
        typename ImageType::Pointer FinalT2ImagePointer,
        std::vector<typename ImageType::Pointer> FinalPerfusionImagePointer,
        std::vector<typename ImageType::Pointer> FinalDTIImagePointer,
        int imagetype, 
        bool conDataPresent, bool dtiDataPresent, 
        bool perfusionDataPresent, bool distanceDataPresent,
        bool useOtherModalities, 
        const std::string &t1cebasefilename,
        const VectorVectorDouble &nearRegionIndices,
        const VectorVectorDouble &farRegionIndices, 
        const std::string &outputDirName)
    {
        mCurrentOutputDir = outputDirName;
        PseudoProgressionEstimateOnGivenSubject<ImageType>(
            edemaMask, tumorMask, 
            FinalT1CEImagePointer, FinalT2FlairImagePointer, 
            FinalT1ImagePointer, FinalT2ImagePointer,
            FinalPerfusionImagePointer, FinalDTIImagePointer,
            imagetype, conDataPresent, dtiDataPresent, 
            perfusionDataPresent, distanceDataPresent,
            useOtherModalities, t1cebasefilename, 
            nearRegionIndices, farRegionIndices
        );
    }

    template <class TImageType>
    typename TImageType::Pointer 
    GetImageWithLabels(std::vector<double> labels, typename TImageType::Pointer inputimage)
    {
        typename TImageType::Pointer output = TImageType::New();
        output->CopyInformation(inputimage);
        output->SetRequestedRegion(inputimage->GetLargestPossibleRegion());
        output->SetBufferedRegion(inputimage->GetBufferedRegion());
        output->Allocate();
        output->FillBuffer(0);

        typedef itk::ImageRegionIteratorWithIndex<TImageType> IteratorType;
        IteratorType maskIt(inputimage, inputimage->GetLargestPossibleRegion());
        IteratorType outputIt(output, output->GetLargestPossibleRegion());
        maskIt.GoToBegin();
        outputIt.GoToBegin();
        while (!maskIt.IsAtEnd())
        {
            for (int j = 0; j < labels.size(); j++)
            {
                if (maskIt.Get() == labels[j])
                {
                    outputIt.Set(CAPTK::VOXEL_STATUS::ON);
                }
            }
            ++maskIt;
            ++outputIt;
        }
        return output;
    }

    template <class ImageType>
    void RunLoadedSubjectOnExistingModel(
        typename ImageType::Pointer edemaMask,
        typename ImageType::Pointer tumorMask,
        typename ImageType::Pointer FinalT1CEImagePointer,
        typename ImageType::Pointer FinalT2FlairImagePointer,
        typename ImageType::Pointer FinalT1ImagePointer,
        typename ImageType::Pointer FinalT2ImagePointer,
        std::vector<typename ImageType::Pointer> FinalPerfusionImagePointer,
        std::vector<typename ImageType::Pointer> FinalDTIImagePointer,
        int imagetype, 
        bool conDataPresent, bool dtiDataPresent, 
        bool perfusionDataPresent, bool distanceDataPresent,
        bool useOtherModalities, 
        const std::string &t1cebasefilename,
        const VectorVectorDouble &nearRegionIndices,
        const VectorVectorDouble &farRegionIndices, 
        const std::string &outputDirName, 
        const std::string &modeldirectory)
    {
        mCurrentOutputDir = outputDirName;
        PseudoProgressionEstimateOnGivenSubjectUsingExistingModel<ImageType>(
            edemaMask, tumorMask, FinalT1CEImagePointer, 
            FinalT2FlairImagePointer, FinalT1ImagePointer, FinalT2ImagePointer,
            FinalPerfusionImagePointer, FinalDTIImagePointer,
            imagetype, conDataPresent, dtiDataPresent, perfusionDataPresent, distanceDataPresent,
            useOtherModalities, t1cebasefilename, nearRegionIndices, farRegionIndices, modeldirectory
        );
    }

    std::string mRecurrenceMapFileName;

private:
    std::string mCurrentOutputDir;
};

}

#include "CaPTkPseudoProgressionEstimator.hxx"

#endif