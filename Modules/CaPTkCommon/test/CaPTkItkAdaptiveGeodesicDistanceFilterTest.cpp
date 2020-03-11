#include "mitkTestingMacros.h"
#include <mitkTestingConfig.h>
#include "mitkTestFixture.h"

#include "itkImageFileReader.h"
#include <itkTestingComparisonImageFilter.h>

#include <iostream>
#include <fstream>

#include "CaPTkData.h"
#include "CaPTkItkAdaptiveGeodesicDistanceFilter.h"

class CaPTkItkAdaptiveGeodesicDistanceFilterTestSuite : public mitk::TestFixture
{

    CPPUNIT_TEST_SUITE(CaPTkItkAdaptiveGeodesicDistanceFilterTestSuite);
    MITK_TEST(Try2D);
    // MITK_TEST(Try2DFloatVsInt);
    // MITK_TEST(Try2DWithMask);
    // MITK_TEST(Try3D);
    // MITK_TEST(Try4D);
    CPPUNIT_TEST_SUITE_END();
    
    using ImageTypeFloat2D = itk::Image<float, 2>; 
    using ImageTypeInt2D   = itk::Image<int, 2>; 
    using ImageTypeFloat3D = itk::Image<float, 3>;
    using ImageTypeInt3D   = itk::Image<int, 3>;

private:
    std::string CAPTK_DATA_DIR;

    /*---- Input and expected output images ----*/

    ImageTypeFloat2D::Pointer   inputImageFloat2D;
    ImageTypeInt2D::Pointer     inputImageInt2D;
    ImageTypeFloat3D::Pointer   inputImageFloat2Das3D;
    ImageTypeFloat3D::Pointer   inputImageFloat3D;
    
    ImageTypeInt2D::Pointer     inputLabelsInt2D;
    ImageTypeFloat2D::Pointer   inputLabelsFloat2D;
    ImageTypeFloat3D::Pointer   inputLabelsFloat2Das3D;
    ImageTypeInt3D::Pointer     inputLabelsInt3D;
    
    ImageTypeFloat2D::Pointer   outputImageFloat2D;
    ImageTypeInt2D::Pointer     outputImageInt2D;
    ImageTypeInt3D::Pointer     outputImageInt2Das3D;
    ImageTypeFloat3D::Pointer   outputImageFloat3D;
    
    ImageTypeFloat2D::Pointer   maskImageFloat2D;
    ImageTypeFloat3D::Pointer   maskImageFloat3D;

    /** \brief Helper function to read an image */
    template <class TImageType>
    static typename TImageType::Pointer ReadImage(const std::string& path)
    {
        using ReaderType = itk::ImageFileReader<TImageType>;
        typename ReaderType::Pointer reader = ReaderType::New();
        reader->SetFileName(path);
        reader->Update();
        return reader->GetOutput();
    }

    /** \brief Helper function to check if every pixel pair of two images are the same */
    template <class TImageType>
    static bool Equal(const typename TImageType::Pointer validImage, 
               const typename TImageType::Pointer testImage)
    {
        using DiffType = itk::Testing::ComparisonImageFilter<TImageType, TImageType>;
        typename DiffType::Pointer diff = DiffType::New();
        diff->SetToleranceRadius( 0 );
        diff->SetDifferenceThreshold( 0.001 ); // So that comparing int to floats is fine
        diff->VerifyInputInformationOn();

        diff->SetValidInput( validImage );
        diff->SetTestInput( testImage );
        diff->UpdateLargestPossibleRegion();
        // averageIntensityDifference kept for reference
        // const double averageIntensityDifference = 
        //         diff->GetTotalDifference();
        const unsigned long numberOfPixelsWithDifferences = 
                diff->GetNumberOfPixelsWithDifferences();

        return (numberOfPixelsWithDifferences == 0);
    }

public:
    void setUp() override
    {
        CAPTK_DATA_DIR = captk::data::absolute_search_dirs[0];

        /*---- Read input and expected output images ----*/

        try
        {
            inputImageFloat2D = ReadImage<ImageTypeFloat2D>(CAPTK_DATA_DIR 
                    + std::string("/datasets/ExampleSmall/Subject_2D/small2D.nii.gz"));
            inputImageInt2D = ReadImage<ImageTypeInt2D>(CAPTK_DATA_DIR 
                    + std::string("/datasets/ExampleSmall/Subject_2D/small2D.nii.gz"));
            inputImageFloat2Das3D = ReadImage<ImageTypeFloat3D>(CAPTK_DATA_DIR 
                    + std::string("/datasets/ExampleSmall/Subject_2D/small2D.nii.gz"));
            inputImageFloat3D = ReadImage<ImageTypeFloat3D>(CAPTK_DATA_DIR 
                    + std::string("/datasets/ExampleSmall/Subject_3D/small3D.nii.gz"));

            inputLabelsInt2D = ReadImage<ImageTypeInt2D>(CAPTK_DATA_DIR 
                    + std::string("/datasets/ExampleSmall/Subject_2D/small2D-labels.nii.gz"));
            inputLabelsFloat2D = ReadImage<ImageTypeFloat2D>(CAPTK_DATA_DIR 
                    + std::string("/datasets/ExampleSmall/Subject_2D/small2D-labels.nii.gz"));
            inputLabelsFloat2Das3D = ReadImage<ImageTypeFloat3D>(CAPTK_DATA_DIR 
                    + std::string("/datasets/ExampleSmall/Subject_2D/small2D-labels.nii.gz"));
            inputLabelsInt3D = ReadImage<ImageTypeInt3D>(CAPTK_DATA_DIR 
                    + std::string("/datasets/ExampleSmall/Subject_3D/small3D-labels.nii.gz"));

            outputImageFloat2D = ReadImage<ImageTypeFloat2D>(CAPTK_DATA_DIR 
                    + std::string("/test-data/Common/AdaptiveGeodesicDistance/"
                      "small2D-agd-result-lb2-limit255-imageasmask.nii.gz"));
            outputImageInt2D = ReadImage<ImageTypeInt2D>(CAPTK_DATA_DIR 
                    + std::string("/test-data/Common/AdaptiveGeodesicDistance/"
                      "small2D-agd-result-lb2-limit255-imageasmask.nii.gz"));
            outputImageInt2Das3D = ReadImage<ImageTypeInt3D>(CAPTK_DATA_DIR 
                    + std::string("/test-data/Common/AdaptiveGeodesicDistance/"
                      "small2D-agd-result-lb2-limit255-imageasmask.nii.gz"));
            outputImageFloat3D = ReadImage<ImageTypeFloat3D>(CAPTK_DATA_DIR 
                    + std::string("/test-data/Common/AdaptiveGeodesicDistance/"
                      "small3D-agd-result-lb1-limit255-imageasmask.nii.gz"));

            // TODO: maskImageFloat2D (also uncomment stuff below)
            // TODO: maskImageFloat3D (also uncomment stuff below)
        }
        catch (const itk::ExceptionObject &err)
        {
            CPPUNIT_FAIL(std::string(
                             "Reading on of the input images failed. "
                             "This might be because you are on an older "
                             "version of the CaPTkData. "
                             "If that is the case, please rebuild your superbuild. "
                             "Exception was: ") +
                         err.what());
        }
    }

    void tearDown() override
    {

    }

    void Try2D()
    {
        /*---- Run filter ----*/

        using FilterType = 
            itk::captk::AdaptiveGeodesicDistanceFilter< ImageTypeFloat2D, ImageTypeFloat2D >;
        FilterType::Pointer filter = FilterType::New();

        filter->SetInput( inputImageFloat2D );
        filter->SetLabels( inputLabelsInt2D );
        filter->SetLabelOfInterest( 2 );
        filter->LimitAt255On();

        try
        {
            filter->Update();
        }
        catch (const itk::ExceptionObject &err)
        {
            CPPUNIT_FAIL(std::string("Updating the filter failed. Exception was: ") +
                         err.what());
        }

        /*---- Compare output image with the "ground truth" ----*/

        CPPUNIT_ASSERT_MESSAGE( "Valid & Test images are not equal", 
                Equal<ImageTypeFloat2D>(outputImageFloat2D, filter->GetOutput())
        );
    }
};

MITK_TEST_SUITE_REGISTRATION(CaPTkItkAdaptiveGeodesicDistanceFilter)
