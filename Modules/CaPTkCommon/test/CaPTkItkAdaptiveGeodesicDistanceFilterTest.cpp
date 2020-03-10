#include "mitkTestingMacros.h"
#include <mitkTestingConfig.h>
#include "mitkTestFixture.h"

#include "itkImageFileReader.h"

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

private:
    std::string CAPTK_DATA_DIR;

    /*---- Input and expected output images ----*/

    ImageTypeFloat2D::Pointer   inputImageFloat2D;
    ImageTypeFloat2D::Pointer   inputImageInt2D;
    ImageTypeFloat3D::Pointer   inputImageFloat2Das3D;
    ImageTypeFloat3D::Pointer   inputImageFloat3D;
    
    ImageTypeInt2D::Pointer     inputLabelsInt2D;
    ImageTypeInt2D::Pointer     inputLabelsFloat2D;
    ImageTypeInt3D::Pointer     inputLabelsFloat2Das3D;
    ImageTypeInt3D::Pointer     inputLabelsInt3D;
    
    ImageTypeFloat2D::Pointer   outputImageFloat2D;
    ImageTypeFloat2D::Pointer   outputImageInt2D;
    ImageTypeFloat3D::Pointer   outputImageInt2Das3D;
    ImageTypeFloat3D::Pointer   outputImageFloat3D;
    
    ImageTypeFloat2D::Pointer   maskImageFloat2D;
    ImageTypeFloat3D::Pointer   maskImageFloat3D;

public:
    void
    setUp() override
    {
        CAPTK_DATA_DIR = captk::data::absolute_search_dirs[0];

        /*---- Read input and expected output images ----*/

        using ReaderTypeFloat2D = itk::ImageFileReader<ImageTypeFloat2D>;
        using ReaderTypeInt2D = itk::ImageFileReader<ImageTypeInt2D>;
        using ReaderTypeFloat3D = itk::ImageFileReader<ImageTypeFloat3D>;
        using ReaderTypeInt3D = itk::ImageFileReader<ImageTypeInt3D>;

        ReaderTypeFloat2D::Pointer reader1  = ReaderTypeFloat2D::New();
        ReaderTypeFloat2D::Pointer reader2  = ReaderTypeFloat2D::New();
        ReaderTypeFloat3D::Pointer reader3  = ReaderTypeFloat3D::New();
        ReaderTypeFloat3D::Pointer reader4  = ReaderTypeFloat3D::New();
        ReaderTypeInt2D::Pointer   reader5  = ReaderTypeInt2D::New();
        ReaderTypeInt2D::Pointer   reader6  = ReaderTypeInt2D::New();
        ReaderTypeInt3D::Pointer   reader7  = ReaderTypeInt3D::New();
        ReaderTypeInt3D::Pointer   reader8  = ReaderTypeInt3D::New();
        ReaderTypeFloat2D::Pointer reader9  = ReaderTypeFloat2D::New();
        ReaderTypeFloat2D::Pointer reader10 = ReaderTypeFloat2D::New();
        ReaderTypeFloat3D::Pointer reader11 = ReaderTypeFloat3D::New();
        ReaderTypeFloat3D::Pointer reader12 = ReaderTypeFloat3D::New();
        ReaderTypeFloat2D::Pointer reader13 = ReaderTypeFloat2D::New();
        ReaderTypeFloat3D::Pointer reader14 = ReaderTypeFloat3D::New();

        reader1->SetFileName(CAPTK_DATA_DIR + "/datasets/ExampleSmall/Subject_2D"
                                              "small2D.nii.gz");
        reader2->SetFileName(CAPTK_DATA_DIR + "/datasets/ExampleSmall/Subject_2D"
                                              "small2D.nii.gz");
        reader3->SetFileName(CAPTK_DATA_DIR + "/datasets/ExampleSmall/Subject_2D"
                                              "small2D.nii.gz");
        reader4->SetFileName(CAPTK_DATA_DIR + "/datasets/ExampleSmall/Subject_3D"
                                              "small3D.nii.gz");

        reader5->SetFileName(CAPTK_DATA_DIR + "/datasets/ExampleSmall/Subject_2D"
                                              "small2D-labels.nii.gz");
        reader6->SetFileName(CAPTK_DATA_DIR + "/datasets/ExampleSmall/Subject_2D"
                                              "small2D-labels.nii.gz");
        reader7->SetFileName(CAPTK_DATA_DIR + "/datasets/ExampleSmall/Subject_2D"
                                              "small2D-labels.nii.gz");
        reader8->SetFileName(CAPTK_DATA_DIR + "/datasets/ExampleSmall/Subject_3D"
                                              "small3D-labels.nii.gz");

        reader9->SetFileName(CAPTK_DATA_DIR + "/test-data/Common/AdaptiveGeodesicDistance/"
                                              "small2D-agd-result-lb2-limit255-imageasmask.nii.gz");
        reader10->SetFileName(CAPTK_DATA_DIR + "/test-data/Common/AdaptiveGeodesicDistance/"
                                               "small2D-agd-result-lb2-limit255-imageasmask.nii.gz");
        reader11->SetFileName(CAPTK_DATA_DIR + "/test-data/Common/AdaptiveGeodesicDistance/"
                                               "small2D-agd-result-lb2-limit255-imageasmask.nii.gz");
        reader12->SetFileName(CAPTK_DATA_DIR + "/test-data/Common/AdaptiveGeodesicDistance/"
                                               "small3D-agd-result-lb1-limit255-imageasmask.nii.gz");
        // TODO: reader13 (also uncomment stuff below)
        // TODO: reader14 (also uncomment stuff below)

        try
        {
            reader1->Update();
            reader2->Update();
            reader3->Update();
            reader4->Update();
            reader5->Update();
            reader6->Update();
            reader7->Update();
            reader8->Update();
            reader9->Update();
            reader10->Update();
            reader11->Update();
            reader12->Update();
            //reader13->Update();
            //reader14->Update();
        }
        catch (const itk::ExceptionObject &err)
        {
            CPPUNIT_FAIL(std::string(
                             "Reading on of the input images failed. "
                             "This might be because you are on an older version of the CaPTkData. "
                             "If that is the case, please rebuild your superbuild. Exception was: ") +
                         err.what());
        }

        inputImageFloat2D = reader1->GetOutput();
        inputImageInt2D = reader2->GetOutput();
        inputImageFloat2Das3D = reader3->GetOutput();
        inputImageFloat3D = reader4->GetOutput();

        inputLabelsInt2D = reader5->GetOutput();
        inputLabelsFloat2D = reader6->GetOutput();
        inputLabelsFloat2Das3D = reader7->GetOutput();
        inputLabelsInt3D = reader8->GetOutput();

        outputImageFloat2D = reader9->GetOutput();
        outputImageInt2D = reader10->GetOutput();
        outputImageInt2Das3D = reader11->GetOutput();
        outputImageFloat3D = reader12->GetOutput();

        //maskImageFloat2D = reader13->GetOutput();
        //maskImageFloat3D = reader14->GetOutput();
    }

    void tearDown() override
    {

    }

    void Try2D()
    {
        /*---- Read image ----*/

        using ReaderType = itk::ImageFileReader<ImageTypeFloat2D>;
        ReaderType::Pointer reader = ReaderType::New();
        reader->SetFileName(CAPTK_DATA_DIR 
            + "/test-data/Common/AdaptiveGeodesicDistance/"
            "small2D-agd-result-lb2-limit255-imageasmask.nii.gz");
        
        try
        {
            reader->Update();
        }
        catch (const itk::ExceptionObject &err)
        {
            CPPUNIT_FAIL(std::string(
                         "Reading the input image failed. "
                         "This might be because you are on an older version of the CaPTkData. "
                         "If that is the case, please rebuild your superbuild. Exception was: ") + 
                     err.what());
        }

        /*---- Run filter ----*/

        using FilterType = 
            itk::captk::AdaptiveGeodesicDistanceFilter< ImageTypeFloat2D, ImageTypeFloat2D >;
        FilterType::Pointer filter = FilterType::New();

        filter->SetInput( reader->GetOutput() );
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

        // CPPUNIT_ASSERT_MESSAGE(
        //     "Always passes",
        //     true);

        //... mitk::Equal() ??

        // CPPUNIT_ASSERT_MESSAGE("Should be equal", ref->Equals(outFib));
    }
};

MITK_TEST_SUITE_REGISTRATION(CaPTkItkAdaptiveGeodesicDistanceFilter)
