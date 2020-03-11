namespace itk
{
namespace captk
{

template <typename TInputImage, typename TOutputImage>
AdaptiveGeodesicDistanceFilter<TInputImage, TOutputImage>::
AdaptiveGeodesicDistanceFilter()
{
    // constructor
}

template <typename TInputImage, typename TOutputImage>
void 
AdaptiveGeodesicDistanceFilter<TInputImage, TOutputImage>::
GenerateData()
{
    // Instructions:
    // https://itk.org/Doxygen/html/Examples_2Filtering_2CompositeFilterExample_8cxx-example.html
    // https://itk.org/ITKSoftwareGuide/html/Book1/ITKSoftwareGuide-Book1ch8.html
    
    // Reset effective mask. Makes the filter reusable.
    if (m_EffectiveMask)
    {
        m_EffectiveMask = NULL;
    }

    /*---- Template dimensionality checks ----*/
    
    if (TInputImage::ImageDimension != 2 && TInputImage::ImageDimension != 3)
    {
        throw itk::ExceptionObject(__FILE__, __LINE__, 
                "Only 2D and 3D input images are supported.");
    }
    
    if (TOutputImage::ImageDimension != 2 && TOutputImage::ImageDimension != 3)
    {
        throw itk::ExceptionObject(__FILE__, __LINE__, 
                "Only 2D and 3D output images are supported.");
    }

    if (TInputImage::ImageDimension != TOutputImage::ImageDimension)
    {
        throw itk::ExceptionObject(__FILE__, __LINE__, 
                "Input and output images should have the same dimensionality.");     
    }

    /*---- Check that the input image is set ----*/

    typename ImageType::Pointer input = ImageType::New();
    input->Graft(const_cast<ImageType *>(this->GetInput()));

    if (!input)
    {
        throw itk::ExceptionObject(__FILE__, __LINE__, 
                "No input image.");
    }

    /*---- Check that label of interest is set ----*/

    if (m_LabelOfInterest == 0)
    {
        throw itk::ExceptionObject(__FILE__, __LINE__, 
                "Label of interest is not set (or is set to zero).");
    }

    /*---- Check that labels image is set ----*/

    if (!m_Labels)
    {
        throw itk::ExceptionObject(__FILE__, __LINE__,
                                   "Labels image is not set.");
    }

    /*---- Create mask if necessary ----*/

    if (m_Mask)
    {
        m_EffectiveMask = m_Mask;
    }
    else
    {
        if (m_UseInputImageAsMask)
        {
            m_EffectiveMask = input;
        }
        else
        {
            // Do calculations for the whole area.
            // We create a mask image that has 1 everywhere
            m_EffectiveMask = ImageType::New();
            m_EffectiveMask->CopyInformation( input );
            m_EffectiveMask->SetRequestedRegion( input->GetLargestPossibleRegion() );
            m_EffectiveMask->SetBufferedRegion( input->GetBufferedRegion() );
            m_EffectiveMask->Allocate();
            m_EffectiveMask->FillBuffer(1);
        }
    }

    /*---- Check that input image, labels image, 
           and mask have the same information ----*/

    // Labels is int, but input can be other types too
    // To compare meta-information we create a "template" labels image
    // that holds the same meta-information as the labels image
    typename ImageType::Pointer labelsJustInformation = ImageType::New();
    labelsJustInformation->CopyInformation( m_Labels );
    labelsJustInformation->SetRequestedRegion( m_Labels->GetLargestPossibleRegion() );
    labelsJustInformation->SetBufferedRegion( m_Labels->GetBufferedRegion() );
    labelsJustInformation->Allocate();
    labelsJustInformation->FillBuffer(1);

    using DiffType = itk::Testing::ComparisonImageFilter<ImageType, ImageType>;
    typename DiffType::Pointer diff = DiffType::New();
    diff->VerifyInputInformationOn();
    diff->SetValidInput( input );

    // Check input image and & "template" labels
    diff->SetTestInput( labelsJustInformation );
    try
    {
        diff->UpdateLargestPossibleRegion();
    }
    catch(const itk::ExceptionObject &err)
    {
        throw itk::ExceptionObject(__FILE__, __LINE__, 
                "Input image and labels image information mismatch.");        
    }

    // memory optimization
    labelsJustInformation = NULL; // ref: https://itk.org/Wiki/ITK/Tutorials/DOs_and_DONTs

    // Check input image & effective mask
    diff->SetTestInput( m_EffectiveMask );
    try
    {
        diff->UpdateLargestPossibleRegion();
    }
    catch(const itk::ExceptionObject &err)
    {
        throw itk::ExceptionObject(__FILE__, __LINE__, 
                "Input image and mask image information mismatch.");    
    }

    // Check that spacing is isotropic
    double spacingDim0 = input->GetSpacing()[0];
    for (unsigned int i = 1; i < ImageType::ImageDimension; i++)
    {
        if (input->GetSpacing()[i] - spacingDim0 > 0.05)
        {
            throw itk::ExceptionObject(__FILE__, __LINE__, 
                    "AGD filter needs isotropic spacing.");
        }
    } 
    
    /*---- Actual calculations ----*/
    
    this->InternalProcessing();
}

template <typename TInputImage, typename TOutputImage>
void
AdaptiveGeodesicDistanceFilter<TInputImage, TOutputImage>::
InternalProcessing()
{
    typename ImageType::Pointer input = ImageType::New();
    input->Graft(const_cast<ImageType *>(this->GetInput()));

    using Iterator          = itk::ImageRegionIterator<ImageType>;
	using IteratorLabels    = itk::ImageRegionIterator<LabelsImageType>;
	using IdxIterator       = itk::ImageRegionIteratorWithIndex<ImageType>;
    using NIterator         = itk::NeighborhoodIterator<ImageType>;

    /*---- Allocate output image ----*/

    typename ImageType::Pointer output = ImageType::New();
    output->SetRegions( input->GetLargestPossibleRegion() );
    output->SetRequestedRegion( input->GetLargestPossibleRegion() );
    output->SetBufferedRegion( input->GetBufferedRegion() );
    output->Allocate();
    output->FillBuffer((m_LimitAt255) ? 
            255 : 
            static_cast<typename ImageType::PixelType>(
                itk::NumericTraits< typename ImageType::PixelType >::max() - 1
            )
    ); // 255 or (maximum_possible_value - 1)
    output->SetDirection(input->GetDirection());
    output->SetOrigin(input->GetOrigin());
    output->SetSpacing(input->GetSpacing());

    /*---- Initialize values for output image 
           (0 at pixels of label of interest, max elsewhere [set above]) ----*/

    Iterator       outIter(output,   output->GetLargestPossibleRegion());
    IteratorLabels labIter(m_Labels, m_Labels->GetLargestPossibleRegion());
    outIter.GoToBegin();
    labIter.GoToBegin();

    while (!labIter.IsAtEnd()) 
    {
        if (labIter.Get() == m_LabelOfInterest) 
        {
            outIter.Set(0);
        }
        ++labIter;
        ++outIter;
    }

    /*---- Variables needed for AGD ----*/

    // We want a radius of 1 in the neighborhood iterators
    typename ImageType::SizeType radius;
    for (unsigned int i = 0; i < ImageType::ImageDimension; i++) 
    {
        radius[i] = 1;
    }

    // Iterators that are used in the loops
    IdxIterator maskIter(m_EffectiveMask, m_EffectiveMask->GetLargestPossibleRegion());
    NIterator outNIter(radius, output, output->GetLargestPossibleRegion());
    NIterator inputNIter(radius, input, input->GetLargestPossibleRegion());

    // Helper variables that are used in the loops
    // Declaring them outside the loops makes an actual difference speed-wise
    typename ImageType::PixelType inpCenterPixel;
    double arr[14];
    double minVal;
    int ii;

    /*---- Forward pass ----*/

    outNIter.GoToBegin();
    inputNIter.GoToBegin();
    maskIter.GoToBegin();

    while (!maskIter.IsAtEnd())
    {
        if (maskIter.Get() != 0)
        {
            // Forward pass
            
            if (ImageType::ImageDimension == 2)
            {
                // 2D

                inpCenterPixel = inputNIter.GetPixel(4);

                arr[4] = outNIter.GetCenterPixel();
                arr[0] = outNIter.GetPixel(0) + sqrt(2.0 + square<ImageType>(
                            inpCenterPixel - inputNIter.GetPixel(0)));
                arr[1] = outNIter.GetPixel(1) + sqrt(1.0 + square<ImageType>(
                            inpCenterPixel - inputNIter.GetPixel(1)));
                arr[2] = outNIter.GetPixel(2) + sqrt(2.0 + square<ImageType>(
                            inpCenterPixel - inputNIter.GetPixel(2)));
                arr[3] = outNIter.GetPixel(3) + sqrt(1.0 + square<ImageType>(
                            inpCenterPixel - inputNIter.GetPixel(3)));
                
                minVal = arr[4];
                for (ii = 0; ii < 4; ii++)
                {
                    if (arr[ii] < minVal) {
                        minVal = arr[ii];
                    }
                }
            }
            else {
                // 3D
                
                inpCenterPixel = inputNIter.GetPixel(13);
                //gamPixel = iterGamma.Get();

                arr[13] = outNIter.GetCenterPixel();
                arr[0]  = outNIter.GetPixel(4)  + sqrt(1.0 + square<ImageType>(
                            inpCenterPixel - inputNIter.GetPixel(4)));
                arr[1]  = outNIter.GetPixel(10) + sqrt(1.0 + square<ImageType>(
                            inpCenterPixel - inputNIter.GetPixel(10)));
                arr[2]  = outNIter.GetPixel(12) + sqrt(1.0 + square<ImageType>(
                            inpCenterPixel - inputNIter.GetPixel(12)));
                arr[3]  = outNIter.GetPixel(1)  + sqrt(2.0 + square<ImageType>(
                            inpCenterPixel - inputNIter.GetPixel(1)));
                arr[4]  = outNIter.GetPixel(3)  + sqrt(2.0 + square<ImageType>(
                            inpCenterPixel - inputNIter.GetPixel(3)));
                arr[5]  = outNIter.GetPixel(9)  + sqrt(2.0 + square<ImageType>(
                            inpCenterPixel - inputNIter.GetPixel(9)));
                arr[6]  = outNIter.GetPixel(0)  + sqrt(3.0 + square<ImageType>(
                            inpCenterPixel - inputNIter.GetPixel(0)));
                arr[7]  = outNIter.GetPixel(7)  + sqrt(2.0 + square<ImageType>(
                            inpCenterPixel - inputNIter.GetPixel(7)));
                arr[8]  = outNIter.GetPixel(6)  + sqrt(3.0 + square<ImageType>(
                            inpCenterPixel - inputNIter.GetPixel(6)));
                arr[9]  = outNIter.GetPixel(15) + sqrt(2.0 + square<ImageType>(
                            inpCenterPixel - inputNIter.GetPixel(15)));
                arr[10] = outNIter.GetPixel(24) + sqrt(3.0 + square<ImageType>(
                            inpCenterPixel - inputNIter.GetPixel(24)));
                arr[11] = outNIter.GetPixel(21) + sqrt(2.0 + square<ImageType>(
                            inpCenterPixel - inputNIter.GetPixel(21)));
                arr[12] = outNIter.GetPixel(18) + sqrt(3.0 + square<ImageType>(
                            inpCenterPixel - inputNIter.GetPixel(18)));

                minVal = arr[13];
                for (ii = 0; ii < 13; ii++)
                {
                    if (arr[ii] < minVal) {
                        minVal = arr[ii];
                    }
                }
            }

            outNIter.SetCenterPixel(minVal);
        }
        ++outNIter;
        ++inputNIter;
        ++maskIter;
    }

    /*---- Backward pass ----*/

    outNIter.GoToEnd();
    inputNIter.GoToEnd();
    --outNIter;
    --inputNIter;
    maskIter.GoToReverseBegin();

    while (!maskIter.IsAtReverseEnd())
    {
        if (maskIter.Get() != 0)
        {
            // Backward pass

            if (ImageType::ImageDimension == 2) 
            {
                // 2D

                inpCenterPixel = inputNIter.GetPixel(4);
                //gamPixelB = iterGamma.Get();

                arr[4] = outNIter.GetCenterPixel();
                arr[0] = outNIter.GetPixel(5) + sqrt(1.0 + square<ImageType>(
                            inpCenterPixel - inputNIter.GetPixel(5)));
                arr[1] = outNIter.GetPixel(6) + sqrt(2.0 + square<ImageType>(
                            inpCenterPixel - inputNIter.GetPixel(6)));
                arr[2] = outNIter.GetPixel(7) + sqrt(1.0 + square<ImageType>(
                            inpCenterPixel - inputNIter.GetPixel(7)));
                arr[3] = outNIter.GetPixel(8) + sqrt(2.0 + square<ImageType>(
                            inpCenterPixel - inputNIter.GetPixel(8)));
                
                minVal = arr[4];
                for (int i = 0; i < 4; i++)
                {
                    if (arr[i] < minVal) {
                        minVal = arr[i];
                    }
                }
            }
            else
            {
                // 3D

                inpCenterPixel = inputNIter.GetPixel(13);
                //gamPixelB = iterGamma.Get();

                arr[13] = outNIter.GetCenterPixel();
                arr[0] = outNIter.GetPixel(22) + sqrt(1.0 + square<ImageType>(
                            inpCenterPixel - inputNIter.GetPixel(22)));
                arr[1] = outNIter.GetPixel(16) + sqrt(1.0 + square<ImageType>(
                            inpCenterPixel - inputNIter.GetPixel(16)));
                arr[2] = outNIter.GetPixel(14) + sqrt(1.0 + square<ImageType>(
                            inpCenterPixel - inputNIter.GetPixel(14)));
                arr[3] = outNIter.GetPixel(25) + sqrt(2.0 + square<ImageType>(
                            inpCenterPixel - inputNIter.GetPixel(25)));
                arr[4] = outNIter.GetPixel(23) + sqrt(2.0 + square<ImageType>(
                            inpCenterPixel - inputNIter.GetPixel(23)));
                arr[5] = outNIter.GetPixel(17) + sqrt(2.0 + square<ImageType>(
                            inpCenterPixel - inputNIter.GetPixel(17)));
                arr[6] = outNIter.GetPixel(26) + sqrt(3.0 + square<ImageType>(
                            inpCenterPixel - inputNIter.GetPixel(26)));
                arr[7] = outNIter.GetPixel(19) + sqrt(2.0 + square<ImageType>(
                            inpCenterPixel - inputNIter.GetPixel(19)));
                arr[8] = outNIter.GetPixel(20) + sqrt(3.0 + square<ImageType>(
                            inpCenterPixel - inputNIter.GetPixel(20)));
                arr[9] = outNIter.GetPixel(11) + sqrt(2.0 + square<ImageType>(
                            inpCenterPixel - inputNIter.GetPixel(11)));
                arr[10] = outNIter.GetPixel(2) + sqrt(3.0 + square<ImageType>(
                            inpCenterPixel - inputNIter.GetPixel(2)));
                arr[11] = outNIter.GetPixel(5) + sqrt(2.0 + square<ImageType>(
                            inpCenterPixel - inputNIter.GetPixel(5)));
                arr[12] = outNIter.GetPixel(8) + sqrt(3.0 + square<ImageType>(
                            inpCenterPixel - inputNIter.GetPixel(8)));
            
                minVal = arr[13];
                for (int i = 0; i < 13; i++)
                {
                    if (arr[i] < minVal) {
                        minVal = arr[i];
                    }
                }
            }

            outNIter.SetCenterPixel(minVal);
        }
        --outNIter;
        --inputNIter;
        --maskIter;
    }

    this->GraftOutput(output);
}

}
}