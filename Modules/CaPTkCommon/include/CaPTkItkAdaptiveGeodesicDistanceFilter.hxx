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
    
    /*---- Actual calculations ----*/
    
    this->InternalProcessing();
}

template <typename TInputImage, typename TOutputImage>
void
AdaptiveGeodesicDistanceFilter<TInputImage, TOutputImage>::
InternalProcessing()
{
    // TODO
    // Delete this:
    this->GraftOutput(m_EffectiveMask);
}

}
}