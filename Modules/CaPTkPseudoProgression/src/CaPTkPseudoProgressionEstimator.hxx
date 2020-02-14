template <class ImageType>
void PseudoProgressionEstimator::PseudoProgressionEstimateOnGivenSubject(typename ImageType::Pointer edemaMask,
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
                                                                         const VectorVectorDouble &farRegionIndices)
{
}

template <class ImageType>
void PseudoProgressionEstimator::PseudoProgressionEstimateOnGivenSubjectUsingExistingModel(typename ImageType::Pointer edemaMask,
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
                                                                                           const VectorVectorDouble &farRegionIndices, const std::string modeldirectory)
{
    cbica::Logging(loggerFile, "mCurrentOutputDir::" + mCurrentOutputDir);
    mOutputLocalPtr.SetOutputDirectoryPath(mCurrentOutputDir);
    VariableSizeMatrixType pca_coefficients;
    VariableLengthVectorType pca_mean;
    VariableLengthVectorType mean;
    VariableLengthVectorType stds;

    //check for the presence of model file
    if (!cbica::fileExists(modeldirectory + "/" + mPseudoTrainedFile) || !cbica::fileExists(modeldirectory + "/" + mRecurrenceTrainedFile))
    {
        mLastEncounteredError = "SVM model file is not present in the directory: " + modeldirectory;
        return;
    }

    //check for the presence of z-score record
    if (!cbica::fileExists(modeldirectory + "/Recurrence_ZScore_Mean.csv") || !cbica::fileExists(modeldirectory + "/Recurrence_ZScore_Std.csv"))
    {
        mLastEncounteredError = "Z-score record is not present in the directory: " + modeldirectory;
        return;
    }

    if (!cbica::fileExists(modeldirectory + "/Recurrence_COEF.csv") || !cbica::fileExists(modeldirectory + "/Recurrence_MR.csv"))
    {
        mLastEncounteredError = "PCA parameters are not present in the directory: " + modeldirectory;
        return;
    }
    else
    {
        mOutputLocalPtr.ReadModelParameters(modeldirectory + "/Recurrence_ZScore_Mean.csv", modeldirectory + "/Recurrence_ZScore_Std.csv", modeldirectory + "/Recurrence_COEF.csv", modeldirectory + "/Recurrence_MR.csv", mean, stds, pca_coefficients, pca_mean);
        mFeatureReductionLocalPtr.SetParameters(pca_coefficients, pca_mean);
        mFeatureScalingLocalPtr.SetParameters(mean, stds);
    }

    cbica::Logging(loggerFile, "Model reading finished.");

    ////------------------------------------------testing data formulation---------------------------------------------------

    VectorVectorDouble perfusionIntensities;
    VectorVectorDouble otherIntensities;
    VectorDouble distanceIntensities;
    std::vector<typename ImageType::IndexType> testindices;

    testindices = mNiftiLocalPtr.LoadTestData(FinalT1CEImagePointer, FinalT2FlairImagePointer, FinalT1ImagePointer, FinalT2ImagePointer, FinalPerfusionImagePointer, FinalDTIImagePointer, tumorMask, edemaMask, perfusionIntensities, otherIntensities, distanceIntensities, CAPTK::ImageExtension::DICOM, conDataPresent, dtiDataPresent, perfusionDataPresent, distanceDataPresent);

    if (testindices.empty())
    {
        cbica::Logging(loggerFile, "{RecurrDebug}testindices is empty");
    }

    VectorVectorDouble reducedTestPerfusionFeatures;
    reducedTestPerfusionFeatures = mFeatureReductionLocalPtr.ApplyPCAOnTestData(perfusionIntensities);

    if (reducedTestPerfusionFeatures.empty())
    {
        cbica::Logging(loggerFile, "{RecurrDebug}reducedTestPerfusionFeatures is empty");
    }

    int NumberOfPCs = 5;
    VectorVectorDouble globaltestintensities;

    for (unsigned int k = 0; k < testindices.size(); k++)
    {
        VectorDouble inten;

        for (int j = 0; j < NumberOfPCs; j++)
            inten.push_back(reducedTestPerfusionFeatures[k][j]);

        for (unsigned int j = 0; j < otherIntensities[0].size(); j++)
            inten.push_back(otherIntensities[k][j]);

        inten.push_back(distanceIntensities[k]);

        if (!inten.empty())
            globaltestintensities.push_back(inten);
    }

    if (globaltestintensities.empty())
    {
        cbica::Logging(loggerFile, "{RecurrDebug}globaltestintensities is empty");
    }

    VariableSizeMatrixType TestingData = mFeatureExtractionLocalPtr.FormulateTestData(globaltestintensities);

    VariableSizeMatrixType ScaledTestingData = mFeatureScalingLocalPtr.ScaleGivenTestingFeatures(TestingData);

    if (TestingData.Cols() == 0)
    {
        cbica::Logging(loggerFile, "{RecurrDebug}TestingData.Cols is empty");
    }

    if (ScaledTestingData.Cols() == 0)
    {
        cbica::Logging(loggerFile, "{RecurrDebug}ScaledTestingData.Cols is empty");
    }

#ifdef APP_BASE_CAPTK_H
    progressUpdate(60);
#endif

    typename ImageType::RegionType region = FinalT1CEImagePointer->GetLargestPossibleRegion();
    typename ImageType::Pointer RecProbabilityMap = ImageType::New();
    RecProbabilityMap->SetRegions(region);
    RecProbabilityMap->Allocate();
    RecProbabilityMap->SetSpacing(FinalT1CEImagePointer->GetSpacing());
    RecProbabilityMap->SetOrigin(FinalT1CEImagePointer->GetOrigin());
    RecProbabilityMap->SetDirection(FinalT1CEImagePointer->GetDirection());

    typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;
    IteratorType imIt(edemaMask, edemaMask->GetLargestPossibleRegion());
    IteratorType RecIt(RecProbabilityMap, RecProbabilityMap->GetLargestPossibleRegion());
    imIt.GoToBegin();
    RecIt.GoToBegin();

    while (!imIt.IsAtEnd())
    {
        if (imIt.Get() == 0)
            RecIt.Set(0);

        ++imIt;
        ++RecIt;
    }
    cbica::Logging(loggerFile, "Before testing.");
    VectorDouble result_modified;

    try
    {
        if (cbica::fileExists(modeldirectory + "/" + mPseudoTrainedFile))
        {
            //cbica::Logging(loggerFile, "Before testing 1.");
            VariableLengthVectorType result;
            result = DistanceFunction(ScaledTestingData, modeldirectory + "/" + mPseudoTrainedFile, RECURRENCE_MODEL_RHO, RECURRENCE_MODEL_G);
            for (unsigned int index = 0; index < result.Size(); index++)
            {
                RecProbabilityMap->SetPixel(testindices[index], result[index] * 1);
                result_modified.push_back(result[index]);
            }
        }
        else
        {
            //cbica::Logging(loggerFile, "Before testing 2.");
            VectorDouble result;
            result = testOpenCVSVM(ScaledTestingData, modeldirectory + "/" + mPseudoTrainedFile);
            for (unsigned int index = 0; index < result.size(); index++)
                RecProbabilityMap->SetPixel(testindices[index], result[index] * 1);
            result_modified = result;
        }
    }
    catch (itk::ExceptionObject &excp)
    {
        cbica::Logging(loggerFile, "Error caught during testing: " + std::string(excp.GetDescription()));
        exit(EXIT_FAILURE);
    }

    VectorDouble result_revised = RecurrenceMapPostprocessing<ImageType>(result_modified, testindices, RecProbabilityMap, edemaMask);
    for (unsigned int index = 0; index < result_modified.size(); index++)
        RecProbabilityMap->SetPixel(testindices[index], result_revised[index] * 1);

    //averaging filter
    typedef itk::MeanImageFilter<ImageType, ImageType> FilterType;
    typename FilterType::Pointer meanFilter = FilterType::New();
    typename FilterType::InputSizeType radius;
    radius.Fill(1);
    meanFilter->SetRadius(radius);
    meanFilter->SetInput(RecProbabilityMap);
    typename ImageType::Pointer RevisedRecurrenceMap = meanFilter->GetOutput();

    //------------------------------------------Writing final output--------------------------------------------------
    mRecurrenceMapFileName = "RecurrenceMap.nii.gz";
    mOutputLocalPtr.WriteRecurrenceOutput<ImageType>(RevisedRecurrenceMap, t1cebasefilename, mRecurrenceMapFileName);

#ifdef APP_BASE_CAPTK_H
    progressUpdate(90);
#endif

#ifdef APP_BASE_CAPTK_H
    progressUpdate(100);
#endif
}

template <class ImageType>
typename ImageType::Pointer PseudoProgressionEstimator::RescaleImageIntensity(const typename ImageType::Pointer &image)
{
    typedef itk::RescaleIntensityImageFilter<ImageType, ImageType> RescaleFilterType;
    typename RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
    rescaleFilter->SetInput(image);
    rescaleFilter->SetOutputMinimum(0);
    rescaleFilter->SetOutputMaximum(255);
    rescaleFilter->Update();
    typename ImageType::Pointer outputimage = rescaleFilter->GetOutput();
    return outputimage;
}

template <class ImageType>
typename ImageType::Pointer PseudoProgressionEstimator::ReadNiftiImage(const std::string &filename)
{
    typedef itk::ImageFileReader<ImageType> ImageReaderType;
    typename ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(filename);

    try
    {
        reader->Update();
    }
    catch (itk::ExceptionObject &e)
    {
        std::cerr << "Error caught: " << e.what() << "\n";
        exit(EXIT_FAILURE);
    }

    return reader->GetOutput();
}

template <class ImageType>
VectorDouble PseudoProgressionEstimator::RecurrenceMapPostprocessing(VectorDouble result, std::vector<typename ImageType::IndexType> testindices, typename ImageType::Pointer RecurrenceProbabilityMap, typename ImageType::Pointer edemaMap)
{
    VectorDouble revised_result;
    VectorDouble positiveDistances;
    VectorDouble negativeDistances;
    for (size_t x = 0; x < result.size(); x++)
    {
        if (result[x] > 0)
            positiveDistances.push_back(result[x]);
        else
            negativeDistances.push_back(result[x]);
    }
    std::sort(positiveDistances.begin(), positiveDistances.end());
    std::sort(negativeDistances.rbegin(), negativeDistances.rend());

    double positivePercentileCutOff = positiveDistances[positiveDistances.size() - 1];
    double negativePercentileCutOff = negativeDistances[negativeDistances.size() - 1];

    for (size_t x = 0; x < positiveDistances.size(); x++)
    {
        double percentile = (100 * (x + 0.5)) / positiveDistances.size();
        if (percentile > 99)
        {
            positivePercentileCutOff = positiveDistances[x];
            break;
        }
    }
    for (size_t x = 0; x < negativeDistances.size(); x++)
    {
        double percentile = (100 * (x + 0.5)) / negativeDistances.size();
        if (percentile > 99)
        {
            negativePercentileCutOff = negativeDistances[x];
            break;
        }
    }
    for (size_t x = 0; x < result.size(); x++)
    {
        if (result[x] > positivePercentileCutOff)
            result[x] = positivePercentileCutOff;
        if (result[x] < negativePercentileCutOff)
            result[x] = negativePercentileCutOff;
    }
    double min = 0;
    double max = positivePercentileCutOff;
    for (size_t x = 0; x < result.size(); x++)
    {
        if (result[x] > 0)
            result[x] = (result[x] - min) / (max - min);
    }
    min = negativePercentileCutOff;
    max = 0;
    for (size_t x = 0; x < result.size(); x++)
    {
        if (result[x] < 0)
            result[x] = (result[x] - min) / (max - min);
    }
    for (int x = 0; x < result[x]; x++)
        result[x] = (result[x] + 1) / 2;

    revised_result = result;
    return revised_result;
}

template <class PerfusionImageType, class ImageType>
VariableSizeMatrixType PseudoProgressionEstimator::LoadPerfusionData(typename ImageType::Pointer maskImagePointerNifti, typename PerfusionImageType::Pointer perfImagePointerNifti, std::vector<typename ImageType::IndexType> &qualifiedIndices)
{
    VectorVectorDouble pIntensities;
    //----------------------find the voxels of the mask------------------------------
    typedef itk::ImageRegionIteratorWithIndex<ImageTypeFloat3D> IteratorType;
    IteratorType maskIt(maskImagePointerNifti, maskImagePointerNifti->GetLargestPossibleRegion());

    maskIt.GoToBegin();
    int mask_counter = 0;
    while (!maskIt.IsAtEnd())
    {
        if (maskIt.Get() == 1)
        {
            mask_counter++;
            qualifiedIndices.push_back(maskIt.GetIndex());
        }
        ++maskIt;
    }
    //--------------------------populate the covariance matrix --------------------------------
    typename ImageTypeFloat4D::RegionType region = perfImagePointerNifti->GetLargestPossibleRegion();
    VariableSizeMatrixType revisedcovariancematrix;
    revisedcovariancematrix.SetSize(mask_counter, region.GetSize()[3]);
    for (size_t i = 0; i < qualifiedIndices.size(); i++)
    {
        typename ImageTypeFloat4D::IndexType regionIndex;
        regionIndex[0] = qualifiedIndices[i][0];
        regionIndex[1] = qualifiedIndices[i][1];
        regionIndex[2] = qualifiedIndices[i][2];
        for (size_t j = 0; j < region.GetSize()[3]; j++)
        {
            regionIndex[3] = j;
            revisedcovariancematrix(i, j) = perfImagePointerNifti->GetPixel(regionIndex);
        }
    }
    return revisedcovariancematrix;
}

template <class ImageType>
std::tuple<VectorDouble, VectorDouble, VectorDouble, VectorDouble, VectorDouble> PseudoProgressionEstimator::GetAllFeaturesPerImagePerROI(typename ImageType::Pointer image, typename ImageType::Pointer mask, std::string modality)
{
    std::vector<typename ImageType::IndexType> roiIndices;
    typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;
    IteratorType imIt(mask, mask->GetLargestPossibleRegion());
    imIt.GoToBegin();
    while (!imIt.IsAtEnd())
    {
        if (imIt.Get() == 1)
            roiIndices.push_back(imIt.GetIndex());
        ++imIt;
    }

    typedef std::vector<float> VectorFloat;
    VectorFloat ROIIntensities;
    for (unsigned int i = 0; i < roiIndices.size(); i++)
        ROIIntensities.push_back(std::round(image.GetPointer()->GetPixel(roiIndices[i])));

    double minvalue = *min_element(ROIIntensities.begin(), ROIIntensities.end());
    double maxvalue = *max_element(ROIIntensities.begin(), ROIIntensities.end());

    VectorDouble HistogramFeatures = GetHistogramFeatures(ROIIntensities, 10);
    VectorDouble HistogramFeatures1 = GetHistogramFeatures(ROIIntensities, 20);
    VectorDouble IntensityFeatures = GetIntensityFeatures(ROIIntensities);
    typename ImageType::Pointer isotropicImageGLCM = cbica::ResampleImage<ImageType>(image, 1.0, "Linear");
    typename ImageType::Pointer isotropicMaskGLCM = cbica::ResampleImage<ImageType>(mask, 1.0, "Nearest");
    auto roundingFilter = itk::RoundImageFilter<ImageType, ImageType>::New();
    roundingFilter->SetInput(isotropicMaskGLCM);
    roundingFilter->Update();
    isotropicMaskGLCM = roundingFilter->GetOutput();
    VectorDouble GLCMFeatures = GetGLCMFeatures<ImageType>(isotropicImageGLCM, isotropicMaskGLCM, minvalue, maxvalue);

    typename ImageType::Pointer isotropicImageGLRLM = cbica::ResampleImage<ImageType>(image, 1.0, "Linear");
    typename ImageType::Pointer isotropicMaskGLRLM = cbica::ResampleImage<ImageType>(mask, 1.0, "Nearest");
    roundingFilter->SetInput(isotropicMaskGLRLM);
    roundingFilter->Update();
    isotropicMaskGLRLM = roundingFilter->GetOutput();
    VectorDouble GLRLMFeatures = GetRunLengthFeatures<ImageType>(isotropicImageGLRLM, isotropicMaskGLRLM, minvalue, maxvalue);

    std::tuple<VectorDouble, VectorDouble, VectorDouble, VectorDouble, VectorDouble> new_tuple(HistogramFeatures, IntensityFeatures, GLCMFeatures, GLRLMFeatures, HistogramFeatures1);
    return new_tuple;
}

template <typename TImageTypeShape>
VectorDouble PseudoProgressionEstimator::GetShapeFeatures(typename TImageTypeShape::Pointer mask)
{
    std::vector<double> features;

    typedef short LabelType;
    typedef itk::Image<LabelType, TImageTypeShape::ImageDimension> OutputImageType;
    typedef itk::ShapeLabelObject<LabelType, TImageTypeShape::ImageDimension> ShapeLabelObjectType;
    typedef itk::LabelMap<ShapeLabelObjectType> LabelMapType;

    typedef itk::ConnectedComponentImageFilter<TImageTypeShape, OutputImageType> ConnectedComponentImageFilterType;
    typedef itk::LabelImageToShapeLabelMapFilter<OutputImageType, LabelMapType> I2LType;
    typename ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
    connected->SetInput(mask);
    connected->Update();

    typename I2LType::Pointer i2l = I2LType::New();
    i2l->SetInput(connected->GetOutput());
    i2l->SetComputePerimeter(true);
    i2l->Update();
    typename LabelMapType::Pointer labelMap = i2l->GetOutput();
    auto spacing = mask->GetSpacing();
    double voxvol = spacing[0] * spacing[1] * spacing[2];

    int numbers = labelMap->GetNumberOfLabelObjects();
    std::vector<double> eccentricity1;
    std::vector<double> eccentricity2;
    std::vector<double> roundness;
    std::vector<double> flatness;
    std::vector<double> elongation;
    std::vector<double> perimeter;

    for (int i = 1; i < numbers; i++)
    {
        typename ShapeLabelObjectType::Pointer labelObject = labelMap->GetNthLabelObject(i);

        if (labelObject->GetNumberOfPixels() < 20)
            continue;
        auto princComps = labelObject->GetPrincipalMoments();
        eccentricity1.push_back(sqrt(1 - std::pow((princComps[0] / princComps[2]), 2)));
        eccentricity2.push_back(sqrt(1 - std::pow((princComps[0] / princComps[1]), 2)));
        roundness.push_back(labelObject->GetRoundness());
        flatness.push_back(labelObject->GetFlatness());
        elongation.push_back(labelObject->GetElongation());
        perimeter.push_back(labelObject->GetPerimeter());
    }
    double mean_ecc1 = 0;
    double mean_ecc2 = 0;
    double mean_round = 0;
    double mean_flat = 0;
    double mean_perim = 0;
    double mean_elong = 0;
    for (unsigned int i = 0; i < eccentricity1.size(); i++)
    {
        mean_ecc1 = mean_ecc1 + eccentricity1[i];
        mean_ecc2 = mean_ecc2 + eccentricity2[i];
        mean_round = mean_round + roundness[i];
        mean_flat = mean_flat + flatness[i];
        mean_perim = mean_perim + perimeter[i];
        mean_elong = mean_elong + elongation[i];
    }

    features.push_back(mean_ecc1 / eccentricity1.size());
    features.push_back(mean_elong / elongation.size());
    features.push_back(mean_perim);
    features.push_back(mean_round / roundness.size());
    features.push_back(mean_flat / flatness.size());

    return features;
}

template <class ImageType>
VectorDouble PseudoProgressionEstimator::GetGLCMFeatures(typename ImageType::Pointer image, typename ImageType::Pointer mask, double minvalue, double maxvalue)
{
    double m_Bins = 16;
    using FeatureextractionImageType = typename ImageType::Pointer;
    using Image2CoOccuranceType = itk::Statistics::ScalarImageToCooccurrenceMatrixFilter<ImageType>;
    using HistogramType = typename Image2CoOccuranceType::HistogramType;
    using Hist2FeaturesType = itk::Statistics::HistogramToTextureFeaturesFilter<HistogramType>;
    using OffsetType = typename ImageType::OffsetType;
    using Offsets = OffsetType;
    using OffsetVector = itk::VectorContainer<unsigned char, OffsetType>;

    double inputRadius = 1;
    double inputDirections = 13;
    itk::Neighborhood<typename ImageType::PixelType, ImageType::ImageDimension> neighborhood;
    neighborhood.SetRadius(inputRadius);
    auto size = neighborhood.GetSize();
    auto directionsToCompute = 1;
    for (size_t sizeCount = 0; sizeCount < ImageType::ImageDimension; sizeCount++)
    {
        directionsToCompute *= size[sizeCount];
    }
    if (inputDirections < directionsToCompute)
    {
        directionsToCompute = inputDirections;
    }

    typename OffsetVector::Pointer offsets = OffsetVector::New();
    auto centerIndex = neighborhood.GetCenterNeighborhoodIndex();

    for (int d = directionsToCompute - 1; d >= 0; d--)
    {
        if (d != static_cast<int>(centerIndex))
        {
            offsets->push_back(neighborhood.GetOffset(d));
        }
    }

    std::vector<double> features;
    double contrast = 0, correl = 0, ener = 0, entro = 0, homo = 0, clustershade = 0, clusterprominance = 0, autocorr = 0;

    for (size_t i = 0; i < offsets->size(); i++)
    {
        auto glcmGenerator = Image2CoOccuranceType::New();
        glcmGenerator->SetNumberOfBinsPerAxis(m_Bins); //reasonable number of bins
        glcmGenerator->SetPixelValueMinMax(minvalue, maxvalue);
        glcmGenerator->SetMaskImage(mask);
        glcmGenerator->SetInput(image);
        auto featureCalc = Hist2FeaturesType::New();

        glcmGenerator->SetOffset(offsets->at(i));
        glcmGenerator->Update();
        featureCalc->SetInput(glcmGenerator->GetOutput());
        featureCalc->Update();

        contrast += static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::Inertia));
        correl += static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::Correlation));
        ener += static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::Energy));
        homo += static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::InverseDifferenceMoment));
        entro += static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::Entropy));
        clustershade += static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::ClusterShade));
        clusterprominance += static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::ClusterProminence));
        autocorr += static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::HaralickCorrelation));
    }

    contrast = contrast / offsets->size();
    correl = correl / offsets->size();
    ener = ener / offsets->size();
    homo = homo / offsets->size();
    entro = entro / offsets->size();
    clusterprominance = clusterprominance / offsets->size();
    clustershade = clustershade / offsets->size();
    autocorr = autocorr / offsets->size();

    features.push_back(correl);
    features.push_back(contrast);
    features.push_back(entro);
    features.push_back(homo);
    features.push_back(clustershade);
    features.push_back(clusterprominance);
    features.push_back(autocorr);
    features.push_back(ener);
    return features;
}

template <class ImageType>
VectorDouble PseudoProgressionEstimator::GetRunLengthFeatures(typename ImageType::Pointer image, typename ImageType::Pointer mask, double minvalue, double maxvalue)
{
    double m_Bins = 16;
    using HistogramFrequencyContainerType = itk::Statistics::DenseFrequencyContainer2;
    using RunLengthFilterType = itk::Statistics::EnhancedScalarImageToRunLengthFeaturesFilter<ImageType, HistogramFrequencyContainerType>;
    using RunLengthMatrixGenerator = typename RunLengthFilterType::RunLengthMatrixFilterType;
    using RunLengthFeatures = typename RunLengthFilterType::RunLengthFeaturesFilterType;
    using OffsetType = typename ImageType::OffsetType;
    using OffsetVector = itk::VectorContainer<unsigned char, OffsetType>;

    //offset calculation
    double inputDirections = 13;
    itk::Neighborhood<typename ImageType::PixelType, ImageType::ImageDimension> neighborhood;
    neighborhood.SetRadius(1);
    auto size = neighborhood.GetSize();
    auto directionsToCompute = 1;
    for (size_t sizeCount = 0; sizeCount < ImageType::ImageDimension; sizeCount++)
    {
        directionsToCompute *= size[sizeCount];
    }
    if (inputDirections < directionsToCompute)
    {
        directionsToCompute = inputDirections;
    }
    typename OffsetVector::Pointer offsets = OffsetVector::New();
    auto centerIndex = neighborhood.GetCenterNeighborhoodIndex();
    for (int d = directionsToCompute - 1; d >= 0; d--)
    {
        if (d != static_cast<int>(centerIndex))
        {
            offsets->push_back(neighborhood.GetOffset(d));
        }
    }

    //matrix generation
    typename RunLengthMatrixGenerator::Pointer matrix_generator = RunLengthMatrixGenerator::New();
    matrix_generator->SetInput(image);
    matrix_generator->SetMaskImage(mask);
    matrix_generator->SetInsidePixelValue(1);
    matrix_generator->SetPixelValueMinMax(minvalue, maxvalue);

    matrix_generator->SetDistanceValueMinMax(0, 10);
    matrix_generator->SetNumberOfBinsPerAxis(m_Bins);
    typename RunLengthFeatures::Pointer runLengthMatrixCalculator = RunLengthFeatures::New();
    typename RunLengthFeatures::Pointer runLengthFeaturesCalculator = RunLengthFeatures::New();
    typename OffsetVector::ConstIterator offsetIt;
    size_t offsetNum = 0;

    double sre = 0, lre = 0, gln = 0, glnn = 0, rln = 0, rlnn = 0, rp = 0, lglre = 0, hglre = 0, srlgle = 0, srhgle = 0, lrlgle = 0, lrhgle = 0;
    std::vector<double> features;

    for (offsetIt = offsets->Begin(); offsetIt != offsets->End(); offsetIt++, offsetNum++)
    {
        matrix_generator->SetOffset(offsetIt.Value());
        matrix_generator->Update();

        //auto outputglcmGenerator = matrix_generator->GetOutput();
        //auto iterator = outputglcmGenerator->Begin();
        //while (iterator != outputglcmGenerator->End())
        //{
        //  std::cout << iterator.GetIndex() << " " << iterator.GetFrequency() << std::endl;
        //  ++iterator;
        //}

        runLengthFeaturesCalculator->SetInput(matrix_generator->GetOutput());
        runLengthFeaturesCalculator->Update();

        sre += runLengthFeaturesCalculator->GetShortRunEmphasis();
        lre += runLengthFeaturesCalculator->GetLongRunEmphasis();
        gln += runLengthFeaturesCalculator->GetGreyLevelNonuniformity();
        rln += runLengthFeaturesCalculator->GetRunLengthNonuniformity();
        lglre += runLengthFeaturesCalculator->GetLowGreyLevelRunEmphasis();
        hglre += runLengthFeaturesCalculator->GetHighGreyLevelRunEmphasis();
        srlgle += runLengthFeaturesCalculator->GetShortRunLowGreyLevelEmphasis();
        srhgle += runLengthFeaturesCalculator->GetShortRunHighGreyLevelEmphasis();
        lrlgle += runLengthFeaturesCalculator->GetLongRunLowGreyLevelEmphasis();
        lrhgle += runLengthFeaturesCalculator->GetLongRunHighGreyLevelEmphasis();
    }
    sre /= offsets->size();
    lre /= offsets->size();
    gln /= offsets->size();
    rln /= offsets->size();
    lglre /= offsets->size();
    hglre /= offsets->size();
    srlgle /= offsets->size();
    srhgle /= offsets->size();
    lrlgle /= offsets->size();
    lrhgle /= offsets->size();

    features.push_back(sre);
    features.push_back(lre);
    features.push_back(gln);
    features.push_back(rln);
    features.push_back(lglre);
    features.push_back(hglre);
    features.push_back(srlgle);
    features.push_back(srhgle);
    features.push_back(lrlgle);
    features.push_back(lrhgle);

    return features;
}
template <class ImageType>
typename ImageType::Pointer PseudoProgressionEstimator::MakeAdditionalModality(typename ImageType::Pointer image1, typename ImageType::Pointer image2)
{
    typename ImageType::Pointer SubtractedImage = ImageType::New();
    SubtractedImage->CopyInformation(image1);
    SubtractedImage->SetRequestedRegion(image1->GetLargestPossibleRegion());
    SubtractedImage->SetBufferedRegion(image1->GetBufferedRegion());
    SubtractedImage->Allocate();
    SubtractedImage->FillBuffer(0);

    typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;
    IteratorType subIt(SubtractedImage, SubtractedImage->GetLargestPossibleRegion());
    IteratorType im1It(image1, image1->GetLargestPossibleRegion());
    IteratorType im2It(image2, image2->GetLargestPossibleRegion());

    im1It.GoToBegin();
    im2It.GoToBegin();
    subIt.GoToBegin();
    while (!im1It.IsAtEnd())
    {
        subIt.Set(im1It.Get() - im2It.Get());
        ++im1It;
        ++im2It;
        ++subIt;
    }
    return SubtractedImage;
}