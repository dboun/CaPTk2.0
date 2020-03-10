#ifndef CaPTkItkAdaptiveGeodesicDistanceFilter_h
#define CaPTkItkAdaptiveGeodesicDistanceFilter_h

// The following header file is generated by CMake and thus it's located in
// the build directory. It provides an export macro for classes and functions
// that you want to be part of the public interface of your module.
#include <MitkCaPTkCommonExports.h>

#include "itkImage.h"
#include "itkImageToImageFilter.h"

namespace itk
{
namespace captk
{
/** \class AdaptiveGeodesicDistanceFilter
 * \brief Applies the Adaptive Geodesic Distance filter to image 
 * 
 * This map reflects a composite of intensity and spatial distance 
 * from the user-drawn regions: 
 * voxels which are far away and/or have very different image intensities 
 * are less likely to be part of the ROI.
 * UseInputImageAsMask is ON by default.
 * In the scenario that mask is set and UseInputImageAsMask is ON, 
 * the set mask takes precedence. If mask is not set and UseInputImageAsMask
 * is off, calculations happen for the largest possible image region.
 * Reference:
 * @inproceedings{gaonkar2014adaptive,
 * title={Adaptive geodesic transform for segmentation of vertebrae on CT images},
 * author={Gaonkar, Bilwaj and Shu, Liao and Hermosillo, Gerardo and Zhan, Yiqiang},
 * booktitle={SPIE Medical Imaging},
 * pages={903516--903516},
 * year={2014},
 * organization={International Society for Optics and Photonics}
 * }
 */
template <typename TInputImage, typename TOutputImage>
class MITKCAPTKCOMMON_EXPORT AdaptiveGeodesicDistanceFilter 
          : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /*---- Code needed for itk filters ----*/

  ITK_DISALLOW_COPY_AND_ASSIGN(AdaptiveGeodesicDistanceFilter);
 
  using Self = AdaptiveGeodesicDistanceFilter;
  using Superclass = ImageToImageFilter<TInputImage, TInputImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  itkNewMacro(Self);
 
  itkTypeMacro(AdaptiveGeodesicDistanceFilter, ImageToImageFilter);
 
  using ImageType = TInputImage;
  using PixelType = typename ImageType::PixelType;

  /*---- Custom types ----*/

  using LabelsImageType = itk::Image<int, TInputImage::ImageDimension>;

  /*---- Setters and Getters ----*/

  itkSetMacro(LabelOfInterest, PixelType);
  itkGetConstMacro(LabelOfInterest, PixelType);

  itkSetMacro(LimitAt255, bool);
  itkGetConstMacro(LimitAt255, bool);
  itkBooleanMacro(LimitAt255); // creates name##On() and name##Off()

  itkSetMacro(UseInputImageAsMask, bool);
  itkGetConstMacro(UseInputImageAsMask, bool);
  itkBooleanMacro(UseInputImageAsMask);

  itkSetMacro(Mask, typename TInputImage::Pointer)

  itkSetMacro(Labels, typename LabelsImageType::Pointer);
 
protected:
  AdaptiveGeodesicDistanceFilter(); // protected to force factory usage

  /** \brief Used internally to separate the actual calculations */
  static
  typename TOutputImage::Pointer
  InternalProcessing(
    const TInputImage* inputImage,
    const TInputImage* mask,
    const LabelsImageType* labels,
    const int labelOfInterest,
    bool limitAt255);

protected:

  void 
  GenerateData() override;

private:

  /// Calculation will ignore other labels present in the "labels" image, apart from this label
  PixelType m_LabelOfInterest = 0;

  /// Limit output values to [0,255] range. Makes it a more useful feature. Unlimited otherwise.
  bool m_LimitAt255 = false;

  /// Calculation will happen only on the non-zero pixels/voxels of the image (on by default)
  bool m_UseInputImageAsMask = true;

  /// Computations will happen only for pixels/voxels that are not zero in this mask
  typename TInputImage::Pointer m_Mask;

  /// Holds the labels indicating the different regions (can be scribbles or a segmentation)
  typename LabelsImageType::Pointer m_Labels;
};

}
}

#ifndef ITK_MANUAL_INSTANTIATION
#  include "CaPTkItkAdaptiveGeodesicDistanceFilter.hxx"
#endif

#endif // ! CaPTkItkAdaptiveGeodesicDistanceFilter_h