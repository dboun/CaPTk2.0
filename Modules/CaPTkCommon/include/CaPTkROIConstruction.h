#ifndef CaPTkROIConstruction_h
#define CaPTkROIConstruction_h

namespace captk
{
class ROIConstructionHelperBase;
}

#include "CaPTkROIConstructionHelper.h"

#include "mitkLabel.h"
#include "mitkLabelSetImage.h"

#include "itkImage.h"

#include <vector>
#include <string>
#include <memory>

// The following header file is generated by CMake and thus it's located in
// the build directory. It provides an export macro for classes and functions
// that you want to be part of the public interface of your module.
#include <MitkCaPTkCommonExports.h>

namespace captk
{
/** \class ROIConstruction 
 * \brief Takes an input label mask and lattice parameters and produces masks
 * 
 * These masks are each at the intersection of input mask and lattice point.
 * The masks can be obtained one at a time, using HasNext() and GetNext().
 * Update(...) should be called prior to attempting to get the masks.
 */
class MITKCAPTKCOMMON_EXPORT ROIConstruction
{
public:
    ROIConstruction();
    ~ROIConstruction();

    /** \brief Find all the different lattice points */
    void Update(
        mitk::LabelSetImage::Pointer input,
        captk::ROIConstructionHelper::MODE mode,
        float radius,
        float step
    );

    /** \brief Shows if there is another lattice mask left */
    bool IsAtEnd();

    /** \brief Get the next lattice mask 
     * \param rMask an empty, but initialized LabelSetImage to be populated
     * \return the weight of the ROI patch
    */
    float Get(mitk::LabelSetImage::Pointer& rMask);

    /** \brief Resets the index to the first lattice ROI */
    void GoToBegin();

    /** \brief ++ overloading to go to the next lattice */
    ROIConstruction &operator++() //suffix
    {
        this->m_CurrentIndex++; // actual operation
        return *this;
    }

    /** \brief ++ overloading to go to the next lattice */
    ROIConstruction operator++(int) //postfix(calls suffix)
    {
        ROIConstruction tmp(*this);
        operator++(); // call suffix
        return tmp;
    }

private:

    template < typename TPixel, unsigned int VImageDimension >
    void CreateHelper(
        typename itk::Image<TPixel,VImageDimension>* mask
    );

    std::shared_ptr<captk::ROIConstructionItkHelperBase> m_Helper;
    
    mitk::LabelSetImage::Pointer m_MaskTemplate;
    
    size_t m_CurrentIndex;
};
}

#include "CaPTkROIConstruction.hxx"

#endif // ! CaPTkROIConstruction_h