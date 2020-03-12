#ifndef CaPTkInteractiveSegmentation_h
#define CaPTkInteractiveSegmentation_h

// The following header file is generated by CMake and thus it's located in
// the build directory. It provides an export macro for classes and functions
// that you want to be part of the public interface of your module.
#include <MitkCaPTkInteractiveSegmentationExports.h>

#include "mitkImage.h"
#include "mitkLabelSetImage.h"

#include <QObject>
#include <QFuture>
#include <QFutureWatcher>

/** \class CaPTkInteractiveSegmentation
 *  \brief Singleton class that runs the interactive segmentation 
 * algorithm and adds the result to the data storage
 */
class MITKCAPTKINTERACTIVESEGMENTATION_EXPORT CaPTkInteractiveSegmentation : 
                                                    public QObject
{
    Q_OBJECT

public:
    CaPTkInteractiveSegmentation(QObject *parent = 0);

    ~CaPTkInteractiveSegmentation() {}

    /** \brief Runs the algorithm
     * 
     * Execute the algorithm in a background thread. When the
     * algorithm finishes, OnAlgorithmFinished() is called.
     * 
     * @param images a list of the co-registered input images
     * @param labels label image that contains the user drawn seeds
    */
    void Run(std::vector<mitk::Image::Pointer>& images, 
             mitk::LabelSetImage::Pointer& seeds);

signals:
    /** \brief This gets emitted when the async execution finishes */
    void finished(bool ok, std::string errorMessage, mitk::LabelSetImage::Pointer result);

    void progressUpdate(int progress);

protected slots:
    /** \brief This function runs in the main thread when 
     * the algorithm is finished to add the result to the data storage
    */
    void OnAlgorithmFinished();

    void OnProgressUpdateInternalReceived(int progress);

protected:

    /** \struct Result
     *  \brief result of the execution of the algorithm
     * 
     * if ok == true, then segmentation is populated, 
     * else errorMessage is populated.
    */
    typedef struct Result 
    {
        mitk::LabelSetImage::Pointer segmentation;
        bool ok = true;
        std::string errorMessage = "";
    } Result;

    /** \brief Runs the algorithm after the operations in Run
     * 
     * This can serve as a background thread. When the
     * algorithm finishes, OnAlgorithmFinished() is called.
     * 
     * @param images a list of the co-registered input images
     * @param seeds label image that contains the user drawn seeds
     * @return the result struct (that contains the output or an errorMessage)
    */
    Result RunThread(std::vector<mitk::Image::Pointer>& images, 
                     mitk::LabelSetImage::Pointer& seeds);

    bool m_IsRunning = false;
    QFutureWatcher<Result> m_Watcher;
    QFuture<Result> m_FutureResult;
};

#endif // ! CaPTkInteractiveSegmentation_h
