#ifndef CaPTkTraining_h
#define CaPTkTraining_h

// The following header file is generated by CMake and thus it's located in
// the build directory. It provides an export macro for classes and functions
// that you want to be part of the public interface of your module.
#include <MitkCaPTkTrainingExports.h>

#include <QObject>
#include <QFuture>
#include <QFutureWatcher>
#include <QProgressBar>

namespace captk {
/** \class captk::Training
 *  \brief Training Module API
 */
class MITKCAPTKTRAINING_EXPORT Training : public QObject
{
    Q_OBJECT

public:
    Training(QObject *parent = 0);

    ~Training();

    /** \brief Runs the algorithm
     * 
     * Execute the algorithm in a background thread. When the
     * algorithm finishes, OnAlgorithmFinished() is called.
     * 
     * @param featuresCsvPath path to the features csv file
     * @param responsesCsvPath path to the responses csv file
     * @param classificationKernelStr 
     * @param configurationStr ("Cross-validation", "Split Train/Test", "Train" or "Test")
     * @param folds number of folds (for configurationStr=="Cross-validation")
     * @param samples number of samples (for configurationStr=="Split Train/Test")
     * @param modelDirPath path to the model directory (for configurationStr=="Test") 
     * @param outputDirPath 
    */
    void Run(
        const QString featuresCsvPath,
        const QString responsesCsvPath,
        const QString classificationKernelStr,
        const QString configurationStr,
        const int     folds,
        const int     samples,
        const QString modelDirPath,
        const QString outputDirPath
    );

    void SetProgressBar(QProgressBar* progressBar);

    /** \struct Result
     *  \brief result of the execution of the algorithm
     * 
     * if ok == true, then everything went fine, 
     * else errorMessage is populated.
    */
    typedef struct Result 
    {
        bool ok = true;
        std::string errorMessage = "";
    } Result;

public slots:
    /** \brief This function runs in the main thread when 
     * the algorithm is finished
    */
    void OnAlgorithmFinished();

protected:

    /** \brief Runs the algorithm after the operations in Run
     * 
     * This can serve as a background thread. When the
     * algorithm finishes, OnAlgorithmFinished() is called.
     * The parameters are the same as Run()
     * 
     * @return the result struct (that contains the output or an errorMessage)
    */
    Result RunThread(
        const QString& featuresCsvPath,
        const QString& responsesCsvPath,
        const QString& classificationKernelStr,
        const QString& configurationStr,
        const int      folds,
        const int      samples,
        const QString& modelDirPath,
        const QString& outputDirPath
    );

    bool m_IsRunning = false;
    QFutureWatcher<Result> m_Watcher;
    QFuture<Result> m_FutureResult;

    QProgressBar* m_ProgressBar;
};

}

#endif // ! CaPTkTraining_h
