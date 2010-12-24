
#include "Image.h"
#include "Transformation.h"
#include "InputCatalog.h"
#include "StarCatalog.h"
#include "FittedPsf.h"
#include "PsfCatalog.h"
#include "ShearCatalog.h"
#include "Log.h"

void doFindStars(
    ConfigFile& params, FindStarsLog& log,
    const Image<double>& im, const Image<double>* weightIm,
    const Transformation& trans, 
    const InputCatalog& inCat,
    std::auto_ptr<StarCatalog>& starCat);

void doMeasurePsf(
    ConfigFile& params, PsfLog& log,
    const Image<double>& im, const Image<double>* weightIm,
    const Transformation& trans, const StarCatalog& starCat,
    std::auto_ptr<PsfCatalog>& psfCat, std::auto_ptr<FittedPsf>& fitPsf,
    double& sigmaP);

void doMeasureShear(
    ConfigFile& params, ShearLog& log,
    const Image<double>& im, const Image<double>* weightIm,
    const Transformation& trans, 
    const InputCatalog& inCat, const FittedPsf& fitPsf,
    std::auto_ptr<ShearCatalog>& shearCat);

void doSplitStars(
    ConfigFile& params, std::string logFile, std::auto_ptr<Log>& log,
    const Image<double>& im, const Image<double>* weightIm,
    const Transformation& trans, 
    const InputCatalog& inCat, const StarCatalog& starCat,
    double sigmaP);

