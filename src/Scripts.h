
#include "Image.h"
#include "Transformation.h"
#include "InputCatalog.h"
#include "StarCatalog.h"
#include "FittedPsf.h"
#include "PsfCatalog.h"
#include "ShearCatalog.h"
#include "Log.h"

void DoFindStars(
    ConfigFile& params, FindStarsLog& log,
    const Image<double>& im, const Image<double>* weight_image,
    const Transformation& trans, 
    const InputCatalog& incat,
    std::auto_ptr<StarCatalog>& starcat);

void DoMeasurePsf(
    ConfigFile& params, PsfLog& log,
    const Image<double>& im, const Image<double>* weight_image,
    const Transformation& trans, const StarCatalog& starcat,
    std::auto_ptr<PsfCatalog>& psfcat, std::auto_ptr<FittedPsf>& fitpsf,
    double& sigmaP);

void DoMeasureShear(
    ConfigFile& params, ShearLog& log,
    const Image<double>& im, const Image<double>* weight_image,
    const Transformation& trans, 
    const InputCatalog& incat, const FittedPsf& fitpsf,
    std::auto_ptr<ShearCatalog>& shearcat);

void DoSplitStars(
    ConfigFile& params, std::string logFile, std::auto_ptr<Log>& log,
    const Image<double>& im, const Image<double>* weight_image,
    const Transformation& trans, 
    const InputCatalog& incat, const StarCatalog& starcat,
    double sigmaP);

