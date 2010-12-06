#ifndef MeasureShearAlgo_H
#define MeasureShearAlgo_H

#include <vector>
#include "BVec.h"
#include "MyMatrix.h"
#include "Transformation.h"
#include "Image.h"
#include "Log.h"
#include "FittedPsf.h"
#include "Pixel.h"

void measureSingleShear(
    Position& cen, const Image<double>& im, double sky,
    const Transformation& trans, const FittedPsf& fitPsf,
    double noise, double gain, const Image<double>* weightIm, 
    double galAperture, double maxAperture,
    int galOrder, int galOrder2,
    double fPsf, double minGalSize, bool fixCen,
    double xOffset, double yOffset,
    bool fixSigma, double fixSigmaValue,
    OverallFitTimes* times, ShearLog& log,
    std::complex<double>& shear, 
    DSmallMatrix22& shearCov, BVec& shapelet,
    double& nu, long& flag);

void measureSingleShear1(
    Position& cen, const Image<double>& im, double sky,
    const Transformation& trans, const std::vector<BVec>& psf,
    double noise, double gain, const Image<double>* weightIm, 
    double galAperture, double maxAperture,
    int galOrder, int galOrder2,
    double fPsf, double minGalSize, bool fixCen,
    double xOffset, double yOffset,
    bool fixSigma, double fixSigmaValue,
    OverallFitTimes* times, ShearLog& log,
    std::complex<double>& shear, 
    DSmallMatrix22& shearCov, BVec& shapelet,
    double& nu, long& flag);

// This is poorly named.  It is still only measuring a single shear,
// but it is doing the multi-exposure version of the measurement.
// It is _not_ measuring multiple shears.
// The algorithm is basically the same as the above function, except that
// we already have the pixel values loaded up with a larger than needed 
// aperture.  Then we take subsets of that give the actual aperture size
// that we want each time.
void measureMultiShear(
    Position& cen, 
    const std::vector<PixelList>& allPix,
    const std::vector<BVec>& psf,
    double galAperture, double maxAperture,
    int galOrder, int galOrder2,
    double fPsf, double minGalSize, bool fixCen,
    OverallFitTimes* times, ShearLog& log,
    std::complex<double>& shear, 
    DSmallMatrix22& shearCov, BVec& shapelet,
    double& nu, long& flag);

#endif

