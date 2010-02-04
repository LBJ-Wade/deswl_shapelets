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
    Position cen, const Image<double>& im, double sky,
    const Transformation& trans, const FittedPsf& fitPsf,
    double noise, double gain, const Image<double>* weightIm, 
    double galAperture, double maxAperture,
    int galOrder, int galOrder2,
    double fPsf, double minGalSize,
    OverallFitTimes* times, ShearLog& log,
    std::complex<double>& shear, 
    DSmallMatrix22& shearCov, BVec& shapelet,
    double& nu, long& flag);

void measureSingleShear1(
    Position cen, const Image<double>& im, double sky,
    const Transformation& trans, const std::vector<BVec>& psf,
    double noise, double gain, const Image<double>* weightIm, 
    double galAperture, double maxAperture,
    int galOrder, int galOrder2,
    double fPsf, double minGalSize,
    OverallFitTimes* times, ShearLog& log,
    std::complex<double>& shear, 
    DSmallMatrix22& shearCov, BVec& shapelet,
    double& nu, long& flag);

void measureMultiShear(
    const Position& cen, 
    const std::vector<PixelList>& allPix,
    const std::vector<BVec>& psf,
    double galAperture, double maxAperture,
    int galOrder, int galOrder2,
    double fPsf, double minGalSize,
    OverallFitTimes* times, ShearLog& log,
    std::complex<double>& shear, 
    DSmallMatrix22& shearCov, BVec& shapelet,
    double& nu, long& flag);

#endif

