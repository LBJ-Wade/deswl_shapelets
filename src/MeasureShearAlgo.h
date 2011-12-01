#ifndef MeasureShearAlgo_H
#define MeasureShearAlgo_H

#include <vector>
#include "dbg.h"
#include "BVec.h"
#include "MyMatrix.h"
#include "Log.h"
#include "Pixel.h"

void DoMeasureShear(
    const std::vector<PixelList>& allpix,
    const std::vector<BVec>& psf,
    double galAperture, double maxAperture,
    int& galOrder, int galOrder2, int maxm, int minGalOrder, bool baseOrderOnNu,
    double minFPsf, double maxFPsf, double minGalSize, bool fixCen,
    bool fixSigma, double fixSigmaValue, bool nativeOnly,
    ShearLog& log, BVec& shapelet, 
    std::complex<double>& gamma, DSmallMatrix22& cov,
    double& nu, long& flag);

#endif

