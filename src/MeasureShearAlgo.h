#ifndef MeasureShearAlgo_H
#define MeasureShearAlgo_H

#include <vector>
#include "dbg.h"
#include "BVec.h"
#include "MyMatrix.h"
#include "Log.h"
#include "Pixel.h"

void measureSingleShear(
    const std::vector<PixelList>& allpix,
    const std::vector<BVec>& psf,
    int& galOrder, const ConfigFile& params,
    ShearLog& log, BVec& shapelet, 
    std::complex<double>& gamma, DSmallMatrix22& cov,
    double& nu, long& flag);

#endif

