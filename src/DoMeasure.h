#ifndef DoMeasure_H
#define DoMeasure_H

#include "ConfigFile.h"
#include "Image.h"
#include "TimeVars.h"
#include "BVec.h"
#include "TMV_Small.h"
#include "Transformation.h"

// Returns how many successful measurements
int DoMeasurePSF(ConfigFile& params);
int DoMeasureShear(ConfigFile& params);

void ReadCatalog(ConfigFile& params,
    std::vector<Position>& all_pos, std::vector<double>& all_sky,
    std::vector<double>& all_noise, double& gain, Image<double>*& weight_im);

void MeasureSingleShear(
    Position cen, const Image<double>& im, double sky, 
    const Transformation& trans,
    const std::vector<BVec>& psf,
    double noise, double gain, const Image<double>* weight_im,
    double gal_aperture, double max_aperture,
    int gal_order, int gal_order2,
    double f_psf, double min_gal_size,
    OverallFitTimes* times,
    bool& success_shear, std::complex<double>& shear, 
    tmv::Matrix<double>& varshear, BVec*& shapelet);

double EstimateSigma(
    const Image<double>& im,
    const std::vector<Position>& all_pos, const std::vector<double>& all_sky,
    const std::vector<double>& all_noise, double gain,
    const Image<double>* weight_im, const Transformation& trans, double psfap);

void MeasureSinglePSF(
    Position cen, const Image<double>& im, double sky,
    const Transformation& trans,
    double noise, double gain, const Image<double>* weight_im,
    double sigma_p, double psf_aperture, int psf_order,
    BVec*& psf, double& nu);

#endif