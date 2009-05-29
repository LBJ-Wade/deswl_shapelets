#ifndef PIXEL_H
#define PIXEL_H

#include <complex>
#include <vector>
#include <string>
#include "Image.h"
#include "Transformation.h"
#include "ConfigFile.h"

struct Pixel { 
  Pixel(double _u, double _v, double _I, double _wt) :
    z(_u,_v), I(_I), wt(_wt) {}
  std::complex<double> z;
  double I,wt;
};

void GetPixList(const Image<double>& im, std::vector<Pixel>& pix,
    const Position cen, double sky, double noise, double gain,
    const Image<double>* wt_im, const Transformation& trans,
    double aperture, long& flag);

void GetSubPixList(std::vector<Pixel>& pix,
    const std::vector<Pixel>& allpix,
    double aperture, long& flag);

#endif
