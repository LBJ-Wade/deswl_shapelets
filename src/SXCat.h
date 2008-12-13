#ifndef _SXCAT_H
#define _SXCAT_H

#include <vector>
#include "BVec.h"
#include "FittedPSF.h"
#include "FitsFile.h"

#include "dbg.h"


// for Position
#include "Bounds.h"

#include "ConfigFile.h"
#include "Name.h"

// for error codes
#include "Params.h"

#include <sstream>

// Dealing with the SExtractor catalog
// This is all tags we think we might need
typedef struct {
  std::vector<long> id;

  std::vector<Position> pos;
  std::vector<float> x;
  std::vector<float> y;
  std::vector<double> local_sky;
  std::vector<float> mag;
  std::vector<float> mag_err;
  std::vector<int> flags;

  // not read from file
  //std::vector<double> sigma;
  //std::vector<int> size_flags;
  //std::vector<int> star_flag;

  // These unused
  std::vector<double> noise;
} SXCAT_STRUCT;

void ResizeSXCat(SXCAT_STRUCT& cat, long n);
void ReadSXCat(ConfigFile& params, SXCAT_STRUCT& cat);


// FindStars 
typedef struct {
  std::vector<long> id;
  std::vector<double> sigma0;
  std::vector<int> size_flags;
  std::vector<int> star_flag;

  // These not output/input, just for convenience
  std::vector<double> local_sky;
  std::vector<Position> pos;
  std::vector<double> noise;
} FINDSTARS_STRUCT;

void ReadFindStarsCat(ConfigFile& params, FINDSTARS_STRUCT& cat);
void WriteFindStarsCat(ConfigFile& params, FINDSTARS_STRUCT& cat);
void ResizeFindStarsCat(FINDSTARS_STRUCT& cat, size_t n);


//MeasurePSF

typedef struct {

  std::vector<long> id;
  std::vector<int> psf_flags;
  std::vector<double> nu;
  std::vector<int> psf_order;
  std::vector<double> sigma_p;
  std::vector<BVec> psf;

} PSF_STRUCT;

void ResizePSFCat(PSF_STRUCT& cat, size_t n, int psf_order, double sigma=0.0);
void WritePSFCat(ConfigFile& params, PSF_STRUCT& cat);
void ReadPSFCat(ConfigFile& params, PSF_STRUCT& cat);


void WriteFittedPSF(ConfigFile& params, FittedPSF& fpsf);
void ReadFittedPSF(ConfigFile& params, FittedPSF& fpsf);


typedef struct {

  std::vector<long> id;
  std::vector<int> shear_flags;
  std::vector<double> shear1;
  std::vector<double> shear2;
  std::vector<double> shear_cov00;
  std::vector<double> shear_cov01;
  std::vector<double> shear_cov11;

  std::vector<int> gal_order;
  std::vector<BVec> shapelets_prepsf;

  // These are extra in order to get around Joe's stupidity
#ifdef SHEXTRA_PARS
  std::vector<double> sigma0;
  std::vector<int> size_flags;
  std::vector<int> star_flag;
#endif
} SHEAR_STRUCT;

void ResizeShearCat(SHEAR_STRUCT& cat, size_t n, int psf_order, 
    double sigma=0.0);
void WriteShearCat(ConfigFile& params, SHEAR_STRUCT& cat);
void ReadShearCat(ConfigFile& params, SHEAR_STRUCT& cat);

#endif
