#ifndef _SXCAT_H
#define _SXCAT_H

#include <vector>
#include "BVec.h"
//#include "FittedPSF.h"
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
  std::vector<long> flags;
  std::vector<float> ra;
  std::vector<float> dec;

  // not read from file
  //std::vector<double> sigma;
  //std::vector<long> size_flags;
  //std::vector<long> star_flag;

  // These unused
  std::vector<double> noise;
} SXCAT_STRUCT;



void WriteMainKeywords(FitsFile& fits, const ConfigFile& params);

void ResizeSXCat(SXCAT_STRUCT& cat, long n);
void ReadSXCat(const ConfigFile& params, SXCAT_STRUCT& cat);


// FindStars 
typedef struct {
  std::vector<long> id;
  std::vector<double> sigma0;
  std::vector<long> size_flags;
  std::vector<long> star_flag;

  // These not output/input, just for convenience
  std::vector<double> local_sky;
  std::vector<Position> pos;
  std::vector<double> noise;
} FINDSTARS_STRUCT;

void ReadFindStarsCat(const ConfigFile& params, FINDSTARS_STRUCT& cat);
void WriteFindStarsKeywords(FitsFile& fits, const ConfigFile& params);
void WriteFindStarsCat(const ConfigFile& params, FINDSTARS_STRUCT& cat);
void ResizeFindStarsCat(FINDSTARS_STRUCT& cat, size_t n);


//MeasurePSF

typedef struct {

  std::vector<long> id;
  std::vector<long> psf_flags;
  std::vector<double> nu;
  std::vector<long> psf_order;
  std::vector<double> sigma_p;
  std::vector<BVec> psf;

} PSF_STRUCT;

void ResizePSFCat(PSF_STRUCT& cat, size_t n, int psf_order, double sigma=0.0);
void WritePSFKeywords(FitsFile& fits, const ConfigFile& params);
void WritePSFCat(const ConfigFile& params, PSF_STRUCT& cat);
void ReadPSFCat(const ConfigFile& params, PSF_STRUCT& cat);


#if 0
void WriteFittedPSF(const ConfigFile& params, FittedPSF& fpsf);
void ReadFittedPSF(const ConfigFile& params, FittedPSF& fpsf);
#endif


typedef struct {

  std::vector<long> id;
  std::vector<long> shear_flags;
  std::vector<double> shear1;
  std::vector<double> shear2;
  std::vector<double> shear_cov00;
  std::vector<double> shear_cov01;
  std::vector<double> shear_cov11;

  std::vector<long> gal_order;
  std::vector<BVec> shapelets_prepsf;

  // These are extra in order to get around Joe's stupidity
#ifdef SHEXTRA_PARS
  std::vector<double> sigma0;
  std::vector<long> size_flags;
  std::vector<long> star_flag;
#endif
} SHEAR_STRUCT;

void ResizeShearCat(SHEAR_STRUCT& cat, size_t n, int psf_order, 
    double sigma=0.0);
void WriteShearKeywords(FitsFile& fits, const ConfigFile& params);
void WriteShearCat(const ConfigFile& params, SHEAR_STRUCT& cat);
void ReadShearCat(const ConfigFile& params, SHEAR_STRUCT& cat);

#endif
