#ifndef _SXCAT_H
#define _SXCAT_H

#include <vector>
#include "fitsio.h"
#include "FitsFile.h"

#include "dbg.h"
// for Position
#include "Bounds.h"

#include "ConfigFile.h"
#include "Name.h"

// for error codes
#include "Params.h"

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
  std::vector<double> sigma;
  std::vector<int> size_flags;
  std::vector<int> star_flag;

  // These unused
  std::vector<double> noise;
} SXCAT_STRUCT;

void ResizeSXCat(SXCAT_STRUCT& cat, long n);
void ReadSXCat(ConfigFile& params, SXCAT_STRUCT& cat);

#endif
