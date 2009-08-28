#ifndef CoaddCatalog_H
#define CoaddCatalog_H

#include <vector>
#include <string>
#include "TMV.h"
#include "Bounds.h"
#include "ConfigFile.h"
#include "Image.h"
#include "Pixel.h"
#include "Transformation.h"
#include "FittedPSF.h"
#include "dbg.h"
#include "Name.h"
#include "WlVersion.h"

class CoaddCatalog 
{

  public :

    CoaddCatalog(ConfigFile& _params);
    ~CoaddCatalog();

    size_t size() const { return id.size(); }

    void ReadCatalog();

    // Leave these public, rather than use Get and Set methods.
    std::vector<long> id;
    std::vector<Position> pos;
    std::vector<Position> skypos;

    std::vector<double> sky;
    std::vector<double> noise;

    std::vector<long> flags;

    std::vector<float> ra;
    std::vector<float> dec;

    std::vector<float> mag;
    std::vector<float> mag_err;

    Bounds skybounds;

  private :

    const ConfigFile& params;

};

#endif
