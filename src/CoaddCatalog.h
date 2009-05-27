#ifndef CoaddCatalog_H
#define CoaddCatalog_H

#include <vector>
#include <string>
#include "Bounds.h"
#include "ConfigFile.h"
#include "Image.h"
#include "Transformation.h"
#include "dbg.h"
#include <CCfits/CCfits>
#include "TMV.h"

class CoaddCatalog 
{

  public :

    CoaddCatalog(ConfigFile& _params);

    size_t size() const { return pos.size(); }

    void ReadCatalog();
    void ReadPixelLists();

    // Leave these public, rather than use Get and Set methods.
    std::vector<long> id;
    std::vector<Position> pos;
    std::vector<Position> skypos;

    std::vector<double> sky;

    std::vector<long> flags;

    std::vector<float> ra;
    std::vector<float> dec;

    std::vector<float> mag;
    std::vector<float> mag_err;

    //std::vector<double> noise;

  private :

    ConfigFile params;

};

#endif
