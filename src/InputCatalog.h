#ifndef InputCatalog_H
#define InputCatalog_H

#include <vector>
#include <string>
#include "Bounds.h"
#include "ConfigFile.h"
#include "Image.h"

class InputCatalog 
{

  public :

    // If you already have the image loaded, you can pass it here.
    // It is only needed for global sky calculation from the image median, 
    // so if you haven't loaded it, and you don't have any bad local
    // sky values, it won't need the global_sky value, so it won't
    // load the image here either.
    InputCatalog(ConfigFile& _params, const Image<double>* im=0);

    size_t size() const { return pos.size(); }

    void Read();
    void ReadFits(std::string file);
    void ReadAscii(std::string file, std::string delim = "  ");

    // Leave these public, rather than use Get and Set methods.
    std::vector<long> id;
    std::vector<Position> pos;
    std::vector<double> sky;
    std::vector<float> mag;
    std::vector<float> mag_err;
    std::vector<double> objsize;
    std::vector<long> flags;
    std::vector<float> ra;
    std::vector<float> dec;
    std::vector<double> noise;

  private :

    const ConfigFile& params;

};

#endif
