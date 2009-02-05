#ifndef InputCatalog_H
#define InputCatalog_H

#include <vector>
#include <string>
#include "Bounds.h"
#include "ConfigFile.h"

class InputCatalog 
{

  public :

    InputCatalog(ConfigFile& _params, std::string key_prefix);

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
    std::string prefix;

};

#endif
