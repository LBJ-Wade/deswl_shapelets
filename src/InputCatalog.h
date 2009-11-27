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
    InputCatalog(ConfigFile& params, const Image<double>* im=0);

    size_t size() const { return _pos.size(); }

    void read();
    void readFits(std::string file);
    void readAscii(std::string file, std::string delim = "  ");

    const std::vector<long>& getIdList() const { return _id; }
    const std::vector<Position>& getPosList() const { return _pos; }
    const std::vector<double>& getSkyList() const { return _sky; }
    const std::vector<float>& getMagList() const { return _mag; }
    const std::vector<float>& getMagErrList() const { return _magErr; }
    const std::vector<double>& getObjSizeList() const { return _objSize; }
    const std::vector<long>& getFlagsList() const { return _flags; }
    const std::vector<float>& getRaList() const { return _ra; }
    const std::vector<float>& getDeclList() const { return _decl; }
    const std::vector<double>& getNoiseList() const { return _noise; }

    long getId(int i) const { return _id[i]; }
    Position getPos(int i) const { return _pos[i]; }
    double getSky(int i) const { return _sky[i]; }
    double getMag(int i) const { return _mag[i]; }
    double getMagErr(int i) const { return _magErr[i]; }
    double getObjSizeList(int i) const { return _objSize[i]; }
    long getFlags(int i) const { return _flags[i]; }
    Position getSkyPos(int i) const { return Position(_ra[i],_decl[i]); }
    double getNoise(int i) const { return _noise[i]; }

    const Bounds& getBounds() const { return _bounds; }

private :

    const ConfigFile& _params;

    std::vector<long> _id;
    std::vector<Position> _pos;
    std::vector<double> _sky;
    std::vector<float> _mag;
    std::vector<float> _magErr;
    std::vector<double> _objSize;
    std::vector<long> _flags;
    std::vector<float> _ra;
    std::vector<float> _decl;
    std::vector<double> _noise;

    Bounds _bounds;

};

#endif
