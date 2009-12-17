#ifndef CoaddCatalog_H
#define CoaddCatalog_H

#include <vector>
#include <string>
#include "Bounds.h"
#include "ConfigFile.h"

class CoaddCatalog 
{

public :

    CoaddCatalog(const ConfigFile& _params);
    ~CoaddCatalog();

    size_t size() const { return _id.size(); }

    void readCatalog();

    const std::vector<long>& getIdList() const { return _id; }
    const std::vector<Position>& getPosList() const { return _pos; }
    const std::vector<Position>& getSkyPosList() const { return _skyPos; }
    const std::vector<double>& getSkyList() const { return _sky; }
    const std::vector<long>& getFlagsList() const { return _flags; }
    const std::vector<float>& getMagList() const { return _mag; }
    const std::vector<float>& getMagErrList() const { return _magErr; }

    long getId(int i) const { return _id[i]; }
    Position getPos(int i) const { return _pos[i]; }
    Position getSkyPos(int i) const { return _skyPos[i]; }
    double getSky(int i) const { return _sky[i]; }
    long getFlags(int i) const { return _flags[i]; }
    double getMag(int i) const { return _mag[i]; }
    double getMagErr(int i) const { return _magErr[i]; }

    const Bounds& getSkyBounds() const { return _skyBounds; }

private :

    std::vector<long> _id;
    std::vector<Position> _pos;
    std::vector<Position> _skyPos;

    std::vector<double> _sky;

    std::vector<long> _flags;

    std::vector<float> _mag;
    std::vector<float> _magErr;

    Bounds _skyBounds;

    const ConfigFile& _params;

};

#endif
