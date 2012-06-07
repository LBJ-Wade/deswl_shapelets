#ifndef CoaddCatalog_H
#define CoaddCatalog_H

#include <vector>
#include <string>
#include "dbg.h"
#include "Bounds.h"
#include "ConfigFile.h"

class CoaddCatalog 
{

public :

    // Just sets up the parameters.  Normally followed by read().
    CoaddCatalog(const ConfigFile& _params);
    ~CoaddCatalog();

    int size() const { return _id.size(); }

    void read();
    void readFits(const std::string& file);

    const std::vector<long>& getIdList() const { return _id; }
    const std::vector<Position>& getPosList() const { return _pos; }
    const std::vector<Position>& getSkyPosList() const { return _skypos; }
    const std::vector<double>& getSkyList() const { return _sky; }
    const std::vector<long>& getFlagsList() const { return _flags; }
    const std::vector<double>& getMagList() const { return _mag; }
    const std::vector<double>& getMagErrList() const { return _mag_err; }

    long getId(int i) const { return _id[i]; }
    Position getPos(int i) const { return _pos[i]; }
    Position getSkyPos(int i) const { return _skypos[i]; }
    double getSky(int i) const { return _sky[i]; }
    long getFlags(int i) const { return _flags[i]; }
    double getMag(int i) const { return _mag[i]; }
    double getMagErr(int i) const { return _mag_err[i]; }

    const Bounds& getSkyBounds() const { return _skybounds; }

private :

    std::vector<long> _id;
    std::vector<Position> _pos;
    std::vector<Position> _skypos;

    std::vector<double> _sky;

    std::vector<long> _flags;

    std::vector<double> _mag;
    std::vector<double> _mag_err;

    Bounds _skybounds;

    const ConfigFile& _params;

};

#endif
