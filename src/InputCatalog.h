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
    // The contructor just sets up the basic parameter information. 
    // It is usually followed by read() or something similar.
    InputCatalog(ConfigFile& params, const Image<double>* im=0);


    // ESS Factor the code a bit.  This allows empty constructor and delayed
    // load of params and data
    InputCatalog() {};
    //void init(ConfigFile& params);
    void init(ConfigFile& params, const Image<double>* im=0);
    void loadParams(std::string file) {
      _params.load(file);
    }
    void loadParams(const ConfigFile& params) {
      // ESS I don't know how else to do this
      std::stringstream s;
      s<<params;
      _params.read(s);
    }
    void determineNoiseMethod();


    size_t size() const { return _pos.size(); }

    void read();
    void readFits(std::string file);
    void readFits(std::string file, int hdu);
    void readAscii(std::string file, std::string delim = "  ");

    const std::vector<long>& getIdList() const { return _id; }
    const std::vector<Position>& getPosList() const { return _pos; }
    const std::vector<double>& getSkyList() const { return _sky; }
    const std::vector<float>& getMagList() const { return _mag; }
    const std::vector<float>& getMagErrList() const { return _magErr; }
    const std::vector<double>& getObjSizeList() const { return _objSize; }
    const std::vector<long>& getFlagsList() const { return _flags; }
    const std::vector<Position>& getSkyPosList() const { return _skyPos; }
    const std::vector<double>& getNoiseList() const { return _noise; }

    long getId(int i) const { return _id[i]; }
    Position getPos(int i) const { return _pos[i]; }
    double getSky(int i) const { return _sky[i]; }
    double getMag(int i) const { return _mag[i]; }
    double getMagErr(int i) const { return _magErr[i]; }
    double getObjSize(int i) const { return _objSize[i]; }
    long getFlags(int i) const { return _flags[i]; }
    Position getSkyPos(int i) const { return _skyPos[i]; }
    double getNoise(int i) const { return _noise[i]; }

    const Bounds& getBounds() const { return _bounds; }
    const Bounds& getSkyBounds() const { return _skyBounds; }

    void printall(int i);

private :

    //ConfigFile& _params;
    ConfigFile _params;

    std::vector<long> _id;
    std::vector<Position> _pos;
    std::vector<double> _sky;
    std::vector<float> _mag;
    std::vector<float> _magErr;
    std::vector<double> _objSize;
    std::vector<long> _flags;
    std::vector<Position> _skyPos;
    std::vector<double> _noise;

    Bounds _bounds;
    Bounds _skyBounds;

    const Image<double>* _im;

    enum NoiseMethod {
        VALUE, CATALOG, CATALOG_SIGMA, GAIN_VALUE, GAIN_FITS, WEIGHT_IMAGE
    };
    NoiseMethod _nm;
    double _noiseValue;
    double _gain;
    double _readNoise;
};

#endif
