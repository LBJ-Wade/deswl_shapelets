#ifndef StarCatalog_H
#define StarCatalog_H

#include <vector>
#include <string>
#include "Bounds.h"
#include "ConfigFile.h"
#include "InputCatalog.h"
#include "Image.h"
#include "Transformation.h"
#include "Log.h"

// This function is also used by PSFCatalog.
void calculateSigma(
    double& sigma, // Initial value -- use <=0 if no initial guess
    const Image<double>& im, const Position& pos, double sky, 
    double noise, double gain, const Image<double>* weightIm,
    const Transformation& trans, double psfAp, double xOffset, double yOffset,
    long& flag, bool shouldUseShapeletSigma);

class StarCatalog
{
public:
    // Make from incat
    // fs_prefix is the prefix of the keywords for the 
    // parameters used by the findstars algorithm.
    StarCatalog(
        const InputCatalog& inCat,
        const ConfigFile& params, std::string fsPrefix = "stars_");

    // copy the input StarCatalog
    StarCatalog(
        const StarCatalog& inStarCat);

    // In this version, we take a subset of the input star catalog.
    // note indices must be sorted!
    StarCatalog(
        const StarCatalog& inStarCat, const std::vector<long> indices);

    // Setup the parameters.  Normally followed by read() or similar.
    StarCatalog(const ConfigFile& params, std::string fsPrefix = "stars_");

    // a deterministic seed is generated based on size of catalog and number
    // of stars
    void splitInTwo(const std::string f1, const std::string f2) const;
    // seed explicitly sent
    void splitInTwo(
        const std::string f1, const std::string f2,
        const unsigned int seed) const;

    size_t size() const { return _id.size(); }
    void read();
    void write() const;

    void readFits(std::string file);
    void readFitsOld(std::string file);
    void readAscii(std::string file, std::string delim = "  ");
    void writeFitsOld(std::string file) const;
    void writeFits(std::string file) const;
    void writeAscii(std::string file, std::string delim = "  ") const;

    void calculateSizes(
        const Image<double>& im, 
        const Image<double>* weight_im, const Transformation& trans);

    int findStars(FindStarsLog& log);

    const std::vector<long>& getIdList() const { return _id; }
    const std::vector<Position>& getPosList() const { return _pos; }
    const std::vector<double>& getSkyList() const { return _sky; }
    const std::vector<double>& getNoiseList() const { return _noise; }
    const std::vector<long>& getFlagsList() const { return _flags; }
    const std::vector<float>& getMagList() const { return _mag; }
    const std::vector<double>& getObjSizeList() const { return _objSize; }
    const std::vector<bool>& getIsAStarList() const { return _isAStar; }


    const ConfigFile& getParams() const { return _params; }
    template <typename T>
    void setPar(const std::string key, const T& val) {
     _params.add(key,val); 
    }

    long     getId(int i)      const { return _id[i]; }
    Position getPos(int i)     const { return _pos[i]; }
    double   getSky(int i)     const { return _sky[i]; }
    double   getNoise(int i)   const { return _noise[i]; }
    long     getFlags(int i)   const { return _flags[i]; }
    float    getMag(int i)     const { return _mag[i]; }
    double   getObjSize(int i) const { return _objSize[i]; }
    bool     isAStar(int i)    const { return _isAStar[i]; }
    bool     getIsAStar(int i) const { return _isAStar[i]; }

    void printall(int i);

private :

    std::vector<long> _id;
    std::vector<Position> _pos;
    std::vector<double> _sky;
    std::vector<double> _noise;
    std::vector<long> _flags;

    std::vector<float> _mag;
    std::vector<double> _objSize;
    std::vector<bool> _isAStar;

    //const ConfigFile& _params;
    // We need to be able to alter this if we want to write alternative
    // names.  So we'll make our own copy
    ConfigFile _params;
    std::string _prefix;
};

#endif
