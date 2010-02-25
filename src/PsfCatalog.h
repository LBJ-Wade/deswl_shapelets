#ifndef PsfCatalog_H
#define PsfCatalog_H

#include <sstream>
#include <vector>
#include "BVec.h"
#include "Bounds.h"
#include "ConfigFile.h"
#include "StarCatalog.h"
#include "Log.h"
#include "Image.h"
#include "Transformation.h"

class PsfCatalog 
{
public :

    // Make from starcat
    PsfCatalog(const StarCatalog& starCat, const ConfigFile& params);

    // Read from file
    PsfCatalog(const ConfigFile& params);
    // this one we don't have to deal with root= stuff
    PsfCatalog(const ConfigFile& params, std::string file);

    size_t size() const { return _id.size(); }

    void read();
    void read(std::string file);
    void readFits(std::string file);
    void readAscii(std::string file, std::string delim = "  ");

    void write() const;
    void writeFits(std::string file) const;
    void writeAscii(std::string file, std::string delim = "  ") const;

    double estimateSigma(
        const Image<double>& im,
        const Image<double>* weightIm, const Transformation& trans);

    int measurePsf(
        const Image<double>& im, const Image<double>* weightIm,
        const Transformation& trans, double sigma_p, PsfLog& log);

    const std::vector<long>& getIdList() const { return _id; }
    const std::vector<Position>& getPosList() const { return _pos; }
    const std::vector<double>& getSkyList() const { return _sky; }
    const std::vector<double>& getNoiseList() const { return _noise; }
    const std::vector<long>& getFlagsList() const { return _flags; }
    const std::vector<double>& getNuList() const { return _nu; }
    const std::vector<BVec>& getPsfList() const { return _psf; }

    long getId(int i) const { return _id[i]; }
    Position getPos(int i) const { return _pos[i]; }
    double getSky(int i) const { return _sky[i]; }
    double getNoise(int i) const { return _noise[i]; }
    long getFlags(int i) const { return _flags[i]; }
    double getNu(int i) const { return _nu[i]; }
    const BVec& getPsf(int i) const { return _psf[i]; }

    void setFlag(int i, long newFlag) { _flags[i] |= newFlag; }
    std::vector<long>& getFlagsList() { return _flags; }

private :

    const ConfigFile& _params;

    std::vector<long> _id;
    std::vector<Position> _pos;
    std::vector<double> _sky;
    std::vector<double> _noise;
    std::vector<long> _flags;
    std::vector<double> _nu;
    EIGEN_mutable std::vector<BVec> _psf;

};

void measureSinglePsf(
    Position cen, const Image<double>& im, double sky,
    const Transformation& trans,
    double noise, double gain, const Image<double>* weightIm,
    double sigmaP, double psfap, int psfOrder,
    PsfLog& log, BVec& psf, double& nu, long& flag);
void measureSinglePsf1(
    Position cen, const Image<double>& im, double sky,
    const Transformation& trans,
    double noise, double gain, const Image<double>* weightIm,
    double sigmaP, double psfap, int psfOrder,
    PsfLog& log, BVec& psf, double& nu, long& flag);

#endif
