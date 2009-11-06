#ifndef PSFCatalog_H
#define PSFCatalog_H

#include <sstream>
#include <vector>
#include "BVec.h"
#include "Bounds.h"
#include "ConfigFile.h"
#include "StarCatalog.h"
#include "Log.h"
#include "Image.h"
#include "Transformation.h"

class PSFCatalog 
{
  public :

    // Make from starcat
    PSFCatalog(const StarCatalog& starcat, const ConfigFile& _params);

    // Read from file
    PSFCatalog(const ConfigFile& _params);
	// this one we don't have to deal with root= stuff
    PSFCatalog(const ConfigFile& _params, std::string file);

    size_t size() const { return id.size(); }
    void Read();
    void Read(std::string file);
    void Write() const;

    void ReadFits(std::string file);
    void ReadAscii(std::string file, std::string delim = "  ");
    void WriteFits(std::string file) const;
    //void WriteFitsOld(std::string file) const;
    void WriteAscii(std::string file, std::string delim = "  ") const;

    double EstimateSigma(const Image<double>& im,
	const Image<double>* weight_im, const Transformation& trans);

    int MeasurePSF(const Image<double>& im, const Image<double>* weight_im,
	const Transformation& trans, double sigma_p, PSFLog& log);

    std::vector<long> id;
    std::vector<Position> pos;
    std::vector<double> sky;
    std::vector<double> noise;
    std::vector<long> flags;
    std::vector<double> nu;
    std::vector<BVec> psf;

  private :

    const ConfigFile& params;

};

void MeasureSinglePSF(
    Position cen, const Image<double>& im, double sky,
    const Transformation& trans,
    double noise, double gain, const Image<double>* weight_im,
    double sigma_p, double psfap, int psforder, bool desqa,
    PSFLog& log, BVec& psf, double& nu, long& flag);
void MeasureSinglePSF1(
    Position cen, const Image<double>& im, double sky,
    const Transformation& trans,
    double noise, double gain, const Image<double>* weight_im,
    double sigma_p, double psfap, int psforder, bool desqa,
    PSFLog& log, BVec& psf, double& nu, long& flag);

#endif
