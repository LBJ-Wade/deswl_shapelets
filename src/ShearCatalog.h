#ifndef ShearCatalog_H
#define ShearCatalog_H

#include <vector>
#include "BVec.h"
#include "Bounds.h"
#include "ConfigFile.h"
#include "Name.h"
#include "TMV_Small.h"
#include "Transformation.h"
#include "InputCatalog.h"
#include "Image.h"
#include "Log.h"
#include "TimeVars.h"
#include "FittedPSF.h"

class ShearCatalog 
{
  public :

    // Make from incat, trans
    ShearCatalog(const InputCatalog& incat, const Transformation& trans,
	const ConfigFile& _params);

    // Read from file
    ShearCatalog(const ConfigFile& _params);

    size_t size() const { return id.size(); }
    void Read();
    void Write() const;

    void ReadFits(std::string file);
    void ReadAscii(std::string file, std::string delim = "  ");
    void WriteFits(std::string file) const;
    void WriteAscii(std::string file, std::string delim = "  ") const;

    int MeasureShears(const Image<double>& im,
	const Image<double>* weight_im, const Transformation& trans,
	const FittedPSF& fitpsf, ShearLog& log);

    std::vector<long> id;
    std::vector<Position> pos;
    std::vector<double> sky;
    std::vector<double> noise;
    std::vector<long> flags;
    std::vector<Position> skypos;
    std::vector<std::complex<double> > shear;
    std::vector<double> nu;
    std::vector<tmv::SmallMatrix<double,2,2> > cov;
    std::vector<BVec> shape;

  private :

    const ConfigFile& params;

};

void MeasureSingleShear(
    Position cen, const Image<double>& im, double sky,
    const Transformation& trans, const FittedPSF& fitpsf,
    double noise, double gain, const Image<double>* weight_im, 
    double gal_aperture, double max_aperture,
    int gal_order, int gal_order2,
    double f_psf, double min_gal_size, bool desqa,
    OverallFitTimes* times, ShearLog& log,
    std::complex<double>& shear, 
    tmv::SmallMatrix<double,2,2>& shearcov, BVec& shapelet,
    long& flag);

void MeasureSingleShear1(
    Position cen, const Image<double>& im, double sky,
    const Transformation& trans, const std::vector<BVec>& psf,
    double noise, double gain, const Image<double>* weight_im, 
    double gal_aperture, double max_aperture,
    int gal_order, int gal_order2,
    double f_psf, double min_gal_size, bool desqa,
    OverallFitTimes* times, ShearLog& log,
    std::complex<double>& shear, 
    tmv::SmallMatrix<double,2,2>& shearcov, BVec& shapelet,
    long& flag);

#endif
