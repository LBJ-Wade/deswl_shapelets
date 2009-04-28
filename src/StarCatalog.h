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
void CalcSigma(
    double& sigma, // Initial value -- use <=0 if no initial guess
    const Image<double>& im, const Position& pos, double sky, 
    double noise, double gain, const Image<double>* weight_im,
    const Transformation& trans, double psfap, long& flag);

class StarCatalog
{
  public:
    // Make from incat
    // fs_prefix is the prefix of the keywords for the 
    // parameters used by the findstars algorithm.
    StarCatalog(const InputCatalog& incat,
	ConfigFile& configfile, std::string fs_prefix = "stars_");

    // Read in from file
    StarCatalog(ConfigFile& configfile, std::string fs_prefix = "stars_");

    size_t size() const { return id.size(); }
    void Read();
    void Write() const;

    void ReadFits(std::string file);
    void ReadFitsOld(std::string file);
    void ReadAscii(std::string file, std::string delim = "  ");
    void WriteFitsOld(std::string file) const;
    void WriteFits(std::string file) const;
    void WriteAscii(std::string file, std::string delim = "  ") const;

    void CalcSizes(const Image<double>& im, 
	const Image<double>* weight_im, const Transformation& trans);

    int FindStars(FindStarsLog& log);

    std::vector<long> id;
    std::vector<Position> pos;
    std::vector<double> sky;
    std::vector<double> noise;
    std::vector<long> flags;

    std::vector<float> mag;
    std::vector<double> objsize;
    std::vector<int> isastar;

  private :

    const ConfigFile& params;
    std::string prefix;
};

#endif
