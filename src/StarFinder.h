#ifndef _star_finder_h
#define _star_finder_h

#include <string>
#include "unistd.h"
#include <string>
using std::string;
#include "Legendre2D.h"
#include "Bounds.h"
#include "PotentialStar.h"
#include "Histogram.h"
#include "ConfigFile.h"

//#include "Errors.hpp"

#include <algorithm>
#include <fstream>
using std::ifstream;
using std::ofstream;
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <functional>
#include <vector>
using std::vector;
#include <sstream>
using std::istringstream;

using std::abs;

/*
   template <class T> inline const T& MIN(const T& a, const T& b)
   { return (a < b) ? a : b; }
   template <class T> inline const T& MAX(const T& a, const T& b)
   { return (a > b) ? a : b; }
   */

class StarFinderException : public std::runtime_error {
  public:
    StarFinderException(const char* m) : std::runtime_error(m) {};
    StarFinderException(std::string m) : std::runtime_error(m.c_str()) {};
    StarFinderException(char* m) : std::runtime_error((const char*) m) {};
};



class StarFinder
{
  public:
    StarFinder() {};
    StarFinder(std::string configfile);
    StarFinder(ConfigFile& configfile);
    ~StarFinder() {};

    void Init(std::string configfile);
    void Init(ConfigFile& configfile);
    void LoadConfig();
    void CopyConfig(ConfigFile& configfile);

    void RunFindStars(
	vector<int>& flags,
	vector<int>& size_flags,
	vector<float>& x,
	vector<float>& y,
	vector<double>& sigma,
	vector<float>& mag,
	vector<int>& starflags);


    void TestConfig();

    std::vector<string> RequiredConfigFields();


    vector<PotentialStar*> FindStars(vector<PotentialStar*>& allobj);
    void FindMinMax(const vector<PotentialStar*>& list, 
	double *min, double *max, const Function2D<double>& f);
    void OutlierReject(vector<PotentialStar*>& list, 
	double nsigma, double minsigma, const Function2D<double>& f);
    vector<PotentialStar*> 
      GetPeakList(const vector<PotentialStar*>& objlist,
	  double binsize, double minsize, double maxsize,
	  size_t startn, int miniter, double magstep, 
	  double maxsignifratio,
	  bool firstpass,
	  const Function2D<double>& f);
    void FitStellarSizes(Function2D<double> *f, size_t order, 
	double sigclip,
	const vector<PotentialStar*>& starlist, double *outsigma);
    void RoughlyFitBrightStars(const vector<PotentialStar*>& objlist,
	Function2D<double> *f,double *outsigma);



  protected:

    std::string mCFile;
    ConfigFile mConfigFile;

    // Default values for the parameters used in this program

    // These only used for initial trim
    double mMinsize;       // The min and max size to consider
    double mMaxsize;
    double mMinmag;       // The min and max magnitude to consider
    double mMaxmag;
    double mMaxoutmag;

    // These are used in search methods.
    size_t mNdivx;           // Divide the image into nx by ny subsections
    size_t mNdivy;

    // Parameters when finding stars in each subdivision
    double mStartn1;        // # objects to start with on first pass 
    double mStarfrac;        // Of these, how many are probably stars?
    double mMagstep1;     // Step size in mag for each successive pass
    double mReject1;       // nsigma rejection.  (actually n x quartile size)
    double mMaxratio1;    // max ratio of count in valley / count in peak
    double mBinsize1;      // binsize of histogram 
    size_t mOkvalcount;      // if valley count <= this, ok no matter what peak is
    size_t mMiniter1;        // min number of mag steps to take
    double mMaxrms;       // in initial linear fit, max rms to allow 

    // Parameters for last pass of whole image
    double mStartn2;       // # objects to start with
    double mMagstep2;     // Step size in mag
    double mMinbinsize;	  // Min width of histogram bins
    double mReject2;       // n quartile rejection
    double mPurityratio;  // max ratio of count in valley / count in peak
    size_t mMiniter2;        // min number of mag steps to take

    // Parameters for fitting results of subsections
    size_t mStarsperbin;    // How many stars per subsection?
    size_t mFitorder;        // order of fit for size(x,y)
    double mFitsigclip;    // nsigma rejection in this fit
    size_t mMaxrefititer;   // max number of times to refit whole image


};

#endif
