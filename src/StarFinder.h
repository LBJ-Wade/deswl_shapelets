
#ifndef StarFinder_H
#define StarFinder_H

#include <vector>
#include "ConfigFile.h"
#include "Function2D.h"
#include "PotentialStar.h"

class StarFinderException : public std::runtime_error {
  public:
    StarFinderException(std::string m) : std::runtime_error(m) {};
};

class StarFinder
{
  public :
    StarFinder(const ConfigFile& params, std::string key_prefix);

    std::vector<PotentialStar*> FindStars(
	std::vector<PotentialStar*>& allobj);

    void FindMinMax(const std::vector<PotentialStar*>& list,
	double *min, double *max, const Function2D<double>& f);

    void OutlierReject(std::vector<PotentialStar*>& list, 
	double nsigma, double minsigma, const Function2D<double>& f);

    std::vector<PotentialStar*> GetPeakList(
	const std::vector<PotentialStar*>& objlist,
	double binsize, double minsize, double maxsize,
	size_t startn, int miniter, double magstep, double maxsignifratio,
	bool firstpass,
	const Function2D<double>& f);

    void FitStellarSizes(Function2D<double> *f, size_t order, 
	double sigclip,
	const std::vector<PotentialStar*>& starlist, double *outsigma);

    void RoughlyFitBrightStars(const std::vector<PotentialStar*>& objlist,
	Function2D<double> *f,double *outsigma);


    double minsize;         // The min and max size to consider
    double maxsize;
    bool logsize;           // 1 if sizes are already log(size)
    double minmag;          // The min and max magnitude to consider
    double maxmag;
    double maxoutmag;
    size_t ndivx;           // Divide the image into nx by ny subsections
    size_t ndivy;

    // Parameters when finding stars in each subdivision
    double startn1;         // # objects to start with on first pass 
    double starfrac;        // Of these, how many are probably stars?
    double magstep1;        // Step size in mag for each successive pass
    double reject1;         // nsigma rejection.  (actually n x quartile size)
    double maxratio1;       // max ratio of count in valley / count in peak
    double binsize1;        // binsize of histogram 
    size_t okvalcount;      // if valley count <= this, ok no matter what peak is
    size_t miniter1;        // min number of mag steps to take
    double maxrms;          // in initial linear fit, max rms to allow 

    // Parameters for last pass of whole image
    double startn2;         // # objects to start with
    double magstep2;        // Step size in mag
    double minbinsize;      // Min width of histogram bins
    double reject2;         // n quartile rejection
    double purityratio;     // max ratio of count in valley / count in peak
    size_t miniter2;        // min number of mag steps to take

    // Parameters for fitting results of subsections
    size_t starsperbin;     // How many stars per subsection?
    size_t fitorder;        // order of fit for size(x,y)
    double fitsigclip;      // nsigma rejection in this fit
    size_t maxrefititer;    // max number of times to refit whole image

#if 0
    // Parameters for reading input file
    size_t xcol;            // which column for x,y,m,ixx,iyy,errcode
    size_t ycol;
    size_t mcol;
    size_t scol1;           // scol1,2 are ixx,iyy usually. 
    size_t scol2;           // if scol2=0, only one column for size
    size_t ecol;            // if ecol = 0, no errcode
    int okerrcode;          // errcodes to allow
#endif

};

#endif
