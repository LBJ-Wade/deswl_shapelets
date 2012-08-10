
#ifndef StarFinder_H
#define StarFinder_H

#include <vector>
#include <stdexcept>
#include "dbg.h"
#include "ConfigFile.h"
#include "Function2D.h"
#include "PotentialStar.h"

class StarFinderException : 
    public std::runtime_error 
{
public:
    StarFinderException(const std::string& m) : std::runtime_error(m) {};
};

class StarFinder
{
public :
    StarFinder(const ConfigFile& params, std::string key_prefix);

    std::vector<PotentialStar*> findStars(
        std::vector<PotentialStar*>& allobj);

    void findMinMax(
        const std::vector<PotentialStar*>& list,
        double *min, double *max, const Function2D& f);

    void rejectOutliers(
        std::vector<PotentialStar*>& list, 
        double nsigma, double min_sigma, const Function2D& f);

    std::vector<PotentialStar*> getPeakList(
        const std::vector<PotentialStar*>& obj_list,
        double binSize, double min_size, double max_size,
        int nstart, int min_iter, double mag_step, double maxSignif_ratio,
        bool first_pass, const Function2D& f);

    void fitStellarSizes(
        Function2D *f, int order, double sig_clip,
        const std::vector<PotentialStar*>& star_list, double *out_sigma);

    void roughlyFitBrightStars(
        const std::vector<PotentialStar*>& obj_list,
        Function2D *f, double *out_sigma);

    void setParams(
        const ConfigFile& params, std::string key_prefix,
        bool mustExist=false);

    bool isOkSize(const double size)
    { return size >= _min_size && size <= _max_size; }

    bool isOkMag(const double mag)
    { return mag >= _min_mag && mag <= _max_mag; }

    bool isOkNu(const double nu)
    { return nu >= _min_nu; }

    bool isOkOutputMag(const double mag)
    { return mag >= _min_mag && mag <= _max_out_mag; }

    bool isOkSg(const double sg)
    { return sg >= _min_sg && sg <= _max_sg; }

    bool isOkSgMag(const double mag)
    { return mag >= _min_sg_mag && mag <= _max_sg_mag; }

    double convertToLogSize(const double size)
    { return _is_size_log ? size : std::log(size); }

    double getMinSize() const { return _min_size; }
    double getMaxSize() const { return _max_size; }
    double getMinMag() const { return _min_mag; }
    double getMaxMag() const { return _max_mag; }
    double getMinSg() const { return _min_sg; }
    double getMaxSg() const { return _max_sg; }
    double getMinNu() const { return _min_nu; }
    double getMinSgFrac() const { return _min_sg_frac; }

private :

    double _min_size;       // The min and max size to consider
    double _max_size;
    bool _is_size_log;       // true if sizes are already log(size)
    double _min_mag;        // The min and max magnitude to consider
    double _max_mag;
    double _min_nu;
    double _max_out_mag;
    double _min_sg;         // The min and max star-galaxy value to consider
    double _max_sg;
    double _min_sg_mag;      // The min and max mag for initial star selection
    double _max_sg_mag;
    double _min_sg_frac;     // The minimum fraction of objects in the star-galaxy range
    int _ndivx;         // Divide the image into nx by ny subsections
    int _ndivy;

    // Parameters when finding stars in each subdivision
    double _nstart1;       // # objects to start with on first pass 
    double _star_frac;      // Of these, how many are probably stars?
    double _mag_step1;      // Step size in mag for each successive pass
    double _reject1;       // nsigma rejection.  (actually n x quartile size)
    double _max_ratio1;     // max ratio of count in valley / count in peak
    double _bin_size1;      // binsize of histogram 
    int _ok_val_count;    // if valley count <= this, ok no matter what peak is
    int _min_iter1;      // min number of mag steps to take
    double _max_rms;        // in initial linear fit, max rms to allow 

    // Parameters for last pass of whole image
    double _nstart2;       // # objects to start with
    double _mag_step2;      // Step size in mag
    double _min_bin_size;    // Min width of histogram bins
    double _reject2;       // n quartile rejection
    double _purity_ratio;   // max ratio of count in valley / count in peak
    int _min_iter2;      // min number of mag steps to take

    // Parameters for fitting results of subsections
    int _stars_per_bin;   // How many stars per subsection?
    int _fit_order;      // order of fit for size(x,y)
    double _fit_sig_clip;    // nsigma rejection in this fit
    int _max_refit_iter;  // max number of times to refit whole image

    bool _des_qa;  // Output warnings in DES QA format

};

#endif
