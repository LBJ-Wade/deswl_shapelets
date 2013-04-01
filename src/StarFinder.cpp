
#include <algorithm>
#include <fstream>
#include <iostream>
#include <functional>
#include <vector>
#include <sstream>
#include <string>
#include <string>
#include <valarray>
#include "unistd.h"

#include "StarCatalog.h"
#include "StarFinder.h"
#include "Legendre2D.h"
#include "Bounds.h"
#include "PotentialStar.h"
#include "Histogram.h"
#include "ConfigFile.h"
#include "fspd.h"

#define SFKeyAssign(var,key) \
    do { \
        if (mustexist || params.keyExists(keyPrefix + key)) { \
            var = params[ keyPrefix + key ];\
        } \
        xdbg << "SFAssign " key " = " << var; \
        xdbg << " from parameter "<<(keyPrefix + key)<<std::endl; \
    } while (false) 

void StarFinder::setParams(
    const ConfigFile& params, std::string keyPrefix, bool mustexist)
{
    SFKeyAssign(_min_size,"minsize");
    SFKeyAssign(_max_size,"maxsize");
    SFKeyAssign(_min_sg,"minsg");
    SFKeyAssign(_max_sg,"maxsg");
    SFKeyAssign(_min_sg_mag,"minsgmag");
    SFKeyAssign(_max_sg_mag,"maxsgmag");
    SFKeyAssign(_min_sg_frac,"minsgfrac");
    SFKeyAssign(_is_size_log,"logsize");



    SFKeyAssign(_min_mag,"minmag");
    SFKeyAssign(_max_mag,"maxmag");
    SFKeyAssign(_max_out_mag,"maxoutmag");
    SFKeyAssign(_min_nu,"minnu");
    SFKeyAssign(_ndivx,"ndivx");
    SFKeyAssign(_ndivy,"ndivy");

    SFKeyAssign(_nstart1,"startn1");
    SFKeyAssign(_star_frac,"starfrac");
    SFKeyAssign(_mag_step1,"magstep1");
    SFKeyAssign(_reject1,"reject1");
    SFKeyAssign(_max_ratio1,"maxratio1");
    SFKeyAssign(_bin_size1,"binsize1");
    SFKeyAssign(_ok_val_count,"okvalcount");
    SFKeyAssign(_min_iter1,"miniter1");
    SFKeyAssign(_max_rms,"maxrms");

    SFKeyAssign(_nstart2,"startn2");
    SFKeyAssign(_mag_step2,"magstep2");
    SFKeyAssign(_min_bin_size,"minbinsize");
    SFKeyAssign(_reject2,"reject2");
    SFKeyAssign(_purity_ratio,"purityratio");
    SFKeyAssign(_min_iter2,"miniter2");

    SFKeyAssign(_stars_per_bin,"starsperbin");
    SFKeyAssign(_fit_order,"fitorder");
    SFKeyAssign(_fit_sig_clip,"fitsigclip");
    SFKeyAssign(_max_refit_iter,"maxrefititer");
}

StarFinder::StarFinder(const ConfigFile& params, std::string keyPrefix)
{
    std::string fspd((const char*)findstars_params_default,
                     findstars_params_default_len);
    std::istringstream is(fspd);
    ConfigFile defaultParams;
    defaultParams.setDelimiter("\t");
    defaultParams.setComment("#");
    defaultParams.read(is);

    setParams(defaultParams,"",true);
    setParams(params,keyPrefix);

    _des_qa = params.read("des_qa",false); 
}

#undef SFKeyAssign

std::vector<PotentialStar*> StarFinder::findStars(
    std::vector<PotentialStar*>& allobj)
{
    dbg<<"starting FindStars, allobj.size() = "<<allobj.size()<<std::endl;

    // sort allobj by magnitude
    std::sort(allobj.begin(),allobj.end(),std::mem_fun(
            &PotentialStar::isBrighterThan));
    dbg<<"sorted allobj\n";


    // First stage:
    // Split area into sections
    // For each secion find as many stars as possible.
    // Use these stars to fit the variation in rsq across the image.

    // Make bounds of whole region called total_bounds
    Bounds total_bounds;
    const int nobj = allobj.size();
    for(int k=0;k<nobj;++k) total_bounds += allobj[k]->getPos();
    dbg<<"totalbounds = \n"<<total_bounds<<std::endl;

    // boundsar is the 3x3 (ndivx x ndivy) array of quadrant bounds
    // call it qbounds even though it won't be quadrants unless 2x2
    std::vector<Bounds> qbounds = total_bounds.divide(_ndivx,_ndivy);
    xdbg<<"made qbounds\n";

    // probstars will be our first pass list of probable stars
    std::vector<PotentialStar*> probstars;

    // For each section, find the stars and add to probstars
    const int nsection = qbounds.size();

    double frac_nobj=0;
    for (int i=0; i<nobj; ++i) {
        if(isOkSg(allobj[i]->getSg()) && 
           isOkSgMag(allobj[i]->getMag())) frac_nobj++;
    }
    dbg<<"  Found "<<frac_nobj<<" with star-galaxy cuts  "<<std::endl;
    frac_nobj/=nobj;

    bool ignore_sg_cut=false;
    if( frac_nobj < _min_sg_frac ) {
        dbg<<"  Fraction of objects in star-galaxy range too small: "<<
            frac_nobj<<std::endl;
        dbg<<"  Ignoring star-galaxy cut selection"<<std::endl;
        ignore_sg_cut=true;
    }


    if(!ignore_sg_cut) {   
        for(int k=0;k<nobj;++k) {
            if (isOkSg(allobj[k]->getSg()) && isOkSgMag(allobj[k]->getMag())) {
                probstars.insert(probstars.end(),allobj[k]);
            }
        }

        // dummy fit to pass to rejection function
        Legendre2D fit(total_bounds);

        rejectOutliersIter(probstars,_reject1,0.1,fit,true,1e-4,8);
        dbg<<"rejected outliers, now have "<<probstars.size()<<" stars\n";
        std::sort(probstars.begin(),probstars.end(),
            std::mem_fun(&PotentialStar::isBrighterThan));

        dbg<<"Current stars mag/size:"<<std::endl;
        for(size_t k=0;k<probstars.size();++k) {
            dbg<<probstars[k]->getMag()<<" "<<probstars[k]->getSize()<<std::endl;
        }
    } else {
        for(int i=0;i<nsection;++i) {
            dbg<<"i = "<<i<<": bounds = "<<qbounds[i]<<std::endl;

            // someobj are the objects in this section
            // Note that someobj is automatically sorted by magnitude, since
            // allobj was sorted.
            std::vector<PotentialStar*> someobj;
            for(int k=0;k<nobj;++k) {
                if (qbounds[i].includes(allobj[k]->getPos())) 
                    someobj.push_back(allobj[k]);
            }
            dbg<<"added "<<someobj.size()<<" obj\n";

            // Does a really quick and dirty fit to the bright stars
            // Basically it takes the 10 smallest of the 50 brightest objects,
            // finds the peakiest 5, then fits their sizes to a 1st order function.
            // It also gives us a rough value for the sigma
            Legendre2D flinear(qbounds[i]);
            double sigma;
            roughlyFitBrightStars(someobj,&flinear,&sigma);
            dbg<<"fit bright stars: sigma = "<<sigma<<std::endl;

            // Calculate the min and max values of the (adjusted) sizes
            double min_size,max_size;
            findMinMax(someobj,&min_size,&max_size,flinear);
            dbg<<"min,max = "<<min_size<<','<<max_size<<std::endl;

            // Find the objects clustered around the stellar peak.
            std::vector<PotentialStar*> qpeak_list =
                getPeakList(
                    someobj,_bin_size1,min_size,max_size,
                    int(_nstart1*someobj.size()),_min_iter1,_mag_step1,_max_ratio1,
                    true,flinear);
            const int npeak = qpeak_list.size();
            dbg<<"peaklist has "<<npeak<<" stars\n";

            // Remove outliers using a median,percentile rejection scheme.
            // _bin_size1/2. is the minimum value of "sigma".
            rejectOutliers(qpeak_list,_reject1,_bin_size1/2.,flinear);
            dbg<<"rejected outliers, now have "<<npeak<<" stars\n";

            // Use at most 10 (stars_per_bin) stars per region to prevent one region 
            // from dominating the fit.  Use the 10 brightest stars to prevent being
            // position biased (as one would if it were based on size
            int nstars_expected = int(_star_frac * someobj.size());
            if (npeak < nstars_expected) {
                if (npeak < int(0.2 * nstars_expected)) {
                    if (_des_qa) {
                        std::cout<<"STATUS3BEG Warning: Only "<<
                            qpeak_list.size()<<" stars found in section "<<
                            i<<". STATUS3END"<<std::endl;
                    }
                    dbg<<"Warning: only "<<qpeak_list.size()<<
                        " stars found in section "<<i<<
                        "  "<<qbounds[i]<<std::endl;
                }
                probstars.insert(probstars.end(),qpeak_list.begin(),qpeak_list.end());
            } else {
                std::sort(qpeak_list.begin(),qpeak_list.end(),
                          std::mem_fun(&PotentialStar::isBrighterThan));
                probstars.insert(probstars.end(),qpeak_list.begin(),
                                 qpeak_list.begin()+nstars_expected);
            }
            dbg<<"added to probstars\n";
        }
        xdbg<<"done qbounds loop\n";
    }
    int nstars = probstars.size();
    dbg<<"nstars = "<<nstars<<std::endl;

    // Now we have a first estimate of which objects are stars.
    // Fit a quadratic function to them to characterize the variation in size
    Legendre2D f(total_bounds);
    double sigma;
    fitStellarSizes(&f,_fit_order,_fit_sig_clip,probstars,&sigma);

    // Second stage:
    // Use the fitted function for the size we just got to adjust
    // the measured sizes according to their position in the image.
    // Make the histogram using measured - predicted sizes (still log
    // size actually) which should then have the peak very close to 0.
    // Set the bin size for this new histogram to be the rms scatter of the
    // stellar peak from pass 1.
    // This time step by 0.25 mag (mag_step2).

    // Find the values of min_size,max_size for the whole thing with the new f
    double min_size,max_size;
    findMinMax(allobj,&min_size,&max_size,f);
    dbg<<"new minmax = "<<min_size<<','<<max_size<<std::endl;

    // Find the objects clustered around the stellar peak.
    // Note that the bin_size is a set fraction of sigma (0.5), so this
    // should make the stellar peak clearer.
    // Also, the f at the end means the functional fitted size will be
    // subtracted off before adding to the histogram.
    probstars = getPeakList(
        allobj,0.5*sigma,min_size,max_size,
        int(_nstart2*allobj.size()),_min_iter2,_mag_step2,_purity_ratio,false,f);
    nstars = probstars.size();
    dbg<<"probstars has "<<nstars<<" stars\n";

    // Remove outliers using a iterative median,percentile rejection scheme.
    rejectOutliersIter(probstars,_reject2,0.1,f,false);
    dbg<<"rejected outliers, now have "<<probstars.size()<<" stars\n";
    std::sort(probstars.begin(),probstars.end(),
              std::mem_fun(&PotentialStar::isBrighterThan));
    

    // If you get a bad fit the first time through, it can take a 
    // few passes to fix it.
    // And always do at least one refit.
    bool refit = true;
    for(int iter=0;refit && iter<_max_refit_iter;++iter) {
        dbg<<"starting refit\n"<<std::endl;
        refit = false;
        std::vector<std::vector<PotentialStar*> > stars_list(nsection);

        // Add each star to the appropriate sublist
        for(int k=0;k<nstars;++k) {
            for(int i=0;i<nsection;++i) {
                if(qbounds[i].includes(probstars[k]->getPos()))
                    stars_list[i].push_back(probstars[k]);
            }
        }

        // Make fit_list, use sg cuts if enough objects otherwise use the 
        // list of the 10 brightest stars per section
        std::vector<PotentialStar*> fit_list;

        if(!ignore_sg_cut) {
            for(int k=0;k<nobj;++k) {

                if (isOkSg(allobj[k]->getSg()) && isOkSgMag(allobj[k]->getMag())) {
                    fit_list.insert(fit_list.end(),allobj[k]);
                }
            }
            rejectOutliersIter(fit_list,_reject1,0.1,f,false);
            dbg<<"rejected outliers, now have "<<fit_list.size()<<" stars\n";
            std::sort(fit_list.begin(),fit_list.end(),
                      std::mem_fun(&PotentialStar::isBrighterThan));
            


        } else {

            dbg<<"  Fraction of objects in star-galaxy range too small "<<std::endl;
            for(int i=0;i<nsection;++i) {
                // if there are still < 10 stars, give a warning, and
                // just add all the stars to fit_list
                if (int(stars_list[i].size()) < _stars_per_bin) {
                    if (_des_qa) {
                        std::cout<<"STATUS3BEG Warning: Only "<<
                            stars_list[i].size()<<" stars in section "<<i<<
                            ". STATUS3END"<<std::endl;
                    }
                    dbg<<"Warning: only "<<stars_list[i].size()<<
                        " stars in section ";
                    dbg<<i<<"  "<<qbounds[i]<<std::endl;
                    fit_list.insert(fit_list.end(),stars_list[i].begin(),
                                    stars_list[i].end());
                    refit = true;
                } else {
                    // sort the sublist by magnitude
                    std::sort(stars_list[i].begin(),stars_list[i].end(),
                              std::mem_fun(&PotentialStar::isBrighterThan));

                    // add the brightest 10 to fit_list
                    fit_list.insert(fit_list.end(),stars_list[i].begin(),
                                    stars_list[i].begin()+_stars_per_bin);
                }
            }
        }

        // Do all the same stuff as before (refit f, get new min,max,
        // find the peak_list again, and reject the outliers)
        nstars = fit_list.size();
        fitStellarSizes(&f,_fit_order,_fit_sig_clip,fit_list,&sigma);
        findMinMax(allobj,&min_size,&max_size,f);
        dbg<<"new minmax = "<<min_size<<','<<max_size<<std::endl;
        probstars = getPeakList(
            allobj,0.5*sigma,min_size,max_size,
            int(_nstart2*allobj.size()),_min_iter2,_mag_step2,
            _purity_ratio,false,f);
        dbg<<"probstars has "<<probstars.size()<<" stars\n";
        rejectOutliersIter(probstars,_reject2,0.1,f,false);
        dbg<<"rejected outliers, now have "<<probstars.size()<<" stars\n";
        if (int(probstars.size()) < nstars) {
            refit = true;
            dbg<<"fewer than "<<nstars<<" - so refit\n";
        }
        nstars = probstars.size();
    }

    dbg<<"done FindStars\n";

    for(int i=0;i<nstars;++i) {
        xdbg<<"stars["<<i<<"]: size = "<<probstars[i]->getSize();
        xdbg<<", size - f(pos) = ";
        xdbg<<probstars[i]->getSize()-f(probstars[i]->getPos())<<std::endl;
        xdbg<<probstars[i]->getLine()<<std::endl;
    }
    return probstars;
}

void StarFinder::findMinMax(
    const std::vector<PotentialStar*>& list, 
    double *min, double *max, const Function2D& f)
{
    double min1,max1;
    min1 = max1 = list[0]->getSize()-f(list[0]->getPos());
    const int nstars = list.size();
    for(int k=1;k<nstars;++k) {
        double size = list[k]->getSize() - f(list[k]->getPos());
        if (size > max1) max1 = size;
        if (size < min1) min1 = size;
    }
    *min = min1;
    *max = max1;
}

void StarFinder::rejectOutliers(
    std::vector<PotentialStar*>& list, 
    double nsigma, double min_sigma, const Function2D& f)
{
    // This rejects outliers by finding the median and the quartile 
    // values, and calls sigma the average of the two quartile deviations.
    // "nsigma" is then how many of these "sigma" away from the median
    //  to consider something an outlier.

    const int nstars = list.size();
    if (nstars <= 4) return;
    // Find median, 1st and 3rd quartile stars;
    std::vector<PotentialStar*> modif_list(nstars);
    for(int k=0;k<nstars;++k) {
        modif_list[k] = new PotentialStar(*list[k]);
        double newSize = list[k]->getSize() - f(list[k]->getPos());
        modif_list[k]->setSize(newSize);
    }
    std::sort(modif_list.begin(),modif_list.end(),
              std::mem_fun(&PotentialStar::isSmallerThan));
    PotentialStar* mStar = modif_list[modif_list.size()/2];
    PotentialStar* q1Star = modif_list[modif_list.size()/4];
    PotentialStar* q3Star = modif_list[modif_list.size()*3/4];
    double median = mStar->getSize();
    double q1 = q1Star->getSize();
    double q3 = q3Star->getSize();
    double sigma = std::max(q3-median,median-q1);
    sigma = std::max(sigma,min_sigma);
    xdbg<<"q1,m,q3 = "<<q1<<" , "<<median<<" , "<<q3<<std::endl;
    xdbg<<"sigma = "<<sigma<<", nsig * sigma = "<<nsigma*sigma<<std::endl;
    // Remove elements which std::abs(x->getSize()-median) > nsigma*sigma
    int j=0;
    for(int k=0; k<nstars;++k) {
        if (std::abs(modif_list[k]->getSize()-median)>nsigma*sigma) {
            xdbg<<"k = "<<k<<", size = "<<modif_list[k]->getSize();
            xdbg<<" might be an outlier (median = "<<median<<")\n";
        }
        list[j] = list[k];  // Note: update original list, not modified version.
        ++j;
    }
    if (j<nstars) list.erase(list.begin()+j,list.end());
    for(int k=0;k<nstars;++k) delete modif_list[k];
}




void StarFinder::rejectOutliersIter(
    std::vector<PotentialStar*>& list,
    double nSigma, double minSigma, const Function2D& f,
    bool use_size,double tol,int max_iter)
{
    // This rejects outliers by finding the median and the quartile
    // values, and calls sigma the average of the two quartile deviations.
    // "nSigma" is then how many of these "sigma" away from the median
    //  to consider something an outlier.


  if (list.size() <= 4) return;

  double old_sigma=1e30,old_median=1e30;
  double sigma=0,median=0;
  int iter=0,nremoved=0;
  bool pass;
  do {

      pass=true;
      const int nStars = list.size();
      if(iter>=1) {
          old_sigma=sigma;
          old_median=median;
      }

      std::vector<double> sizes;
      std::map<int,double> msizes;
      for(int k=0;k<nStars;++k) {
          double newsize;
           if(!use_size) {
               newsize = list[k]->getSize() - f(list[k]->getPos());
           } else {
               newsize = list[k]->getSize();
           }
           msizes[list[k]->getIndex()]=newsize;
           sizes.push_back(newsize);
      }

      // Find median, 1st and 3rd quartile stars;
      std::sort(sizes.begin(),sizes.end());

      median = sizes[sizes.size()/2];
      double q1 = sizes[sizes.size()/4];
      double q3 = sizes[sizes.size()*3/4];

      sigma = std::max(q3-median,median-q1);

      if(sigma > minSigma) {
          dbg<<"Sigma too small reducing by a factor of 2"<<std::endl;
          sigma/=2;

          // Don't allow loop to stop if sigma is too large
          // set previous values very high
          old_sigma=1e30;
          old_median=1e30;
          pass=false;

          // if we had to reduce two times in a row further reduce the sigma
          // to avoid infinite loop
          if(nremoved==10000) sigma/=2;
      }

      dbg<<"q1,m,q3 = "<<q1<<" , "<<median<<" , "<<q3<<std::endl;
      dbg<<"sigma = "<<sigma<<", nsig * sigma = "<<nSigma*sigma<<std::endl;

      int j=0;
      // Remove elements which std::abs(x->getSize()-median) > nSigma*sigma
      for(int k=0; k<nStars;++k) {
          if (std::abs(msizes[list[k]->getIndex()]-median)>nSigma*sigma) {
              xdbg<<"k = "<<k<<",index= "<<list[k]->getIndex()
                  <<"  size = "<<msizes[list[k]->getIndex()];
              xdbg<<" rejecting (median = "<<median<<")\n";
          } else {
              xdbg<<"k = "<<k<<",index= "<<list[k]->getIndex()
                  <<", size = "<<msizes[list[k]->getIndex()]
                  <<", mag = "<<list[k]->getMag()<<" keeping\n";
              list[j++]=list[k];
          }
      }
      if (j<nStars) list.erase(list.begin()+j,list.end());
      nremoved=nStars-j;
    
      // modify nremoved if sigma was too small so we continue looping
      if(!pass) nremoved=10000;
      dbg<<"Remove Iteration: "<<iter+1<<"  Sigma "<<old_sigma<<" "<<sigma
           <<"  Median "<<old_median<<" "<<median
           <<" removed "<<nremoved<<" stars"
           <<" max_iter "<<max_iter<<std::endl;
      ++iter;
    } while( (fabs(old_sigma-sigma)>tol &&
              fabs(old_median-median)>tol &&
              (iter-1)<max_iter && nremoved>0));
}

  
std::vector<PotentialStar*> StarFinder::getPeakList(
    const std::vector<PotentialStar*>& obj_list,
    double bin_size, double min_size, double max_size,
    int nstart, int min_iter, double mag_step, double max_signif_ratio,
    bool first_pass, const Function2D& f)
{
    if (bin_size < _min_bin_size) bin_size = _min_bin_size;

    // Make a histogram to find the stellar peak
    Histogram<PotentialStar*> hist(bin_size,min_size,max_size);

    // Start with the nstart brightest objects, then step up in
    // magnitude by mag_step at a time until histogram peak
    // gets dilluted.
    int k;
    Assert(nstart<=int(obj_list.size()));
    for(k=0;k<nstart;++k) {
        hist.add(obj_list[k]->getSize()-f(obj_list[k]->getPos()),obj_list[k]);
    }
    xdbg<<"added "<<nstart<<" objects\n";
    if(XDEBUG && dbgout) hist.print(*dbgout,min_size,max_size);

    // Get first estimates of the stellar and (first) galaxy peaks
    // along with the corresponding valleys.
    double peak1 = hist.findFirstPeakAfter(min_size);
    xdbg<<"peak1 = "<<peak1<<" ("<<hist[peak1]<<")\n";
    double valley1 = hist.findFirstValleyAfter(peak1);
    xdbg<<"valley1 = "<<valley1<<" ("<<hist[valley1]<<")\n";
    while (valley1 < 0.0) {
        // Since we subtract the fit, stars should be close to 0
        peak1 = hist.findFirstPeakAfter(valley1);
        xdbg<<"new peak1 = "<<peak1<<" ("<<hist[peak1]<<")\n";

        // Old code, did not do well with stupid no-galaxy images created
        // by DES pipeline.
        //valley1 = hist.findFirstValleyAfter(peak1);
        //xdbg<<"new valley1 = "<<valley1<<" ("<<hist[valley1]<<")\n";

        double new_valley1 = hist.findFirstValleyAfter(peak1);
        if (!(new_valley1 > valley1)) {
            std::ostringstream err;
            err<<"Couldn't find a stellar peak.  ";
            err<<"All the objects seem to be basically the same size.";
            throw StarFinderException(err.str());
        }
        valley1 = new_valley1;
        xdbg<<"new valley1 = "<<valley1<<" ("<<hist[valley1]<<")\n";


    }
    double peak2 = hist.findFirstPeakAfter(valley1);
    xdbg<<"peak2 = "<<peak2<<" ("<<hist[peak2]<<")\n";
    // Sometimes the stellar peak is horned (eg. 212), in which case, find
    // the next peak
    double nbins_for_horn = first_pass ? 2.5 : 4.5;
    if (peak2-peak1 < nbins_for_horn*bin_size && peak2+bin_size < max_size) {
        valley1 = hist.findFirstValleyAfter(peak2);
        peak2 = hist.findFirstPeakAfter(peak2+bin_size);
        xdbg<<"horn pattern: next peak2 = "<<peak2<<" ("<<hist[peak2]<<")\n";
    }
    double valley2 = hist.findFirstValleyAfter(peak2);
    xdbg<<"valley2 = "<<valley2<<" ("<<hist[valley2]<<")\n";
    double valley0 = hist.findFirstValleyBefore(peak1);
    if (!first_pass && valley0 > 0.) valley0 = hist.findFirstValleyBefore(0.);
    xdbg<<"valley0 = "<<valley0<<" ("<<hist[valley0]<<")\n";

    // Get the stars in the first peak for the initial value of peak_list
    // When we're all done it will be the list of objects at the stellar peak.
    std::vector<PotentialStar*> peak_list=hist.getRefsInRange(valley0,valley1);

    // Define the highest mag so far in list
    double highMag = obj_list[nstart-1]->getMag();

    // Loop at least min_iter times.  Done is set to false whenever
    // we've just found a new best peak.
    bool done=false;
    double prev_valley=valley1;
    double prev_peak=0.;
    for(int count = 0; count < min_iter || !done; ++count) {

        // Increment highMag
        highMag += mag_step;
        xdbg<<"count = "<<count<<", highmag = "<<highMag<<std::endl;

        // Add new objects to histogram.
        const int nobj = obj_list.size();
        for(;k<nobj && obj_list[k]->getMag()<highMag;++k) {
            hist.add(obj_list[k]->getSize()-f(obj_list[k]->getPos()),obj_list[k]);
        }
        xdbg<<"hist has "<<k<<" objects\n";
        if(XDEBUG && dbgout) hist.print(*dbgout,valley0,valley2);

        xdbg<<"valley0 = "<<valley0<<std::endl;
        // Find the peak on the stellar side
        double peak = hist.findFirstPeakAfter(valley0,!first_pass);
        valley0 = hist.findFirstValleyBefore(peak); // for next pass
        xdbg<<"done findpeak: "<<peak<<std::endl;

        // Find the valley 
        // The !first_pass allows poisson noise fluctuations for the final pass
        double valley = hist.findFirstValleyAfter(peak,!first_pass);
        xdbg<<"done findvalley: "<<valley<<std::endl;

        // If valley ends up negative, find next peak and valley
        if (valley < 0.0) {
            double nextPeak = hist.findFirstPeakAfter(valley,!first_pass);
            // if next peak is closer to 0 than current one, use new values.
            while (std::abs(nextPeak) < std::abs(peak)) {
                peak = nextPeak;
                xdbg<<"negative valley - new peak: "<<peak<<std::endl;
                valley = hist.findFirstValleyAfter(peak,!first_pass);
                xdbg<<"new valley: "<<valley<<std::endl;
                if (valley < 0.0) {
                    nextPeak = hist.findFirstPeakAfter(valley,!first_pass);
                } else {
                    break;
                }
            }
        }

        // For the first pass, be extra conservative and don't let the valley
        // position get much larger than a previous good value.
        // If the peak is also larger than prev_valley, something weird is
        // going on, so bail on this pass.
        if (first_pass && valley > prev_valley+bin_size) {
            valley = prev_valley;
            xdbg<<"moved valley to "<<valley<<std::endl;
            if (peak >= valley) { 
                done = true;
                continue; 
            }
        }

        // If new peak is larger than the previous valley, and the new peak 
        // isn't closer to 0, then reject it (probably added enough objects to
        // the valley to fill it in completely).
        if (peak > prev_valley && std::abs(peak) > std::abs(prev_peak)) {
            xdbg<<"New peak passed previous valley.\n";
            done = true;
            continue;
        }

        // The refined counts allow the bins to be centered anywhere near
        // the respective peak and valley and finds how many objects
        // would fall in the optimally centered bin.
        int peak_count = hist.getRefinedPeakCount(&peak);
        xdbg<<"peakcount = "<<peak_count<<" centered at "<<peak<<std::endl;
        xdbg<<"before refined: valley = "<<valley<<std::endl;
        int valley_count = hist.getRefinedValleyCount(&valley);
        xdbg<<"valleycount = "<<valley_count<<" centered at "<<valley<<std::endl;

        // Check for the horned pattern and use the next valley if it is there.
        // But only if the new valley is at least as good as the first one.
        //double nbins_for_horn = first_pass ? 1.5 : 3.5;
        nbins_for_horn = 1.5;
        if ((valley-peak < nbins_for_horn*bin_size) && 
            (hist.findFirstPeakAfter(valley) - valley) < nbins_for_horn*bin_size) {
            double new_valley = hist.findFirstValleyAfter(
                hist.findFirstPeakAfter(valley));
            int new_valley_count = hist.getRefinedValleyCount(&new_valley);
            if (((first_pass && new_valley_count <= valley_count) ||
                 (!first_pass && valley_count > 0 && new_valley_count == 0)) && 
                new_valley <= prev_valley + nbins_for_horn*bin_size) {

                xdbg<<"horn pattern detected.  new valley = "<<
                    new_valley<<std::endl;
                xdbg<<"new valley count = "<<new_valley_count<<std::endl;
                valley = new_valley;
                valley_count = new_valley_count;
            }
        }

        // The significance of the peak is the ratio of the number in the
        // valley to the number in the peak
        // So _small_ significances are good.
        // Has to be this way, since valley_count can be 0,
        // so peak_count/valley_count can be undefined.
        double signif = (double)valley_count/(double)peak_count;
        dbg<<"signif = "<<valley_count<<"/"<<peak_count<<" = "<<signif<<std::endl;

        // If the valley count is <= ok_val_count then drop the significance
        // to 0 regardless of the peak count, since sometimes the peak is fairly
        // broad so peak_count isn't high enough to get a good ratio.
        if (first_pass && valley_count <= _ok_val_count) {
            signif = 0.;
            dbg<<"reset signif to 0, since valleycount "<<valley_count<<" <= ";
            dbg<<_ok_val_count<<std::endl;
        }

        // If we have a new good peak, get all the stars in it for peak_list
        // Also make sure we repeat the loop by setting done = false
        if (signif <= max_signif_ratio) {
            dbg<<"signif = "<<signif<<" <= "<<max_signif_ratio<<std::endl;
            xdbg<<"good signif value\n";
            // Range is symmetrical about peakvalue with upper limit at valley
            double peak_start = std::min(2.*peak-valley,peak-bin_size*0.5);
            std::vector<PotentialStar*> temp = 
                hist.getRefsInRange(peak_start,valley);
            xdbg<<"get refs from "<<peak_start<<" to "<<valley<<std::endl;
            xdbg<<"got refs: "<<temp.size()<<std::endl;
            // If this list has fewer objects than a previous list,
            // and it doesn't just have a tighter range,
            // don't take it.  Surprisingly, this is actually possible, and it
            // means something weird happened
            if (temp.size() >= peak_list.size() || 
                (valley<prev_valley && valley>prev_peak) || 
                (peak > prev_peak) ) {

                if (temp.size() >= peak_list.size()) {
                    xdbg<<"tempsize = "<<temp.size()<<" >= "<<peak_list.size();
                } else if ((valley<prev_valley && valley>prev_peak)) {
                    xdbg<<"valley = "<<valley<<" < "<<prev_valley;
                } else {
                    xdbg<<"peak = "<<peak<<" > "<<prev_peak;
                }
                xdbg<<" so keep adding...\n";
                peak_list = temp;
                done = false;
                prev_valley = valley;
                prev_peak = peak;
            } else {
                xdbg<<"tempsize < "<<peak_list.size();
                xdbg<<" and peak,valley didn't contract, so stop adding now.\n";
                done=true;
            }
        } else {
            dbg<<"signif = "<<signif<<" > "<<max_signif_ratio<<std::endl;
            done=true; // maybe.  still loop if count < min_iter
        }

        // If we've already used all the stars, don't loop anymore
        if (k == nobj) {
            xdbg<<"no more to add\n";
            count = min_iter; done = true; 
        }
    } // for count
    xdbg<<"done finding peaklist\n";
    return peak_list;
}

void StarFinder::fitStellarSizes(
    Function2D *f, int order, double sig_clip,
    const std::vector<PotentialStar*>& star_list, double *out_sigma)
{
    // Make list of positions
    std::vector<Position> pos_list(star_list.size());
    std::transform(star_list.begin(),star_list.end(),pos_list.begin(),
                   std::mem_fun(&PotentialStar::getPos));

    // Make list of sizes
    std::vector<double> sizelist(star_list.size());
    std::transform(star_list.begin(),star_list.end(),sizelist.begin(),
                   std::mem_fun(&PotentialStar::getSize));

    // Use all stars in list (to start with anyway), so all true
    std::vector<bool> use_list(star_list.size(),true);

    if (int(star_list.size()) <= (order+1)*(order+2)/2) {
        std::ostringstream err;
        err<<"Not enough stars to do fit.  ";
        err<<"Increase stars_per_bin or decrease fit_order.";
        throw StarFinderException(err.str());
    }

    // Do the fit.
    double chisq;
    int dof;
    xdbg<<"before outlier fit\n";
    f->outlierFit(order,sig_clip,pos_list,sizelist,&use_list,0,&chisq,&dof,0);
    xdbg<<"after outlier fit\n";

    xdbg<<"chisq,dof,sigma = "<<chisq<<','<<dof<<','<<
        sqrt(chisq/dof)<<std::endl;
    *out_sigma = sqrt(chisq/dof);
}

void StarFinder::roughlyFitBrightStars(
    const std::vector<PotentialStar*>& obj_list,
    Function2D *f,double *out_sigma)
{
    // obj_list is already sorted by magnitude
    // make a new list with just the 50 brightest objects:
    std::vector<PotentialStar*> bright_list(
        obj_list.begin(),
        obj_list.begin()+int(_nstart1*obj_list.size()));

    // now sort this list by size
    std::sort(bright_list.begin(),bright_list.end(),
              std::mem_fun(&PotentialStar::isSmallerThan));

    xdbg<<"brightlist is:\n";
    const int twenty = std::min(20,int(bright_list.size()));
    for(int i=0; i<twenty; ++i) {
        xdbg<<bright_list[i]->getMag()<<" , "<<
            bright_list[i]->getSize()<<std::endl;
    }
    if (bright_list.size() > 20) {
        xdbg<<"...  (total of "<<bright_list.size()<<" elements)\n";
    }
    // Of the smallest 20, find the 5 that form the tightest peak
    // Originally I hardcoded this number as 5, but now it is calculated
    // from _star_frac.
    int five = int(0.5*_star_frac*bright_list.size());
    xdbg<<"'five' = "<<five<<std::endl;
    if (five < 5) {
        five = 5;
        xdbg<<"'five' => "<<five<<std::endl; 
    }
    if (3*five-1 >= int(bright_list.size())) {
        std::ostringstream err;
        err<<"Too few objects in bright_list.  Increase startn1.";
        throw StarFinderException(err.str());
    }
    int peak_start = 0;
    double peakWidth = bright_list[five-1]->getSize()-bright_list[0]->getSize();
    for(int k=1;k<=2*five;++k) {
        Assert(k+five-1<int(bright_list.size()));
        double sizeRange = bright_list[k+five-1]->getSize() -
            bright_list[k]->getSize();
        if (sizeRange < peakWidth) {
            peakWidth = bright_list[k+five-1]->getSize() -
                bright_list[k]->getSize();
            peak_start = k;
        }
    }
    xdbg<<"found peak starting at "<<peak_start<<
        " with width "<<peakWidth<<std::endl;

    int k1=peak_start, k2=peak_start+five;
    // Now expand list to be 2*bin_size1 wide
    while (
        k1 > 0 && bright_list[peak_start+2]->getSize() -
        bright_list[k1-1]->getSize() < _bin_size1) {
        --k1;
    }
    const int nbright = bright_list.size();
    while (
        k2 < nbright && 
        bright_list[k2]->getSize() - bright_list[k1]->getSize() < 2.*_bin_size1) {
        ++k2;
    }
    xdbg<<"expanded to "<<k1<<','<<k2<<std::endl;

    // Call these objects in the expanded peak stars
    std::vector<PotentialStar*> star_list(
        bright_list.begin()+k1, bright_list.begin()+k2);

    // Fit f using a linear fit with 2 sigma clipping
    fitStellarSizes(f,0,2.,star_list,out_sigma);

    // if sigma is too big, drop either the first or last "star" and try again.
    while(*out_sigma > _max_rms && star_list.size() > 4) {
        double diff1 = star_list.front()->getSize() -
            (*f)(star_list.front()->getPos());
        double diff2 = star_list.back()->getSize() -
            (*f)(star_list.back()->getPos());

        if (std::abs(diff1) > std::abs(diff2)) {
            star_list.erase(star_list.begin());
            xdbg<<"Erased first star.  New size = "<<star_list.size()<<std::endl;
        } else {
            star_list.erase(star_list.end()-1);
            xdbg<<"Erased last star.  New size = "<<star_list.size()<<std::endl;
        }
        fitStellarSizes(f,0,2.,star_list,out_sigma);
    }
}

