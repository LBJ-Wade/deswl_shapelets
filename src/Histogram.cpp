#include <algorithm>
#include "dbg.h"
#include "Form.h"
#include "Histogram.h"

template <typename T>
Histogram<T>::Histogram(double bin_size, double min_value, double max_value) :
    _bin_size(bin_size)
{
    // Add an extra bin to prevent rounding errors from giving you a
    // value outsize the allowed range.
    int nBins = int(ceil((max_value-min_value)/_bin_size))+1;
    _min_value = min_value - bin_size/2.;
    _max_value = max_value + bin_size/2.;
    _refs = std::vector<std::vector<T> >(nBins);
    _values = std::vector<std::vector<double> >(nBins);
    if (value(nBins-1) > _max_value) _max_value = value(nBins-1)+bin_size/4.;

    dbg<<"made histogram:\n";
    dbg<<"minvalue="<<_min_value<<
        " index(minvalue)="<<index(_min_value)<<std::endl;
    dbg<<"maxvalue="<<_max_value<<
        " index(maxvalue)="<<index(_max_value)<<std::endl;
    dbg<<"mini=0  value(mini)="<<value(0)<<std::endl;
    dbg<<"maxi="<<_refs.size()-1;
    dbg<<" value(maxi)="<<value(_refs.size()-1)<<std::endl;
}

template <typename T>
void Histogram<T>::add(double value,const T& ref)
{
    _refs[index(value)].push_back(ref);
    _values[index(value)].push_back(value);
}

template <typename T>
int Histogram<T>::getCount(int i) const
{
    Assert(i<int(_refs.size()));
    Assert(_values[i].size() == _refs[i].size());
    return _refs[i].size();
}

template <typename T>
double Histogram<T>::findPeak(double val1, double val2) const
{
    if (val1<_min_value) val1 = _min_value;
    if (val2>_max_value) val2 = _max_value;
    int i1=index(val1), i2=index(val2);
    Assert(i1<int(_refs.size()));
    Assert(i2<int(_refs.size()));
    int ipeak=i1;
    int maxcount = getCount(i1);

    for(int i=i1+1;i<=i2;++i) {
        if (getCount(i) > maxcount) {
            maxcount = getCount(i);
            ipeak = i;
        }
    }
    return value(ipeak);
}

template <typename T>
bool Histogram<T>::hasSinglePeak(double val1, double val2) const
// Finds the first peak from the left... set to ipeak1
// and the first peak from the right... set to ipeak2
// Returns whether they are the same peak
{
    if (val1<_min_value) val1 = _min_value;
    if (val2>_max_value) val2 = _max_value;
    int i1=index(val1), i2=index(val2);
    Assert(i1<int(_refs.size()));
    Assert(i2<int(_refs.size()));
    int ipeak1=i1,ipeak2=i2;
    while(ipeak1+1<=i2 && getCount(ipeak1+1)>=getCount(ipeak1)) ++ipeak1;
    while(ipeak2-1>=i1 && getCount(ipeak2-1)>=getCount(ipeak2)) --ipeak2;
    return ipeak1 >= ipeak2; // Usually ==, but if plateau, they swap
}

template <typename T>
double Histogram<T>::findValley(double val1, double val2) const
{
    if (val1<_min_value) val1 = _min_value;
    if (val2>_max_value) val2 = _max_value;
    int i1=index(val1), i2=index(val2);
    Assert(i1<int(_refs.size()));
    Assert(i2<int(_refs.size()));
    int ivalley=i1;
    int minCount = getCount(i1);

    for(int i=i1+1;i<=i2;++i) {
        if (getCount(i) < minCount) {
            minCount = getCount(i);
            ivalley = i;
        }
    }
    // Make sure you're not still sloping down.  If so, keep going.
    if (ivalley==i2) {
        int nrefs = _refs.size();
        while(ivalley+1 < nrefs && getCount(ivalley+1) < getCount(ivalley)) 
            ++ivalley;
    }

    return value(ivalley);
}

template <typename T>
double Histogram<T>::findFirstValleyAfter(double val1, bool poisson_noise) const
{
    dbg<<"Start FindFirstValleyAfter "<<val1<<std::endl;
    if (val1<_min_value) val1 = _min_value;
    int i1=index(val1);
    Assert(i1<int(_refs.size()));
    // Keep going until it stops dropping
    // The sqrt bit allows for Poisson noise in the counts
    // Otherwise you can get caught early by ledges.
    // eg. 15 9 10 4 1 3 3 ... would stop at the 9 instead of at the 1
    int ivalley=i1;
    double minCount = getCount(i1);
    int icheck=i1+1;
    const int nrefs = _refs.size();
    while(icheck < nrefs && getCount(icheck) <= minCount +
          (poisson_noise ? int(sqrt(float(getCount(icheck)))) : 0)) {
        if (getCount(icheck) < minCount) {
            minCount = getCount(icheck);
            ivalley = icheck;
        }
        if (++icheck == nrefs) break;
    }
    dbg<<"valley = "<<value(ivalley)<<std::endl;
    return value(ivalley);
}

template <typename T>
double Histogram<T>::findFirstValleyBefore(double val1, bool poisson_noise) const
{
    if (val1>_max_value) val1 = _max_value;
    int i1=index(val1);
    Assert(i1<int(_refs.size()));
    // Keep going until it stops dropping
    // The sqrt bit allows for Poisson noise in the counts
    // Otherwise you can get caught early by ledges.
    // eg. 15 9 10 4 1 3 3 ... would stop at the 9 instead of at the 1
    int ivalley=i1;
    double minCount = getCount(i1);
    if (i1 > 0) {
        int icheck=i1-1;
        while(getCount(icheck) <= minCount +
              (poisson_noise ? int(sqrt(float(getCount(icheck)))) : 0)) {
            if (getCount(icheck) < minCount) {
                minCount = getCount(icheck);
                ivalley = icheck;
            }
            if (icheck == 0) break;
            --icheck;
        }
    }
    return value(ivalley);
}

template <typename T>
double Histogram<T>::findFirstValueAfter(double start) const
{
    int i = (start<_min_value) ? index(_min_value) : index(start)+1;
    const int nvalues = _values.size();
    while (i<nvalues && _values[i].size() == 0) ++i;
    if (i == nvalues) return _max_value;
    else return *std::min_element(_values[i].begin(),_values[i].end());
}

template <typename T>
double Histogram<T>::findFirstPeakAfter(double val1, bool poisson_noise) const
{
    if (val1<_min_value) val1 = _min_value;
    int i1=index(val1);
    const int nrefs = _refs.size();
    Assert(i1 < nrefs);
    // Keep going until it stops rising
    int ipeak=i1;
    double maxcount = getCount(i1);
    int icheck=i1+1;
    while(icheck < nrefs && maxcount <= getCount(icheck) +
          (poisson_noise ? int(sqrt(float(getCount(icheck)))) : 0)) {
        if (getCount(icheck) > maxcount) {
            maxcount = getCount(icheck);
            ipeak = icheck;
        }
        if (++icheck == nrefs) break;
    }
    return value(ipeak);
}

template <typename T>
double Histogram<T>::findFirstPeakBefore(double val1, bool poisson_noise) const
{
    if (val1>_max_value) val1 = _max_value;
    int i1=index(val1);
    Assert(i1 < int(_refs.size()));
    // Keep going until it stops rising
    int ipeak=i1;
    double maxcount = getCount(i1);
    if (i1 > 0) {
        int icheck=i1-1;
        while(maxcount <= getCount(icheck) +
              (poisson_noise ? int(sqrt(float(getCount(icheck)))) : 0)) {
            if (getCount(icheck) > maxcount) {
                maxcount = getCount(icheck);
                ipeak = icheck;
            }
            if (icheck == 0) break;
            --icheck;
        }
    }
    return value(ipeak);
}

template <typename T>
int Histogram<T>::getTotalCountBefore(double val1) const
{
    int count=0;
    int i1 = index(val1);
    Assert(i1 < int(_refs.size()));
    for(int i=0;i<i1;++i) count += getCount(i);
    return count;
}

template <typename T>
int Histogram<T>::getTotalCountAfter(double val1) const
{
    int count=0;
    int i1 = index(val1);
    const int nrefs = _refs.size();
    Assert(i1 < nrefs);
    for(int i=i1+1;i<nrefs;++i) count += getCount(i);
    return count;
}

template <typename T>
int Histogram<T>::getRefinedPeakCount(double* peak) const
// Takes the three bins around the peak and find the maximum number of
// objects which fit into a bin width.  
// It does this by sliding a window of width bin_size along a sorted list
// of values and counting how many objects fit into the window. 
{
    std::vector<double> vals = 
        getValuesInRange(*peak-1.5*_bin_size,
                         *peak+1.5*_bin_size);
    Assert(vals.size()>=1);
    std::sort(vals.begin(),vals.end());
    int best_count=0;
    double best_val=*peak;
    // i is the object at the left side of the virtual bin
    // j is the first object after the end of the virtual bin
    const int nvals = vals.size();
    for(int i=0,j=0;j<nvals;++i) {
        // Determine j for this value of i
        while (j<nvals && vals[j] < vals[i]+_bin_size) ++j;
        // If the corresponding count is better than best_count, update it
        Assert(j>i);
        if (j-i > best_count) {
            best_count = j-i;
            best_val = (vals[i]+vals[j-1])/2.;
        }
    }
    *peak = best_val;
    return best_count;
}

template <typename T>
int Histogram<T>::getRefinedValleyCount(double* valley) const
// Just like getRefinedPeakCount, but for a valley.
// One slight difference is that for the best peak, you want to 
// include as many objects as possible, so the "i" value is included
// For the best valley, we imagine that the "i" value is just slightly
// to the left of the start of the bin.  So the count is j-i-1 rather
// than j-i.
{
    std::vector<double> vals = 
        getValuesInRange(*valley-1.5*_bin_size,
                         *valley+1.5*_bin_size);
    const int nvals = vals.size();
    if (nvals == 0) return 0;
    std::sort(vals.begin(),vals.end());
    // start high since want to find the lowest count
    int best_count=nvals;
    double best_val=*valley;
    // i is the last object before the start of the virtual bin
    // j is the first object after the end of the virtual bin
    for(int i=0,j=0; i<nvals && vals[i]<*valley+0.5*_bin_size; ++i) {
        // Determine j for this value of i
        while (j<nvals && vals[j] < vals[i]+_bin_size) ++j;
        // If we can get rid of any i values, given this j, do so.
        double next_bin_start = j<nvals ? vals[j] : *valley+1.5*_bin_size;
        while (i+1<nvals && vals[i+1]+_bin_size < next_bin_start) ++i;
        // If the corresponding count is better than best_count, update it
        Assert(j>i);
        if (j-i-1 < best_count) {
            best_count = j-i-1;
            best_val = vals[i]+_bin_size/2.;
        }
    }
    *valley = best_val;
    return best_count;
}

template <typename T>
int Histogram<T>::operator[](double value) const
{
    if (index(value) == int(_refs.size())) return 0;
    Assert(index(value) < int(_refs.size()));
    return getCount(index(value));
}

template <typename T>
double Histogram<T>::findThresh(double min_val, double max_val) const
// This function finds the threshold point which maximizes the 
// "Normalized Between Group Variance" (var).
// For a given threshold point, we have two groups, the left side and 
// the right side. 
// Define fL and fR to be the fraction of points in each group.
// xL and xR are the centroids of each group (the mean value)
// vL and vR are the variances for each group.
// Then the Between Group Variance is just:
// bgv = fL fR (xR - xL)^2
// Maximizing this works well for distributions which may have different
// numbers of elements so long as the widths of the two distributions
// are similar.
// When the widths are not similar, it is better to normalize this by the
// effective sigmas, or sqrt(variance).
// var = fL fR (xR-xL)^2 / sqrt(vL vR)
// Maximizing this seems to work well for separating stars from galaxies.
//
// To calculate this efficiently, we want to be able to keep track of all
// the values as we go along the histogram to make it order N, rather than
// N^2.
//
// ntot = Sumtot(N(i))
// meantot = Sumtot(N(i)*i)/Ntot
// meanSqntt = Sumtot(N(i)*i*i)/Ntot
//
// (All sums below are from i=0..i1, "sumtot" above means from i=0..maxi)
//
// fL = Sum(N(i))/ntot
// fR = 1-fL
// xL = Sum(N(i)*i)/ntot/fL
// xR = (meantot - fL*xL)/fR
// vL = Sum(N(i)*i*i)/ntot/fL - xL^2
// vR = (meansqtot - fL*(vL+xL^2))/fR
{
    const double EPSILON = 1.e-6;

    dbg<<"starting FindThresh\n";
    if (dbgout) print(*dbgout);

    double ntot=0.,meantot=0.,meansqtot=0.;

    if (min_val < _min_value) min_val = _min_value;
    if (max_val > _max_value) max_val = _max_value;
    int i1 = index(min_val);
    int i2 = index(max_val);
    Assert(i1 < int(_refs.size()));
    Assert(i2 < int(_refs.size()));

    // First calculate the "tot" values:
    for(int i=i1;i<=i2;++i) {
        double dcount = getCount(i);
        ntot += dcount;
        meantot += dcount*i;
        meansqtot += dcount*i*i;
    }
    meantot /= ntot;
    meansqtot /= ntot;
    dbg<<"ntot="<<ntot<<", meantot="<<meantot<<
        ", meansqtot="<<meansqtot<<std::endl;

    double sumN=0.,sumNi=0.,sumNii=0.,bestvar=-1.;
    int besti=-1;

    const double min_variance = 0.25;

    for(int i=i1;i<=i2;++i) {
        double dcount = getCount(i);
        sumN += dcount;
        sumNi += dcount*i;
        sumNii += dcount*i*i;
        double fL = sumN/ntot;
        double fR = 1.-fL;
        Assert(fL >=0. && fR >=0.);
        if (fL == 0. || fR == 0.) continue; // var = 0
        double xL = sumNi/ntot;  // actuall fL*xL now
        double xR = (meantot-xL)/fR;
        xL /= fL;
        double vL = sumNii/ntot; // actually fL*(vL+xL^2)
        double vR = (meansqtot-vL)/fR - xR*xR;
        vL = vL/fL - xL*xL;
        if (vL<min_variance) vL = min_variance;
        if (vR<min_variance) vR = min_variance;
        double var = fL*fR*std::pow(xR-xL,2)/sqrt(vL*vR);
        Assert(var >= 0.);
        if (var > bestvar+EPSILON) {
            besti = i;
            bestvar = var;
        }
    }

    dbg<<"besti = "<<besti<<", bestvar = "<<bestvar<<std::endl;
    Assert(besti >= 0);
    Assert(bestvar > 0.);
    // +1 below forces thresh to be at least 1 bin away from first peak.
    // otherwise if peak is only 1 bin, can get besti = that bin.
    dbg<<"returning thresh = "<<value(besti)<<std::endl;
    return value(besti);
}

template <typename T>
std::vector<T> Histogram<T>::getRefsInRange(double min, double max) const
{
    if (min < _min_value) min = _min_value;
    if (max > _max_value) max = _max_value;
    int i1=index(min);
    int i2=index(max);

    dbg<<"in getrefs: min,max = "<<min<<','<<max<<std::endl;

    std::vector<T> temp;
    const int nrefs1 = _refs[i1].size();
    for(int k=0; k<nrefs1; ++k) {
        if(_values[i1][k]>=min) {
            dbg<<"i1 - add ref @ "<<_values[i1][k]<<std::endl;
            temp.push_back(_refs[i1][k]);
        }
    }
    for(int i=i1+1;i<i2;++i) {
        if (_refs[i].size() > 0) 
            dbg<<"i - add all ("<<_refs[i].size()<<") refs near "<<
                _values[i].front()<<std::endl;
        temp.insert(temp.end(),_refs[i].begin(),_refs[i].end());
    }
    const int nrefs2 = _refs[i2].size();
    for(int k=0; k<nrefs2; ++k) {
        if(_values[i2][k] <= max) {
            dbg<<"i2 - add ref @ "<<_values[i2][k]<<std::endl;
            temp.push_back(_refs[i2][k]);
        }
    }
    return temp;
}

template <typename T>
std::vector<double> Histogram<T>::getValuesInRange(
    double min, double max) const
{
    if (min < _min_value) min = _min_value;
    if (max > _max_value) max = _max_value;
    int i1=index(min);
    int i2=index(max);

    std::vector<double> temp;
    const int nvals1 = _values[i1].size();
    for(int k=0;k<nvals1;++k)
        if(_values[i1][k]>=min) temp.push_back(_values[i1][k]);
    for(int i=i1+1;i<i2;++i)
        temp.insert(temp.end(),_values[i].begin(),_values[i].end());
    const int nvals2 = _values[i2].size();
    for(int k=0;k<nvals2;++k)
        if(_values[i2][k]<=max) temp.push_back(_values[i2][k]);
    return temp;
}

template <typename T>
int Histogram<T>::index(double value) const
{
    Assert (value >= _min_value && value <= _max_value);
    return int(floor((value - _min_value)/_bin_size));
}

template <typename T>
double Histogram<T>::value(int i) const
{
    Assert(i<int(_refs.size()));
    return _bin_size*(i+0.5)+_min_value;
}

template <typename T>
void Histogram<T>::print(std::ostream& fout,double val1,double val2) const
{
    if (val1 < _min_value) val1 = _min_value;
    if (val2 > _max_value) val2 = _max_value;
    int i1=index(val1);
    int i2=index(val2);
    Assert(i1 < int(_refs.size()));
    Assert(i2 < int(_refs.size()));
    Form sci2; sci2.sci().width(10).prec(2);

    for(int i=i1;i<=i2;++i) {
        //if ((i-i1)%5 == 0) 
        fout << sci2(value(i));
        //else fout << "           ";
        for(int j=0;j<getCount(i);++j) fout<<'*';
        fout<<std::endl;
    }
}

#include "PotentialStar.h"

template class Histogram<PotentialStar*>;
