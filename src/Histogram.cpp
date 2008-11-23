#include "Histogram.h"
#include <algorithm>
#include "Form.h"
#include "dbg.h"

#define EPSILON 1.e-6

template <class T> inline T SQR(const T& x) { return x*x; }

template <class T>
Histogram<T>::Histogram(double binsize, double minvalue, double maxvalue) :
    itsbinsize(binsize)
{
  // Add an extra bin to prevent rounding errors from giving you a
  // value outsize the allowed range.
  size_t nbins = (size_t) ceil((maxvalue-minvalue)/itsbinsize)+1;
  itsminvalue = minvalue - binsize/2.;
  itsmaxvalue = maxvalue + binsize/2.;
  itsrefs = std::vector<std::vector<T> >(nbins);
  itsvalues = std::vector<std::vector<double> >(nbins);
  if (value(nbins-1) > itsmaxvalue) itsmaxvalue = value(nbins-1)+binsize/4.;

  dbg<<"made histogram:\n";
  dbg<<"minvalue="<<itsminvalue<<" index(minvalue)="<<index(itsminvalue)<<std::endl;
  dbg<<"maxvalue="<<itsmaxvalue<<" index(maxvalue)="<<index(itsmaxvalue)<<std::endl;
  dbg<<"mini=0  value(mini)="<<value(0)<<std::endl;
  dbg<<"maxi="<<itsrefs.size()-1;
  dbg<<" value(maxi)="<<value(itsrefs.size()-1)<<std::endl;
}

template <class T>
void Histogram<T>::Add(double value,const T& ref)
{
  itsrefs[index(value)].push_back(ref);
  itsvalues[index(value)].push_back(value);
}

template <class T>
size_t Histogram<T>::GetCount(size_t i) const
{
  Assert(i<itsrefs.size());
  Assert(itsvalues[i].size() == itsrefs[i].size());
  return itsrefs[i].size();
}
  
template <class T>
double Histogram<T>::FindPeak(double val1, double val2) const
{
  if (val1<itsminvalue) val1 = itsminvalue;
  if (val2>itsmaxvalue) val2 = itsmaxvalue;
  size_t i1=index(val1), i2=index(val2);
  Assert(i1<itsrefs.size());
  Assert(i2<itsrefs.size());
  size_t ipeak=i1;
  size_t maxcount = GetCount(i1);

  for(size_t i=i1+1;i<=i2;i++) {
    if (GetCount(i) > maxcount) {
      maxcount = GetCount(i);
      ipeak = i;
    }
  }
  return value(ipeak);
}

template <class T>
bool Histogram<T>::SinglePeak(double val1, double val2) const
// Finds the first peak from the left... set to ipeak1
// and the first peak from the right... set to ipeak2
// Returns whether they are the same peak
{
  if (val1<itsminvalue) val1 = itsminvalue;
  if (val2>itsmaxvalue) val2 = itsmaxvalue;
  size_t i1=index(val1), i2=index(val2);
  Assert(i1<itsrefs.size());
  Assert(i2<itsrefs.size());
  size_t ipeak1=i1,ipeak2=i2;
  while(ipeak1+1<=i2 && GetCount(ipeak1+1)>=GetCount(ipeak1)) 
    ipeak1++;
  while(ipeak2-1>=i1 && GetCount(ipeak2-1)>=GetCount(ipeak2))
    ipeak2--;
  return ipeak1 >= ipeak2; // Usually ==, but if plateau, they swap
}

template <class T>
double Histogram<T>::FindValley(double val1, double val2) const
{
  if (val1<itsminvalue) val1 = itsminvalue;
  if (val2>itsmaxvalue) val2 = itsmaxvalue;
  size_t i1=index(val1), i2=index(val2);
  Assert(i1<itsrefs.size());
  Assert(i2<itsrefs.size());
  size_t ivalley=i1;
  size_t mincount = GetCount(i1);

  for(size_t i=i1+1;i<=i2;i++) {
    if (GetCount(i) < mincount) {
      mincount = GetCount(i);
      ivalley = i;
    }
  }
  // Make sure you're not still sloping down.  If so, keep going.
  if (ivalley==i2) 
    while(ivalley+1<itsrefs.size()&&GetCount(ivalley+1)<GetCount(ivalley)) 
      ivalley++;
   
  return value(ivalley);
}

template <class T>
double Histogram<T>::FindFirstValleyAfter(double val1,bool poissonnoise) const
{
  dbg<<"Start FindFirstValleyAfter "<<val1<<std::endl;
  if (val1<itsminvalue) val1 = itsminvalue;
  size_t i1=index(val1);
  Assert(i1<itsrefs.size());
  // Keep going until it stops dropping
  // The sqrt bit allows for Poisson noise in the counts
  // Otherwise you can get caught early by ledges.
  // eg. 15 9 10 4 1 3 3 ... would stop at the 9 instead of at the 1
  size_t ivalley=i1;
  double mincount = GetCount(i1);
  size_t icheck=i1+1;
  while(icheck < itsrefs.size() && GetCount(icheck) <= mincount +
      (poissonnoise ? (size_t)sqrt((float)GetCount(icheck)) : 0)) {
    if (GetCount(icheck) < mincount) {
      mincount = GetCount(icheck);
      ivalley = icheck;
    }
    if (++icheck == itsrefs.size()) break;
  }
  dbg<<"valley = "<<value(ivalley)<<std::endl;
  return value(ivalley);
}

template <class T>
double Histogram<T>::FindFirstValleyBefore(double val1,bool poissonnoise) const
{
  if (val1>itsmaxvalue) val1 = itsmaxvalue;
  size_t i1=index(val1);
  Assert(i1<itsrefs.size());
  // Keep going until it stops dropping
  // The sqrt bit allows for Poisson noise in the counts
  // Otherwise you can get caught early by ledges.
  // eg. 15 9 10 4 1 3 3 ... would stop at the 9 instead of at the 1
  size_t ivalley=i1;
  double mincount = GetCount(i1);
  if (i1 > 0) {
    size_t icheck=i1-1;
    while(GetCount(icheck) <= mincount +
	(poissonnoise ? (size_t)sqrt((float)GetCount(icheck)) : 0)) {
      if (GetCount(icheck) < mincount) {
	mincount = GetCount(icheck);
	ivalley = icheck;
      }
      if (icheck == 0) break;
      --icheck;
    }
  }
  return value(ivalley);
}

template <class T>
double Histogram<T>::FindFirstValueAfter(double start) const
{
  size_t i = (start<itsminvalue) ? index(itsminvalue) : index(start)+1;
  while (i<itsvalues.size() && itsvalues[i].size() == 0) i++;
  if (i == itsvalues.size()) return itsmaxvalue;
  else return *std::min_element(itsvalues[i].begin(),itsvalues[i].end());
}
  
template <class T>
double Histogram<T>::FindFirstPeakAfter(double val1,bool poissonnoise) const
{
  if (val1<itsminvalue) val1 = itsminvalue;
  size_t i1=index(val1);
  Assert(i1 < itsrefs.size());
  // Keep going until it stops rising
  size_t ipeak=i1;
  double maxcount = GetCount(i1);
  size_t icheck=i1+1;
  while(icheck < itsrefs.size() && maxcount <= GetCount(icheck) +
      (poissonnoise ? (size_t)sqrt((float)GetCount(icheck)) : 0)) {
    if (GetCount(icheck) > maxcount) {
      maxcount = GetCount(icheck);
      ipeak = icheck;
    }
    if (++icheck == itsrefs.size()) break;
  }
  return value(ipeak);
}

template <class T>
double Histogram<T>::FindFirstPeakBefore(double val1,bool poissonnoise) const
{
  if (val1>itsmaxvalue) val1 = itsmaxvalue;
  size_t i1=index(val1);
  Assert(i1 < itsrefs.size());
  // Keep going until it stops rising
  size_t ipeak=i1;
  double maxcount = GetCount(i1);
  if (i1 > 0) {
    size_t icheck=i1-1;
    while(maxcount <= GetCount(icheck) +
	(poissonnoise ? (size_t)sqrt((float)GetCount(icheck)) : 0)) {
      if (GetCount(icheck) > maxcount) {
	maxcount = GetCount(icheck);
	ipeak = icheck;
      }
      if (icheck == 0) break;
      --icheck;
    }
  }
  return value(ipeak);
}

template <class T>
size_t Histogram<T>::GetTotalCountBefore(double val1) const
{
  size_t count=0;
  size_t i1 = index(val1);
  Assert(i1 < itsrefs.size());
  for(size_t i=0;i<i1;i++) count += GetCount(i);
  return count;
}

template <class T>
size_t Histogram<T>::GetTotalCountAfter(double val1) const
{
  size_t count=0;
  size_t i1 = index(val1);
  Assert(i1 < itsrefs.size());
  for(size_t i=i1+1;i<itsrefs.size();i++) count += GetCount(i);
  return count;
}

template <class T>
size_t Histogram<T>::GetRefinedPeakCount(double* peak) const
// Takes the three bins around the peak and find the maximum number of
// objects which fit into a bin width.  
// It does this by sliding a window of width binsize along a sorted list
// of values and counting how many objects fit into the window. 
{
  std::vector<double> vals = GetValuesInRange(*peak-1.5*itsbinsize,
    *peak+1.5*itsbinsize);
  Assert(vals.size()>=1);
  std::sort(vals.begin(),vals.end());
  size_t bestcount=0;
  double bestval=*peak;
  // i is the object at the left side of the virtual bin
  // j is the first object after the end of the virtual bin
  for(size_t i=0,j=0;j<vals.size();i++) {
    // Determine j for this value of i
    while (j<vals.size() && vals[j] < vals[i]+itsbinsize) j++;
    // If the corresponding count is better than bestcount, update it
    Assert(j>i);
    if (j-i > bestcount) {
      bestcount = j-i;
      bestval = (vals[i]+vals[j-1])/2.;
      //dbg<<"bestcount = "<<bestcount<<", bestval = "<<bestval<<std::endl;
      //dbg<<"i,j = "<<i<<','<<j<<std::endl;
    }
  }
  *peak = bestval;
  return bestcount;
}

template <class T>
size_t Histogram<T>::GetRefinedValleyCount(double* valley) const
// Just like GetRefinedPeakCount, but for a valley.
// One slight difference is that for the best peak, you want to 
// include as many objects as possible, so the "i" value is included
// For the best valley, we imagine that the "i" value is just slightly
// to the left of the start of the bin.  So the count is j-i-1 rather
// than j-i.
{
  std::vector<double> vals = GetValuesInRange(*valley-1.5*itsbinsize,
    *valley+1.5*itsbinsize);
  if (vals.size() == 0) return 0;
  std::sort(vals.begin(),vals.end());
  // start high since want to find the lowest count
  size_t bestcount=vals.size();
  double bestval=*valley;
  // i is the last object before the start of the virtual bin
  // j is the first object after the end of the virtual bin
  for(size_t i=0,j=0;i<vals.size()&&vals[i]<*valley+0.5*itsbinsize;i++) {
    // Determine j for this value of i
    while (j<vals.size() && vals[j] < vals[i]+itsbinsize) j++;
    // If we can get rid of any i values, given this j, do so.
    double nextbinstart = j<vals.size() ? vals[j] : *valley+1.5*itsbinsize;
    while (i+1<vals.size() && vals[i+1]+itsbinsize < nextbinstart) i++;
    // If the corresponding count is better than bestcount, update it
    Assert(j>i);
    if (j-i-1 < bestcount) {
      bestcount = j-i-1;
      bestval = vals[i]+itsbinsize/2.;
      //dbg<<"bestcount = "<<bestcount<<", bestval = "<<bestval<<std::endl;
      //dbg<<"i,j = "<<i<<','<<j<<std::endl;
    }
  }
  *valley = bestval;
  return bestcount;
}

template <class T>
size_t Histogram<T>::operator[](double value) const
{
  if (index(value) == itsrefs.size()) return 0;
  Assert(index(value) < itsrefs.size());
  return GetCount(index(value));
}

template <class T>
double Histogram<T>::FindThresh(double minval, double maxval) const
// This function finds the threshold point which maximizes the 
// "Normalized Between Group Variance" (nbgv).
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
// nbgv = fL fR (xR-xL)^2 / sqrt(vL vR)
// Maximizing this seems to work well for separating stars from galaxies.
//
// To calculate this efficiently, we want to be able to keep track of all
// the values as we go along the histogram to make it order N, rather than
// N^2.
//
// ntot = Sumtot(N(i))
// meantot = Sumtot(N(i)*i)/Ntot
// meansqtot = Sumtot(N(i)*i*i)/Ntot
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
  dbg<<"starting FindThresh\n";
  if (dbgout) Print(*dbgout);

  double ntot=0.,meantot=0.,meansqtot=0.;

  if (minval < itsminvalue) minval = itsminvalue;
  if (maxval > itsmaxvalue) maxval = itsmaxvalue;
  size_t i1 = index(minval);
  size_t i2 = index(maxval);
  Assert(i1 < itsrefs.size());
  Assert(i2 < itsrefs.size());

  // First calculate the "tot" values:
  for(size_t i=i1;i<=i2;i++) {
    double dcount = GetCount(i);
    ntot += dcount;
    meantot += dcount*i;
    meansqtot += dcount*i*i;
  }
  meantot /= ntot;
  meansqtot /= ntot;
  dbg<<"ntot="<<ntot<<", meantot="<<meantot<<", meansqtot="<<meansqtot<<std::endl;

  double sumn=0.,sumni=0.,sumnii=0.,bestnbgv=-1.;
  int besti=-1;

  const double minvariance = 0.25;

  for(size_t i=i1;i<=i2;i++) {
    double dcount = GetCount(i);
    sumn += dcount;
    sumni += dcount*i;
    sumnii += dcount*i*i;
    double fL = sumn/ntot;
    double fR = 1.-fL;
    Assert(fL >=0. && fR >=0.);
    if (fL == 0. || fR == 0.) continue; // nbgv = 0
    double xL = sumni/ntot;  // actuall fL*xL now
    double xR = (meantot-xL)/fR;
    xL /= fL;
    double vL = sumnii/ntot; // actually fL*(vL+xL^2)
    double vR = (meansqtot-vL)/fR - xR*xR;
    vL = vL/fL - xL*xL;
    if (vL<minvariance) vL = minvariance;
    if (vR<minvariance) vR = minvariance;
    double nbgv = fL*fR*SQR(xR-xL)/sqrt(vL*vR);
    Assert(nbgv >= 0.);
    if (nbgv > bestnbgv+EPSILON) {
      besti = i;
      bestnbgv = nbgv;
    }
  }

  dbg<<"besti = "<<besti<<", bestnbgv = "<<bestnbgv<<std::endl;
  Assert(besti >= 0);
  Assert(bestnbgv > 0.);
  // +1 below forces thresh to be at least 1 bin away from first peak.
  // otherwise if peak is only 1 bin, can get besti = that bin.
  dbg<<"returning thresh = "<<value(besti)<<std::endl;
  return value(besti);
}
 
template <class T>
std::vector<T> Histogram<T>::GetRefsInRange(double min, double max) const
{
  if (min < itsminvalue) min = itsminvalue;
  if (max > itsmaxvalue) max = itsmaxvalue;
  size_t i1=index(min);
  size_t i2=index(max);

  dbg<<"in getrefs: min,max = "<<min<<','<<max<<std::endl;

  std::vector<T> temp;
  for(size_t k=0;k<itsrefs[i1].size();k++) {
    if(itsvalues[i1][k]>=min) {
      dbg<<"i1 - add ref @ "<<itsvalues[i1][k]<<std::endl;
      temp.push_back(itsrefs[i1][k]);
    }
  }
  for(size_t i=i1+1;i<i2;i++) {
    if (itsrefs[i].size() > 0) 
      dbg<<"i - add all ("<<itsrefs[i].size()<<") refs near "<<itsvalues[i].front()<<std::endl;
    temp.insert(temp.end(),itsrefs[i].begin(),itsrefs[i].end());
  }
  for(size_t k=0;k<itsrefs[i2].size();k++) {
    if(itsvalues[i2][k]<=max) {
      dbg<<"i2 - add ref @ "<<itsvalues[i2][k]<<std::endl;
      temp.push_back(itsrefs[i2][k]);
    }
  }
  return temp;
}

template <class T>
std::vector<double> Histogram<T>::GetValuesInRange(double min, double max) const
{
  if (min < itsminvalue) min = itsminvalue;
  if (max > itsmaxvalue) max = itsmaxvalue;
  size_t i1=index(min);
  size_t i2=index(max);

  std::vector<double> temp;
  for(size_t k=0;k<itsvalues[i1].size();k++)
    if(itsvalues[i1][k]>=min) temp.push_back(itsvalues[i1][k]);
  for(size_t i=i1+1;i<i2;i++)
    temp.insert(temp.end(),itsvalues[i].begin(),itsvalues[i].end());
  for(size_t k=0;k<itsvalues[i2].size();k++)
    if(itsvalues[i2][k]<=max) temp.push_back(itsvalues[i2][k]);
  return temp;
}

template <class T>
size_t Histogram<T>::index(double value) const
{
  Assert (value >= itsminvalue && value <= itsmaxvalue);
  return (size_t) floor((value - itsminvalue)/itsbinsize);
}

template <class T>
double Histogram<T>::value(size_t i) const
{
  Assert(i<itsrefs.size());
  return itsbinsize*(i+0.5)+itsminvalue;
}

template <class T>
void Histogram<T>::Print(std::ostream& fout,double val1,double val2) const
{
  if (val1 < itsminvalue) val1 = itsminvalue;
  if (val2 > itsmaxvalue) val2 = itsmaxvalue;
  size_t i1=index(val1);
  size_t i2=index(val2);
  Assert(i1 < itsrefs.size());
  Assert(i2 < itsrefs.size());
  Form sci2; sci2.sci().width(10).prec(2);

  for(size_t i=i1;i<=i2;i++) {
    //if ((i-i1)%5 == 0) 
    fout << sci2(value(i));
    //else fout << "           ";
    for(size_t j=0;j<GetCount(i);j++) fout<<'*';
    fout<<std::endl;
  }
}

#include "PotentialStar.h"

template class Histogram<PotentialStar*>;
