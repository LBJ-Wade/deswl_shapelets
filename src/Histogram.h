//---------------------------------------------------------------------------
#ifndef HistogramH
#define HistogramH
//---------------------------------------------------------------------------

#include <vector>
#include <iostream>

template <class T>
class Histogram {

public:

  Histogram(double binsize, double minvalue, double maxvalue);
  ~Histogram() {}

  void Add(double value,const T& ref);
  double FindPeak(double minval, double maxval) const;
  bool SinglePeak(double minval, double maxval) const;
  double FindValley(double minval, double maxval) const;
  double FindFirstValueAfter(double start) const;
  double FindFirstValleyAfter(double val1,bool poissonnoise=false) const;
  double FindFirstValleyBefore(double val1,bool poissonnoise=false) const;
  double FindFirstPeakAfter(double val1,bool poissonnoise=false) const;
  double FindFirstPeakBefore(double val1,bool poissonnoise=false) const;
  size_t GetTotalCountBefore(double val1) const;
  size_t GetTotalCountAfter(double val1) const;
  size_t operator[](double value) const;
  double FindThresh(double minval, double maxval) const;
  std::vector<T> GetRefsInRange(double min, double max) const;
  std::vector<double> GetValuesInRange(double min, double max) const;
  size_t GetRefinedPeakCount(double* peak) const;
  size_t GetRefinedValleyCount(double* valley) const;
  void Print(std::ostream& fout,double val1=-1.e10,double val2=1.e10) const;

private:

  size_t GetCount(size_t i) const;
  size_t index(double value) const;
  double value(size_t i) const;

  double itsbinsize,itsminvalue,itsmaxvalue;
  std::vector<std::vector<T> > itsrefs;
  std::vector<std::vector<double> > itsvalues;
};



#endif
