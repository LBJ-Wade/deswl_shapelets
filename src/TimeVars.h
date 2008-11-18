#ifndef TIME_VARS_H
#define TIME_VARS_H

#include <sys/time.h>

struct OverallFitTimes {
  OverallFitTimes() :
    ts_native_integ(0.), ts_native_centroid(0.),
    ts_native_gamma(0.), ts_native_mu(0.),
    ts_native_fixflux(0.), ts_native_final(0.),
    tf_native_integ(0.), tf_native_centroid(0.),
    tf_native_gamma(0.), tf_native_mu(0.),
    tf_native_fixflux(0.), tf_native_final(0.),
    ts_mu_integ(0.), ts_mu_centroid(0.),
    ts_mu_gamma(0.), ts_mu_mu(0.),
    ts_mu_fixflux(0.), ts_mu_final(0.),
    tf_mu_integ(0.), tf_mu_centroid(0.),
    tf_mu_gamma(0.), tf_mu_mu(0.),
    tf_mu_fixflux(0.), tf_mu_final(0.),
    //ts_full_integ(0.), ts_full_centroid(0.),
    //ts_full_gamma(0.), ts_full_mu(0.),
    //ts_full_fixflux(0.), ts_full_final(0.),
    //tf_full_integ(0.), tf_full_centroid(0.),
    //tf_full_gamma(0.), tf_full_mu(0.),
    //tf_full_fixflux(0.), tf_full_final(0.),
    ts_gamma_integ(0.), ts_gamma_centroid(0.),
    ts_gamma_gamma(0.), ts_gamma_mu(0.),
    ts_gamma_fixflux(0.), ts_gamma_final(0.),
    tf_gamma_integ(0.), tf_gamma_centroid(0.),
    tf_gamma_gamma(0.), tf_gamma_mu(0.),
    tf_gamma_fixflux(0.), tf_gamma_final(0.),
    nf_range1(0), nf_range2(0),
    nf_pixflag1(0), nf_npix1(0), nf_small(0),
    nf_pixflag2(0), nf_npix2(0),
    ns_native(0), nf_native(0),
    ns_mu(0), nf_mu(0),
    //ns_full(0), nf_full(0),
    ns_gamma(0), nf_gamma(0) {}

  OverallFitTimes& operator+=(const OverallFitTimes& rhs)
  {
    ts_native_integ += rhs.ts_native_integ;
    ts_native_centroid += rhs.ts_native_centroid;
    ts_native_gamma += rhs.ts_native_gamma;
    ts_native_mu += rhs.ts_native_mu;
    ts_native_fixflux += rhs.ts_native_fixflux;
    ts_native_final += rhs.ts_native_final;
    tf_native_integ += rhs.tf_native_integ;
    tf_native_centroid += rhs.tf_native_centroid;
    tf_native_gamma += rhs.tf_native_gamma;
    tf_native_mu += rhs.tf_native_mu;
    tf_native_fixflux += rhs.tf_native_fixflux;
    tf_native_final += rhs.tf_native_final;
    ts_mu_integ += rhs.ts_mu_integ;
    ts_mu_centroid += rhs.ts_mu_centroid;
    ts_mu_gamma += rhs.ts_mu_gamma;
    ts_mu_mu += rhs.ts_mu_mu;
    ts_mu_fixflux += rhs.ts_mu_fixflux;
    ts_mu_final += rhs.ts_mu_final;
    tf_mu_integ += rhs.tf_mu_integ;
    tf_mu_centroid += rhs.tf_mu_centroid;
    tf_mu_gamma += rhs.tf_mu_gamma;
    tf_mu_mu += rhs.tf_mu_mu;
    tf_mu_fixflux += rhs.tf_mu_fixflux;
    tf_mu_final += rhs.tf_mu_final;
    //ts_full_integ += rhs.ts_full_integ;
    //ts_full_centroid += rhs.ts_full_centroid;
    //ts_full_gamma += rhs.ts_full_gamma;
    //ts_full_mu += rhs.ts_full_mu;
    //ts_full_fixflux += rhs.ts_full_fixflux;
    //ts_full_final += rhs.ts_full_final;
    //tf_full_integ += rhs.tf_full_integ;
    //tf_full_centroid += rhs.tf_full_centroid;
    //tf_full_gamma += rhs.tf_full_gamma;
    //tf_full_mu += rhs.tf_full_mu;
    //tf_full_fixflux += rhs.tf_full_fixflux;
    //tf_full_final += rhs.tf_full_final;
    ts_gamma_integ += rhs.ts_gamma_integ;
    ts_gamma_centroid += rhs.ts_gamma_centroid;
    ts_gamma_gamma += rhs.ts_gamma_gamma;
    ts_gamma_mu += rhs.ts_gamma_mu;
    ts_gamma_fixflux += rhs.ts_gamma_fixflux;
    ts_gamma_final += rhs.ts_gamma_final;
    tf_gamma_integ += rhs.tf_gamma_integ;
    tf_gamma_centroid += rhs.tf_gamma_centroid;
    tf_gamma_gamma += rhs.tf_gamma_gamma;
    tf_gamma_mu += rhs.tf_gamma_mu;
    tf_gamma_fixflux += rhs.tf_gamma_fixflux;
    tf_gamma_final += rhs.tf_gamma_final;
    nf_range1 += rhs.nf_range1;
    nf_range2 += rhs.nf_range2;
    nf_pixflag1 += rhs.nf_pixflag1;
    nf_npix1 += rhs.nf_npix1;
    nf_small += rhs.nf_small;
    nf_pixflag2 += rhs.nf_pixflag2;
    nf_npix2 += rhs.nf_npix2;
    ns_native += rhs.ns_native;
    nf_native += rhs.nf_native;
    ns_mu += rhs.ns_mu;
    nf_mu += rhs.nf_mu;
    //ns_full += rhs.ns_full;
    //nf_full += rhs.nf_full;
    ns_gamma += rhs.ns_gamma;
    nf_gamma += rhs.nf_gamma;
    return *this;
  }

  void Write(std::ostream& os) const
  {
    os<<"N_Rejected for range error - distortion = "<<nf_range1<<std::endl;
    os<<"N_Rejected for range error - psf interp = "<<nf_range2<<std::endl;
    os<<"N_Rejected for pixel flag (#1) = "<<nf_pixflag1<<std::endl;
    os<<"N_Rejected for too few pixels (#1) = "<<nf_npix1<<std::endl;

    os<<"Native fits:\n";
    os<<"  N_Success = "<<ns_native<<std::endl;
    os<<"  Times = "<<ts_native_integ<<"  "<<ts_native_centroid<<"  "<<
      ts_native_gamma<<"  "<<ts_native_mu<<"  "<<
      ts_native_fixflux<<"  "<<ts_native_final<<std::endl;
    os<<"  N_Fail = "<<nf_native<<std::endl;
    os<<"  Times = "<<tf_native_integ<<"  "<<tf_native_centroid<<"  "<<
      tf_native_gamma<<"  "<<tf_native_mu<<"  "<<
      tf_native_fixflux<<"  "<<tf_native_final<<std::endl;

    os<<"N_Rejected for being too small = "<<nf_small<<std::endl;
    os<<"N_Rejected for pixel flag (#2) = "<<nf_pixflag2<<std::endl;
    os<<"N_Rejected for too few pixels (#2) = "<<nf_npix2<<std::endl;

    os<<"Mu fits:\n";
    os<<"  N_Success = "<<ns_mu<<std::endl;
    os<<"  Times = "<<ts_mu_integ<<"  "<<ts_mu_centroid<<"  "<<
      ts_mu_gamma<<"  "<<ts_mu_mu<<"  "<<
      ts_mu_fixflux<<"  "<<ts_mu_final<<std::endl;
    os<<"  N_Fail = "<<nf_mu<<std::endl;
    os<<"  Times = "<<tf_mu_integ<<"  "<<tf_mu_centroid<<"  "<<
      tf_mu_gamma<<"  "<<tf_mu_mu<<"  "<<
      tf_mu_fixflux<<"  "<<tf_mu_final<<std::endl;

    /*
    os<<"Full fits:\n";
    os<<"  N_Success = "<<ns_full<<std::endl;
    os<<"  Times = "<<ts_full_integ<<"  "<<ts_full_centroid<<"  "<<
      ts_full_gamma<<"  "<<ts_full_mu<<"  "<<
      ts_full_fixflux<<"  "<<ts_full_final<<std::endl;
    os<<"  N_Fail = "<<nf_full<<std::endl;
    os<<"  Times = "<<tf_full_integ<<"  "<<tf_full_centroid<<"  "<<
      tf_full_gamma<<"  "<<tf_full_mu<<"  "<<
      tf_full_fixflux<<"  "<<tf_full_final<<std::endl;
      */

    os<<"Gamma fits:\n";
    os<<"  N_Success = "<<ns_gamma<<std::endl;
    os<<"  Times = "<<ts_gamma_integ<<"  "<<ts_gamma_centroid<<"  "<<
      ts_gamma_gamma<<"  "<<ts_gamma_mu<<"  "<<
      ts_gamma_fixflux<<"  "<<ts_gamma_final<<std::endl;
    os<<"  N_Fail = "<<nf_gamma<<std::endl;
    os<<"  Times = "<<tf_gamma_integ<<"  "<<tf_gamma_centroid<<"  "<<
      tf_gamma_gamma<<"  "<<tf_gamma_mu<<"  "<<
      tf_gamma_fixflux<<"  "<<tf_gamma_final<<std::endl;
  }

  double ts_native_integ;
  double ts_native_centroid;
  double ts_native_gamma;
  double ts_native_mu;
  double ts_native_fixflux;
  double ts_native_final;
  double tf_native_integ;
  double tf_native_centroid;
  double tf_native_gamma;
  double tf_native_mu;
  double tf_native_fixflux;
  double tf_native_final;
  double ts_mu_integ;
  double ts_mu_centroid;
  double ts_mu_gamma;
  double ts_mu_mu;
  double ts_mu_fixflux;
  double ts_mu_final;
  double tf_mu_integ;
  double tf_mu_centroid;
  double tf_mu_gamma;
  double tf_mu_mu;
  double tf_mu_fixflux;
  double tf_mu_final;
  //double ts_full_integ;
  //double ts_full_centroid;
  //double ts_full_gamma;
  //double ts_full_mu;
  //double ts_full_fixflux;
  //double ts_full_final;
  //double tf_full_integ;
  //double tf_full_centroid;
  //double tf_full_gamma;
  //double tf_full_mu;
  //double tf_full_fixflux;
  //double tf_full_final;
  double ts_gamma_integ;
  double ts_gamma_centroid;
  double ts_gamma_gamma;
  double ts_gamma_mu;
  double ts_gamma_fixflux;
  double ts_gamma_final;
  double tf_gamma_integ;
  double tf_gamma_centroid;
  double tf_gamma_gamma;
  double tf_gamma_mu;
  double tf_gamma_fixflux;
  double tf_gamma_final;
  int nf_range1;
  int nf_range2;
  int nf_pixflag1;
  int nf_npix1;
  int nf_small;
  int nf_pixflag2;
  int nf_npix2;
  int ns_native;
  int nf_native;
  int ns_mu;
  int nf_mu;
  //int ns_full;
  //int nf_full;
  int ns_gamma;
  int nf_gamma;
};

inline std::ostream& operator<<(std::ostream& os, const OverallFitTimes& t)
{ t.Write(os); return os; }

#endif
