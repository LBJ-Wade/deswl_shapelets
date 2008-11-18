#ifndef LOG_H
#define LOG_H

#include "Params.h"
#include <fstream>

struct ShearLog {

  ShearLog(std::string logfile) :
    exitcode(SUCCESS), ngals(0),
    nf_range1(0), nf_range2(0),
    nf_edge1(0), nf_npix1(0), nf_small(0),
    nf_edge2(0), nf_npix2(0),
    nf_tmverror(0), nf_othererror(0),
    ns_native(0), nf_native(0),
    ns_mu(0), nf_mu(0),
    ns_gamma(0), nf_gamma(0),
    logout(logfile.c_str(),std::ios_base::app) 
  {
    if (!logout) {
       throw std::runtime_error(
	   std::string("Error: Unable to open logfile ") + logfile
	   + " for append output");
    }
  }

  ~ShearLog() { WriteLog(); logout.close(); }

  void WriteLog() const
  {
    logout << exitcode <<"  "<<ngals<<"  "<<
      ns_gamma <<"  "<<ns_native<<"  "<<
      nf_range1<<"  "<<nf_range2<<"  "<<nf_edge1<<"  "<<nf_npix1<<"  "<<
      nf_native<<"  "<<nf_small<<"  "<<nf_edge2<<"  "<<nf_npix2<<"  "<<
      nf_tmverror<<"  "<<nf_othererror<<"  "<<
      nf_mu<<"  "<<nf_gamma<<std::endl;
  }

  void Write(std::ostream& os) const
  {
    os<<"N_Rejected for range error - distortion = "<<nf_range1<<std::endl;
    os<<"N_Rejected for range error - psf interp = "<<nf_range2<<std::endl;
    os<<"N_Rejected for pixel flag (#1) = "<<nf_edge1<<std::endl;
    os<<"N_Rejected for too few pixels (#1) = "<<nf_npix1<<std::endl;

    os<<"Native fits:\n";
    os<<"  N_Success = "<<ns_native<<std::endl;
    os<<"  N_Fail = "<<nf_native<<std::endl;

    os<<"N_Rejected for being too small = "<<nf_small<<std::endl;
    os<<"N_Rejected for pixel flag (#2) = "<<nf_edge2<<std::endl;
    os<<"N_Rejected for too few pixels (#2) = "<<nf_npix2<<std::endl;

    os<<"Mu fits:\n";
    os<<"  N_Success = "<<ns_mu<<std::endl;
    os<<"  N_Fail = "<<nf_mu<<std::endl;

    os<<"Gamma fits:\n";
    os<<"  N_Success = "<<ns_gamma<<std::endl;
    os<<"  N_Fail = "<<nf_gamma<<std::endl;

    os<<"N_Error: TMV Error caught = "<<nf_tmverror<<std::endl;
    os<<"N_Error: Other caught = "<<nf_othererror<<std::endl;
  }

  ExitCode exitcode;

  int ngals;
  int nf_range1;
  int nf_range2;
  int nf_edge1;
  int nf_npix1;
  int nf_small;
  int nf_edge2;
  int nf_npix2;
  int nf_tmverror;
  int nf_othererror;
  int ns_native;
  int nf_native;
  int ns_mu;
  int nf_mu;
  int ns_gamma;
  int nf_gamma;

  mutable std::ofstream logout;

};

struct PSFLog {

  PSFLog(std::string logfile) :
    exitcode(SUCCESS), nstars(0),
    nf_range(0), nf_edge(0), nf_npix(0),
    nf_tmverror(0), nf_othererror(0),
    ns_psf(0), nf_psf(0),
    logout(logfile.c_str(),std::ios_base::app) 
  {
    if (!logout) {
       throw std::runtime_error(
	   std::string("Error: Unable to open logfile ") + logfile
	   + " for append output");
    }
  }

  ~PSFLog() { WriteLog(); logout.close(); }

  void WriteLog() const
  {
    logout << exitcode <<"  "<<nstars<<"  "<<
      ns_psf <<"  "<<
      nf_range<<"  "<<nf_edge<<"  "<<nf_npix<<"  "<<
      nf_tmverror<<"  "<<nf_othererror<<"  "<<
      nf_psf<<std::endl;
  }

  void Write(std::ostream& os) const
  {
    os<<"N_Rejected for range error - distortion = "<<nf_range<<std::endl;
    os<<"N_Rejected for edge error = "<<nf_edge<<std::endl;
    os<<"N_Rejected for too few pixels = "<<nf_npix<<std::endl;

    os<<"PSF measurement:\n";
    os<<"  N_Success = "<<ns_psf<<std::endl;
    os<<"  N_Fail = "<<nf_psf<<std::endl;

    os<<"N_Error: TMV Error caught = "<<nf_tmverror<<std::endl;
    os<<"N_Error: Other caught = "<<nf_othererror<<std::endl;
  }

  ExitCode exitcode;

  int nstars;
  int nf_range;
  int nf_edge;
  int nf_npix;
  int nf_tmverror;
  int nf_othererror;
  int ns_psf;
  int nf_psf;

  mutable std::ofstream logout;

};

inline std::ostream& operator<<(std::ostream& os, const ShearLog& log)
{ log.Write(os); return os; }

inline std::ostream& operator<<(std::ostream& os, const PSFLog& log)
{ log.Write(os); return os; }

#endif
