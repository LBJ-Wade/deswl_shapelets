#ifndef LOG_H
#define LOG_H

#include "Params.h"
#include <fstream>

struct ShearLog {

  ShearLog(std::string logfile, std::string _delim="  ");
  ~ShearLog();

  void WriteLog() const;
  void Write(std::ostream& os) const;

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

  private :

  mutable std::ofstream logout;
  std::string delim;

};

struct PSFLog {

  PSFLog(std::string logfile, std::string _delim="  ");
  ~PSFLog();

  void WriteLog() const;
  void Write(std::ostream& os) const;

  ExitCode exitcode;

  int nstars;
  int nf_range;
  int nf_edge;
  int nf_npix;
  int nf_tmverror;
  int nf_othererror;
  int ns_psf;
  int nf_psf;

  private :

  mutable std::ofstream logout;
  std::string delim;

};

inline std::ostream& operator<<(std::ostream& os, const ShearLog& log)
{ log.Write(os); return os; }

inline std::ostream& operator<<(std::ostream& os, const PSFLog& log)
{ log.Write(os); return os; }

#endif
