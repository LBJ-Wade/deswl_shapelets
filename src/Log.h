#ifndef LOG_H
#define LOG_H

#include "Params.h"
#include "FitsFile.h"
#include <fstream>

struct ShearLog {

  ShearLog(std::string logfile="", std::string fits_file="");
  ~ShearLog();

  ShearLog& operator+=(const ShearLog& rhs)
  {
    nf_range1 += rhs.nf_range1;
    nf_range2 += rhs.nf_range2;
    nf_edge1 += rhs.nf_edge1;
    nf_npix1 += rhs.nf_npix1;
    nf_small += rhs.nf_small;
    nf_edge2 += rhs.nf_edge2;
    nf_npix2 += rhs.nf_npix2;
    nf_tmverror += rhs.nf_tmverror;
    nf_othererror += rhs.nf_othererror;
    ns_native += rhs.ns_native;
    nf_native += rhs.nf_native;
    ns_mu += rhs.ns_mu;
    nf_mu += rhs.nf_mu;
    ns_gamma += rhs.ns_gamma;
    nf_gamma += rhs.nf_gamma;
    
    return *this;
  }

  void NoWriteLog(); // Don't write the log on deletion.
  void WriteLog() const;
  void WriteLogToFitsHeader() const;
  void Write(std::ostream& os) const;

  ExitCode exitcode;
  std::string extraexitinfo;

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

  std::ostream* logout;
  std::string fits_file;

};

struct PSFLog {

  PSFLog(std::string logfile="", std::string fits_file="");
  ~PSFLog();

  PSFLog& operator+=(const PSFLog& rhs)
  {
    nf_range += rhs.nf_range;
    nf_edge += rhs.nf_edge;
    nf_npix += rhs.nf_npix;
    nf_tmverror += rhs.nf_tmverror;
    nf_othererror += rhs.nf_othererror;
    ns_psf += rhs.ns_psf;
    nf_psf += rhs.nf_psf;

    return *this;
  }

  void NoWriteLog(); // Don't write the log on deletion.
  void WriteLog() const;
  void WriteLogToFitsHeader() const;
  void Write(std::ostream& os) const;

  ExitCode exitcode;
  std::string extraexitinfo;

  int nstars;
  int nf_range;
  int nf_edge;
  int nf_npix;
  int nf_tmverror;
  int nf_othererror;
  int ns_psf;
  int nf_psf;

  private :

  std::ostream* logout;
  std::string fits_file;

};

struct FindStarsLog {

  FindStarsLog(std::string _logfile="", std::string _fits_file="");
  ~FindStarsLog();

  void NoWriteLog(); // Don't write the log on deletion.
  void WriteLog() const;
  void WriteLogToFitsHeader() const;
  void Write(std::ostream& os) const;

  ExitCode exitcode;
  std::string extraexitinfo;

  int nobj;
  int nstars;

  private :

  std::ostream* logout;
  std::string fits_file;

};

inline std::ostream& operator<<(std::ostream& os, const ShearLog& log)
{ log.Write(os); return os; }

inline std::ostream& operator<<(std::ostream& os, const PSFLog& log)
{ log.Write(os); return os; }

inline std::ostream& operator<<(std::ostream& os, const FindStarsLog& log)
{ log.Write(os); return os; }

#endif
