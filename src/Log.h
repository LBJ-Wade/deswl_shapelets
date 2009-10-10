#ifndef LOG_H
#define LOG_H

#include "Params.h"
#include <fstream>
#include "ConfigFile.h"

class Log 
{
  public :

    Log(const ConfigFile& params,
	std::string logfile="", std::string fits_file="");
    virtual ~Log();

    void NoWriteLog(); // Don't write the log on deletion.

	virtual void WriteMessage(std::string mess) const=0;
    virtual void WriteLog() const = 0;
    virtual void WriteLogToFitsHeader() const = 0;
    virtual void Write(std::ostream& os) const = 0;

    ExitCode exitcode;
    std::string extraexitinfo;

  protected :

    const ConfigFile& params;
    std::ostream* logout;
    std::string fits_file;

};

class ShearLog : public Log 
{
  public :

    ShearLog(const ConfigFile& params,
	std::string logfile="", std::string fits_file="");
    virtual ~ShearLog();

	virtual void WriteMessage(std::string mess) const;
    virtual void WriteLog() const;
    virtual void WriteLogToFitsHeader() const;
    virtual void Write(std::ostream& os) const;

    ShearLog& operator+=(const ShearLog& rhs)
    {
      nf_range1 += rhs.nf_range1;
      nf_range2 += rhs.nf_range2;
      nf_small += rhs.nf_small;
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

    int ngals;
    int ngoodin;
    int ngood;
    int nf_range1;
    int nf_range2;
    int nf_small;
    int nf_tmverror;
    int nf_othererror;
    int ns_native;
    int nf_native;
    int ns_mu;
    int nf_mu;
    int ns_gamma;
    int nf_gamma;
    int ns_good;

};

class PSFLog : public Log 
{
  public :

    PSFLog(const ConfigFile& params,
	std::string logfile="", std::string fits_file="");
    virtual ~PSFLog();

	virtual void WriteMessage(std::string mess) const;
    virtual void WriteLog() const;
    virtual void WriteLogToFitsHeader() const;
    virtual void Write(std::ostream& os) const;

    PSFLog& operator+=(const PSFLog& rhs)
    {
      nf_range += rhs.nf_range;
      nf_tmverror += rhs.nf_tmverror;
      nf_othererror += rhs.nf_othererror;
      ns_psf += rhs.ns_psf;
      nf_psf += rhs.nf_psf;

      return *this;
    }

    int nstars;
    int ngoodin;
    int ngood;
    int nf_range;
    int nf_tmverror;
    int nf_othererror;
    int ns_psf;
    int nf_psf;
    int ns_good;

};

class FindStarsLog : public Log 
{
  public :

    FindStarsLog(const ConfigFile& params,
	std::string logfile="", std::string fits_file="");
    virtual ~FindStarsLog();

	virtual void WriteMessage(std::string mess) const;
    virtual void WriteLog() const;
    virtual void WriteLogToFitsHeader() const;
    virtual void Write(std::ostream& os) const;

    int ntot;
    int nr_flag;
    int nr_mag;
    int nr_size;
    int nobj;
    int nallstars;
    int nstars;

};

inline std::ostream& operator<<(std::ostream& os, const Log& log)
{ log.Write(os); return os; }

#endif
