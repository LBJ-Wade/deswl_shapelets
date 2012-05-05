#ifndef LOG_H
#define LOG_H

#include <fstream>
#include "dbg.h"
#include "Params.h"
#include "ConfigFile.h"

class Log 
{
public :

    // Default Log constructor doesn't write anything anywhere.
    Log() : _params(0), _logout(0), _fits_file("") {}

    Log(const ConfigFile& params,
        std::string log_file="", std::string fits_file="");
    virtual ~Log();

    void noWriteLog(); // Don't write the log on deletion.

    virtual void writeMessage(std::string message) const;
    virtual void writeLog() const = 0;
    virtual void writeLogToFitsHeader() const = 0;
    virtual void write(std::ostream& os) const = 0;

    void setExitCode(ExitCode ec, const std::string& s="")
    { 
        _exit_code = ec; 
        _extra_exit_info = s;
    }

protected :

    ExitCode _exit_code;
    std::string _extra_exit_info;

    const ConfigFile* _params;
    std::ostream* _logout;
    std::string _fits_file;

};

class ShearLog : public Log 
{
public :

    ShearLog() {}
    ShearLog(const ConfigFile& params,
             std::string log_file="", std::string fits_file="");
    virtual ~ShearLog();

    virtual void writeLog() const;
    virtual void writeLogToFitsHeader() const;
    virtual void write(std::ostream& os) const;

    ShearLog& operator+=(const ShearLog& rhs)
    {
        _nf_range1 += rhs._nf_range1;
        _nf_range2 += rhs._nf_range2;
        _nf_small += rhs._nf_small;
        _nf_tmv_error += rhs._nf_tmv_error;
        _nf_other_error += rhs._nf_other_error;
        _ns_centroid += rhs._ns_centroid;
        _nf_centroid += rhs._nf_centroid;
        _ns_native += rhs._ns_native;
        _nf_native += rhs._nf_native;
        _ns_mu += rhs._ns_mu;
        _nf_mu += rhs._nf_mu;
        _ns_gamma += rhs._ns_gamma;
        _nf_gamma += rhs._nf_gamma;

        return *this;
    }

    int _ngals;
    int _ngoodin;
    int _ngood;
    int _nf_range1;
    int _nf_range2;
    int _nf_small;
    int _nf_tmv_error;
    int _nf_other_error;
    int _ns_centroid;
    int _nf_centroid;
    int _ns_native;
    int _nf_native;
    int _ns_mu;
    int _nf_mu;
    int _ns_gamma;
    int _nf_gamma;

};

class PsfLog : public Log 
{
public :

    PsfLog() {}
    PsfLog(const ConfigFile& params,
           std::string log_file="", std::string fits_file="");
    virtual ~PsfLog();

    virtual void writeLog() const;
    virtual void writeLogToFitsHeader() const;
    virtual void write(std::ostream& os) const;

    PsfLog& operator+=(const PsfLog& rhs)
    {
        _nf_range += rhs._nf_range;
        _nf_tmv_error += rhs._nf_tmv_error;
        _nf_other_error += rhs._nf_other_error;
        _ns_psf += rhs._ns_psf;
        _nf_psf += rhs._nf_psf;

        return *this;
    }

    int _nstars;
    int _ngoodin;
    int _ngood;
    int _nf_range;
    int _nf_tmv_error;
    int _nf_other_error;
    int _ns_psf;
    int _nf_psf;
    int _noutliers;
    int _nfit;
    int _npc;
    double _chisq_fit;
    int _dof_fit;
};

class FindStarsLog : public Log 
{
public :

    FindStarsLog() {}
    FindStarsLog(const ConfigFile& params,
                 std::string log_file="", std::string fits_file="");
    virtual ~FindStarsLog();

    virtual void writeLog() const;
    virtual void writeLogToFitsHeader() const;
    virtual void write(std::ostream& os) const;

    int _ntot;
    int _nr_flag;
    int _nr_mag;
    int _nsg;
    int _nr_size;
    int _nobj;
    int _nallstars;
    int _nstars;
};

inline std::ostream& operator<<(std::ostream& os, const Log& log)
{ log.write(os); return os; }

#endif
