#ifndef LOG_H
#define LOG_H

#include "Params.h"
#include <fstream>
#include "ConfigFile.h"

class Log 
{
public :

    Log(const ConfigFile& params,
        std::string logFile="", std::string fitsFile="");
    virtual ~Log();

    void noWriteLog(); // Don't write the log on deletion.

    virtual void writeMessage(std::string message) const;
    virtual void writeLog() const = 0;
    virtual void writeLogToFitsHeader() const = 0;
    virtual void write(std::ostream& os) const = 0;

    void setExitCode(ExitCode ec, const std::string& s="")
    { 
        _exitCode = ec; 
        _extraExitInfo = s;
    }

protected :

    ExitCode _exitCode;
    std::string _extraExitInfo;

    const ConfigFile& _params;
    std::ostream* _logOut;
    std::string _fitsFile;

};

class ShearLog : public Log 
{
public :

    ShearLog(const ConfigFile& params,
             std::string logFile="", std::string fitsFile="");
    virtual ~ShearLog();

    virtual void writeLog() const;
    virtual void writeLogToFitsHeader() const;
    virtual void write(std::ostream& os) const;

    ShearLog& operator+=(const ShearLog& rhs)
    {
        _nfRange1 += rhs._nfRange1;
        _nfRange2 += rhs._nfRange2;
        _nfSmall += rhs._nfSmall;
        _nfTmvError += rhs._nfTmvError;
        _nfOtherError += rhs._nfOtherError;
        _nsCentroid += rhs._nsCentroid;
        _nfCentroid += rhs._nfCentroid;
        _nsNative += rhs._nsNative;
        _nfNative += rhs._nfNative;
        _nsMu += rhs._nsMu;
        _nfMu += rhs._nfMu;
        _nsGamma += rhs._nsGamma;
        _nfGamma += rhs._nfGamma;

        return *this;
    }

    int _nGals;
    int _nGoodIn;
    int _nGood;
    int _nfRange1;
    int _nfRange2;
    int _nfSmall;
    int _nfTmvError;
    int _nfOtherError;
    int _nsCentroid;
    int _nfCentroid;
    int _nsNative;
    int _nfNative;
    int _nsMu;
    int _nfMu;
    int _nsGamma;
    int _nfGamma;
    int _nsGood;

};

class PsfLog : public Log 
{
public :

    PsfLog(const ConfigFile& params,
           std::string logFile="", std::string fitsFile="");
    virtual ~PsfLog();

    virtual void writeLog() const;
    virtual void writeLogToFitsHeader() const;
    virtual void write(std::ostream& os) const;

    PsfLog& operator+=(const PsfLog& rhs)
    {
        _nfRange += rhs._nfRange;
        _nfTmvError += rhs._nfTmvError;
        _nfOtherError += rhs._nfOtherError;
        _nsPsf += rhs._nsPsf;
        _nfPsf += rhs._nfPsf;

        return *this;
    }

    int _nStars;
    int _nGoodIn;
    int _nGood;
    int _nfRange;
    int _nfTmvError;
    int _nfOtherError;
    int _nsPsf;
    int _nfPsf;
    int _nsGood;

};

class FindStarsLog : public Log 
{
public :

    FindStarsLog(const ConfigFile& params,
                 std::string logFile="", std::string fitsFile="");
    virtual ~FindStarsLog();

    virtual void writeLog() const;
    virtual void writeLogToFitsHeader() const;
    virtual void write(std::ostream& os) const;

    int _nTot;
    int _nrFlag;
    int _nrMag;
    int _nrSize;
    int _nObj;
    int _nAllStars;
    int _nStars;
};

inline std::ostream& operator<<(std::ostream& os, const Log& log)
{ log.write(os); return os; }

#endif
