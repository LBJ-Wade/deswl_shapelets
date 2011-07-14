
#include <iostream>
#include <CCfits/CCfits>
#include "Log.h"
#include "Params.h"
#include "Name.h"

Log::Log(
    const ConfigFile& params,
    std::string logFile, std::string fitsFile
) :
    _exitCode(SUCCESS), _extraExitInfo(""),
    _params(params), _logOut(0), _fitsFile(fitsFile)
{
    bool shouldOutputDesQa = params.read("des_qa",false);
    if (shouldOutputDesQa) { 
        if (logFile == "") {
            _logOut = &std::cout;
        } else {
            _logOut = new std::ofstream(logFile.c_str(),std::ios_base::app);
            if (!_logOut) {
                throw WriteException(
                    "Error: Unable to open logfile " + logFile +
                    " for append output");
            }
        }
    }
}

Log::~Log() 
{
#if 0
    this->writeLog();
    this->writeLogToFitsHeader();
    // I thought this would work, but it doesn't. The reason is a bit subtle.  
    // Basically, the destructors get called in order from most derived 
    // up to the base.
    // So by this point, the derived portions have already been deleted,
    // so the derived implementations of these functions are already gone.
    // The solution is to put this same thing in each derived destructor.
#endif
    noWriteLog(); // deletes _logOut
}


void Log::writeMessage(std::string message) const
{
    if (_logOut) {
        *_logOut << message;
    }
}

void Log::noWriteLog() 
{
    if (_logOut && _logOut != &std::cout) {
        delete _logOut;
    }
    _logOut = 0;
}

ShearLog::ShearLog(
    const ConfigFile& params,
    std::string logFile, std::string fitsFile
) : 
    Log(params,logFile,fitsFile),
    _nGals(0), _nGoodIn(0), _nGood(0), _nfRange1(0), _nfRange2(0), 
    _nfSmall(0), _nfTmvError(0), _nfOtherError(0), 
    _nsCentroid(0), _nfCentroid(0), _nsNative(0), _nfNative(0),
    _nsMu(0), _nfMu(0), _nsGamma(0), _nfGamma(0)
{}

ShearLog::~ShearLog() 
{
    this->writeLog();
    this->writeLogToFitsHeader();
}

void ShearLog::writeLogToFitsHeader() const
{
    if ("" == _fitsFile) return;
    // Also return if the provided file isn't actually a fits file.
    if (_fitsFile.find(".fits") == std::string::npos) return;

    try {
        int hdu = getHdu(this->_params,"shear",_fitsFile,2);
        CCfits::FITS fits(_fitsFile, CCfits::Write, hdu-1);

        CCfits::ExtHDU& table=fits.extension(hdu-1);

        // the enumerated type ExitCode sometimes defaults to long
        table.addKey("msexit", (int) _exitCode, "Exit code for MeasureShear");
        table.addKey("msnobj", _nGals, 
                     "# of objects processed by MeasureShear");
        table.addKey("ms_ngoodin", _nGoodIn,
                     "# of objects with no input error flags");
        table.addKey("ms_ngood", _nGood,
                     "# of measurements with no error flags");
        table.addKey("ns_centr", _nsCentroid, 
                     "# of objects successfully centroided");
        table.addKey("ns_gamma", _nsGamma, 
                     "# of successful shear measurements");
        table.addKey("ns_nativ", _nsNative, 
                     "# of successful native shapelet measurements");

        table.addKey("nf_rnge1", _nfRange1, 
                     "# of MeasureShear failures range1");
        table.addKey("nf_rnge2", _nfRange2, 
                     "# of MeasureShear failures range2");

        table.addKey("nf_centr", _nfCentroid, 
                     "# of MeasureShear failures during centroiding");
        table.addKey("nf_nativ", _nfNative, 
                     "# of MeasureShear failures native calculations");
        table.addKey("nf_small", _nfSmall, 
                     "# of MeasureShear failures too small");

        // prefix ms to make unique
        table.addKey("nf_mstmv", _nfTmvError, 
                     "# of MeasureShear failures TMV errors");
        table.addKey("nf_msoth", _nfOtherError, 
                     "# of MeasureShear failures unclassified errors");

        table.addKey("nf_mu", _nfMu, 
                     "# of MeasureShear failures calculating mu");
        table.addKey("nf_gamma", _nfGamma, 
                     "# of MeasureShear failures calculating shear");
    } catch (std::exception& e) {
        xdbg<<"Caught exception during Log write:\n";
        xdbg<<e.what()<<std::endl;
        xdbg<<"Ignoring, and moving on... "<<std::endl;
    } catch (...) {
        xdbg<<"Caught exception during Log write.\n";
        xdbg<<"Ignoring, and moving on... "<<std::endl;
    }
}

inline std::string StatusText1(ExitCode exitCode, const ConfigFile& params)
{
    return std::string("STATUS") +
        char('0'+Status(exitCode,params)) +
        std::string("BEG"); 
}
inline std::string StatusText2(ExitCode exitCode, const ConfigFile& params)
{
    return std::string("STATUS") +
        char('0'+Status(exitCode,params)) +
        std::string("END");
}

void ShearLog::writeLog() const
{
    if (_logOut) {
        // Emit logging information
        if (_exitCode) {
            *_logOut << 
                StatusText1(_exitCode,this->_params)<<" "<<
                Text(_exitCode)<<" "<<
                _extraExitInfo<<" ";
            if (_fitsFile != "")
                *_logOut << " (Name="<<_fitsFile<<") ";
            *_logOut<< StatusText2(_exitCode,this->_params)<<std::endl;
        } else {
            std::string name = "measureshear";
            if (_fitsFile != "") name = _fitsFile;
            *_logOut << 
                "QA2BEG "<<
                "Name="<<name <<" & "<<
                "ngals="<<_nGals <<" & "<<
                "ngoodin="<<_nGoodIn <<" & "<<
                "ngood="<<_nGood <<" & "<<
                "ns_centroid="<<_nsCentroid <<" & "<<
                "ns_gamma="<<_nsGamma <<" & "<<
                "ns_native="<<_nsNative <<" & "<<
                "nf_range1="<<_nfRange1 <<" & "<<
                "nf_range2="<<_nfRange2 <<" & "<<
                "nf_centroid="<<_nfCentroid <<" & "<<
                "nf_native="<<_nfNative <<" & "<<
                "nf_small="<<_nfSmall <<" & "<<
                "nf_tmverror="<<_nfTmvError <<" & "<<
                "nf_othererror="<<_nfOtherError <<" & "<<
                "nf_mu="<<_nfMu <<" & "<<
                "nf_gamma="<<_nfGamma <<
                " QA2END"<<std::endl;
        }
    }
}

void ShearLog::write(std::ostream& os) const
{
    os<<"N Total objects = "<<_nGals<<std::endl;
    os<<"N objects with no input error flags = "<<_nGoodIn<<std::endl;
    os<<"N objects with no measurement error flags = "<<_nGood<<std::endl;

    os<<"N_Rejected for range error - distortion = "<<_nfRange1<<std::endl;
    os<<"N_Rejected for range error - psf interp = "<<_nfRange2<<std::endl;

    os<<"Centroid step:\n";
    os<<"  N_Success = "<<_nsCentroid<<std::endl;
    os<<"  N_Fail = "<<_nfCentroid<<std::endl;

    os<<"Native fits:\n";
    os<<"  N_Success = "<<_nsNative<<std::endl;
    os<<"  N_Fail = "<<_nfNative<<std::endl;

    os<<"N_Rejected for being too small = "<<_nfSmall<<std::endl;

    os<<"Mu fits:\n";
    os<<"  N_Success = "<<_nsMu<<std::endl;
    os<<"  N_Fail = "<<_nfMu<<std::endl;

    os<<"Gamma fits:\n";
    os<<"  N_Success = "<<_nsGamma<<std::endl;
    os<<"  N_Fail = "<<_nfGamma<<std::endl;

    os<<"N_Error: TMV Error caught = "<<_nfTmvError<<std::endl;
    os<<"N_Error: Other caught = "<<_nfOtherError<<std::endl;
}

PsfLog::PsfLog(
    const ConfigFile& params,
    std::string logFile, std::string fitsFile
) :
    Log(params,logFile,fitsFile),
    _nStars(0), _nGoodIn(0), _nGood(0), _nfRange(0),
    _nfTmvError(0), _nfOtherError(0), _nsPsf(0), _nfPsf(0),
    _nOutliers(0), _nFit(0), _nPC(0), _chisqFit(0.), _dofFit(0) {}

PsfLog::~PsfLog() 
{
    this->writeLog();
    this->writeLogToFitsHeader();
}

void PsfLog::writeLogToFitsHeader() const
{
    if ("" == _fitsFile) return;
    // Also return if the provided file isn't actually a fits file.
    if (_fitsFile.find(".fits") == std::string::npos) return;

    try {
        int hdu = getHdu(this->_params,"psf",_fitsFile,2);
        CCfits::FITS fits(_fitsFile, CCfits::Write, hdu-1);

        CCfits::ExtHDU& table=fits.extension(hdu-1);

        // the enumerated type ExitCode sometimes defaults to long
        table.addKey("mpexit", (int) _exitCode, "Exit code for MeasurePSF");
        table.addKey("mpnstars", _nStars, 
                     "# of PSF star candidates used");
        table.addKey("mpgoodin", _nGoodIn,
                     "# of objects with no input error flags");
        table.addKey("mp_ngood", _nGood,
                     "# of measurements with no error flags");
        table.addKey("ns_psf", _nsPsf, 
                     "# of successful PSF decompositions");
        table.addKey("nf_range", _nfRange, 
                     "# PSF failures due to range error");
        table.addKey("nf_mptmv", _nfTmvError, 
                     "# PSF failures due to TMV errors");
        table.addKey("nf_mpoth", _nfOtherError, 
                     "# PSF failures due to unclassified errors");
        table.addKey("nf_psf", _nfPsf, 
                     "# PSF failures");
        table.addKey("noutlier", _nOutliers,
                     "# outliers rejected during fit");
        table.addKey("n_fit", _nFit,
                     "# final stars used for fit");
        table.addKey("n_pc", _nPC,
                     "# principal components used for fit");
        table.addKey("chisqfit", _chisqFit,
                     "# chisq of fit");
        table.addKey("dof_fit", _dofFit,
                     "# degrees of freedom of fit");
    } catch (std::exception& e) {
        xdbg<<"Caught exception during Log write:\n";
        xdbg<<e.what()<<std::endl;
        xdbg<<"Ignoring, and moving on... "<<std::endl;
    } catch (...) {
        xdbg<<"Caught exception during Log write.\n";
        xdbg<<"Ignoring, and moving on... "<<std::endl;
    }
}

void PsfLog::writeLog() const
{
    if (_logOut) {
        // Emit logging information
        if (_exitCode) {
            *_logOut << 
                StatusText1(_exitCode,this->_params)<<" "<<
                Text(_exitCode)<<" "<<
                _extraExitInfo<<" ";
            if (_fitsFile != "")
                *_logOut << " (Name="<<_fitsFile<<") ";
            *_logOut<< StatusText2(_exitCode,this->_params)<<std::endl;
        } else {
            std::string name = "measurepsf";
            if (_fitsFile != "") name = _fitsFile;
            *_logOut << 
                "QA2BEG "<<
                "Name="<<name <<" & "<<
                "nstars="<<_nStars <<" & "<<
                "ngoodin="<<_nGoodIn <<" & "<<
                "ngood="<<_nGood <<" & "<<
                "ns_psf="<<_nsPsf <<" & "<<
                "nf_range="<<_nfRange <<" & "<<
                "nf_tmverror="<<_nfTmvError <<" & "<<
                "nf_othererror="<<_nfOtherError <<" & "<<
                "nf_psf="<<_nfPsf <<" & "<<
                "n_outlier="<<_nOutliers <<" & "<<
                "n_fit="<<_nFit <<" & "<<
                "n_pc="<<_nPC <<" & "<<
                "chisq_fit="<<_chisqFit <<" & "<<
                "dof_fit="<<_dofFit <<
                " QA2END"<<std::endl;
        }
    }
}

void PsfLog::write(std::ostream& os) const
{
    os<<"N Total stars = "<<_nStars<<std::endl;
    os<<"N stars with no input error flags = "<<_nGoodIn<<std::endl;
    os<<"N stars with no measurement error flags = "<<_nGood<<std::endl;

    os<<"N Rejected for range error - distortion = "<<_nfRange<<std::endl;

    os<<"PSF measurement:\n";
    os<<"  N Success = "<<_nsPsf<<std::endl;
    os<<"  N Fail = "<<_nfPsf<<std::endl;

    os<<"PSF interpolation:\n";
    os<<"  N outliers found = "<<_nOutliers<<std::endl;
    os<<"  N stars used in fit = "<<_nFit<<std::endl;
    os<<"  N principal components used for fit = "<<_nPC<<std::endl;
    os<<"  chisq of fit = "<<_chisqFit<<std::endl;
    os<<"  degrees of freedom of fit = "<<_dofFit<<std::endl;

    os<<"N_Error: TMV Error caught = "<<_nfTmvError<<std::endl;
    os<<"N_Error: Other caught = "<<_nfOtherError<<std::endl;
}


FindStarsLog::FindStarsLog(
    const ConfigFile& params,
    std::string logFile, std::string fitsFile
) : 
    Log(params,logFile,fitsFile),
    _nTot(0), _nrFlag(0), _nrMag(0), _nSg(0), _nrSize(0),
    _nObj(0), _nAllStars(0), _nStars(0) 
{}

FindStarsLog::~FindStarsLog() 
{
    this->writeLog();
    this->writeLogToFitsHeader();
}

void FindStarsLog::writeLogToFitsHeader() const
{
    if ("" == _fitsFile) return;
    // Also return if the provided file isn't actually a fits file.
    if (_fitsFile.find(".fits") == std::string::npos) return;

    try {
        int hdu = getHdu(this->_params,"stars",_fitsFile,2);
        CCfits::FITS fits(_fitsFile, CCfits::Write, hdu-1);

        CCfits::ExtHDU& table=fits.extension(hdu-1);

        // the enumerated type ExitCode sometimes defaults to long
        table.addKey("fsexit", (int) _exitCode, "Exit code for FindStars");
        table.addKey("fsntot", _nTot, 
                     "# of total objects processed in FindStars");
        table.addKey("fsnrflag", _nrFlag,
                     "# of objects rejected because of input flag");
        table.addKey("fsnrmag", _nrMag,
                     "# of objects rejected because of mag limits");
        table.addKey("fsnsg", _nSg,
                     "# of objects in inital star-galaxy,mag limits");
        table.addKey("fsnrsize", _nrSize,
                     "# of objects rejected because of size limits");
        table.addKey("fsnobj", _nObj, 
                     "# of objects processed in FindStars");
        table.addKey("fsnall", _nAllStars,
                     "# of total stars found by FindStars");
        table.addKey("fsnstars", _nStars, 
                     "# of good stars found by FindStars");
    } catch (std::exception& e) {
        xdbg<<"Caught exception during Log write:\n";
        xdbg<<e.what()<<std::endl;
        xdbg<<"Ignoring, and moving on... "<<std::endl;
    } catch (...) {
        xdbg<<"Caught exception during Log write.\n";
        xdbg<<"Ignoring, and moving on... "<<std::endl;
    }
}

void FindStarsLog::writeLog() const
{
    if (_logOut) {
        // Emit logging information
        if (_exitCode) {
            *_logOut << 
                StatusText1(_exitCode,this->_params)<<" "<<
                Text(_exitCode)<<" "<<
                _extraExitInfo<<" ";
            if (_fitsFile != "")
                *_logOut << " (Name="<<_fitsFile<<") ";
            *_logOut<< StatusText2(_exitCode,this->_params)<<std::endl;
        } else {
            std::string name = "findstars";
            if (_fitsFile != "") name = _fitsFile;
            *_logOut << 
                "QA2BEG "<<
                "Name="<<name <<" & "<<
                "ntot="<<_nTot <<" & "<<
                "nr_flag="<<_nrFlag <<" & "<<
                "nr_mag="<<_nrMag <<" & "<<
                "nsg="<<_nSg <<" & "<<
                "nr_size="<<_nrSize <<" & "<<
                "nobj="<<_nObj <<" & "<<
                "nallstars="<<_nAllStars <<" & "<<
                "nstars="<<_nStars <<
                " QA2END"<<std::endl;
        }
    }
}

void FindStarsLog::write(std::ostream& os) const
{
    os<<"N total input objects = "<<_nTot<<std::endl;
    os<<"N rejected because of input flag = "<<_nrFlag<<std::endl;
    os<<"N rejected because of mag limits = "<<_nrMag<<std::endl;
    os<<"N rejected because of size limits = "<<_nrSize<<std::endl;
    os<<"N objects processed in FindStars = "<<_nObj<<std::endl;
    os<<"N total stars found by FindStars = "<<_nAllStars<<std::endl;
    os<<"N good stars found by FindStars = "<<_nStars<<std::endl;
}


