
#include <iostream>
#include <CCfits/CCfits>
#include "Log.h"
#include "Params.h"
#include "Name.h"

Log::Log(
    const ConfigFile& params,
    std::string log_file, std::string fits_file
) :
    _exit_code(SUCCESS), _extra_exit_info(""),
    _params(&params), _logout(0), _fits_file(fits_file)
{
    bool des_qa = params.read("des_qa",false);
    if (des_qa) { 
        if (log_file == "") {
            _logout = &std::cout;
        } else {
            _logout = new std::ofstream(log_file.c_str(),std::ios_base::app);
            if (!_logout) {
                throw WriteException(
                    "Error: Unable to open logfile " + log_file +
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
    noWriteLog(); // deletes _logout
}


void Log::writeMessage(std::string message) const
{
    if (_logout) {
        *_logout << message;
    }
}

void Log::noWriteLog() 
{
    if (_logout && _logout != &std::cout) {
        delete _logout;
    }
    _logout = 0;
}

ShearLog::ShearLog( const ConfigFile& params, std::string log_file, std::string fits_file) : 
    Log(params,log_file,fits_file),
    _ngals(0), _ngoodin(0), _ngood(0), _nf_range1(0), _nf_range2(0), 
    _nf_small(0), _nf_tmv_error(0), _nf_other_error(0), 
    _ns_centroid(0), _nf_centroid(0), _ns_native(0), _nf_native(0),
    _ns_mu(0), _nf_mu(0), _ns_gamma(0), _nf_gamma(0)
{}

ShearLog::~ShearLog() 
{
    this->writeLog();
    this->writeLogToFitsHeader();
}

void ShearLog::writeLogToFitsHeader() const
{
    if ("" == _fits_file) return;
    // Also return if the provided file isn't actually a fits file.
    if (_fits_file.find(".fits") == std::string::npos) return;

    try {
        int hdu = GetHdu(*this->_params,"shear",_fits_file,2);
        CCfits::FITS fits(_fits_file, CCfits::Write, hdu-1);

        CCfits::ExtHDU& table=fits.extension(hdu-1);

        // the enumerated type ExitCode sometimes defaults to long
        table.addKey("msexit", (int) _exit_code, "Exit code for MeasureShear");
        table.addKey("msnobj", _ngals, 
                     "# of objects processed by MeasureShear");
        table.addKey("ms_ngoodin", _ngoodin,
                     "# of objects with no input error flags");
        table.addKey("ms_ngood", _ngood,
                     "# of measurements with no error flags");
        table.addKey("ns_centr", _ns_centroid, 
                     "# of objects successfully centroided");
        table.addKey("ns_gamma", _ns_gamma, 
                     "# of successful shear measurements");
        table.addKey("ns_nativ", _ns_native, 
                     "# of successful native shapelet measurements");

        table.addKey("nf_rnge1", _nf_range1, 
                     "# of MeasureShear failures range1");
        table.addKey("nf_rnge2", _nf_range2, 
                     "# of MeasureShear failures range2");

        table.addKey("nf_centr", _nf_centroid, 
                     "# of MeasureShear failures during centroiding");
        table.addKey("nf_nativ", _nf_native, 
                     "# of MeasureShear failures native calculations");
        table.addKey("nf_small", _nf_small, 
                     "# of MeasureShear failures too small");

        // prefix ms to make unique
        table.addKey("nf_mstmv", _nf_tmv_error, 
                     "# of MeasureShear failures TMV errors");
        table.addKey("nf_msoth", _nf_other_error, 
                     "# of MeasureShear failures unclassified errors");

        table.addKey("nf_mu", _nf_mu, 
                     "# of MeasureShear failures calculating mu");
        table.addKey("nf_gamma", _nf_gamma, 
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

inline std::string StatusText1(ExitCode exit_code, const ConfigFile& params)
{
    return std::string("STATUS") +
        char('0'+Status(exit_code,params)) +
        std::string("BEG"); 
}
inline std::string StatusText2(ExitCode exit_code, const ConfigFile& params)
{
    return std::string("STATUS") +
        char('0'+Status(exit_code,params)) +
        std::string("END");
}

void ShearLog::writeLog() const
{
    if (_logout) {
        // Emit logging information
        if (_exit_code) {
            *_logout << 
                StatusText1(_exit_code,*this->_params)<<" "<<
                Text(_exit_code)<<" "<<
                _extra_exit_info<<" ";
            if (_fits_file != "")
                *_logout << " (Name="<<_fits_file<<") ";
            *_logout<< StatusText2(_exit_code,*this->_params)<<std::endl;
        } else {
            std::string name = "measureshear";
            if (_fits_file != "") name = _fits_file;
            *_logout << 
                "QA2BEG "<<
                "Name="<<name <<" & "<<
                "ngals="<<_ngals <<" & "<<
                "ngoodin="<<_ngoodin <<" & "<<
                "ngood="<<_ngood <<" & "<<
                "ns_centroid="<<_ns_centroid <<" & "<<
                "ns_gamma="<<_ns_gamma <<" & "<<
                "ns_native="<<_ns_native <<" & "<<
                "nf_range1="<<_nf_range1 <<" & "<<
                "nf_range2="<<_nf_range2 <<" & "<<
                "nf_centroid="<<_nf_centroid <<" & "<<
                "nf_native="<<_nf_native <<" & "<<
                "nf_small="<<_nf_small <<" & "<<
                "nf_tmverror="<<_nf_tmv_error <<" & "<<
                "nf_othererror="<<_nf_other_error <<" & "<<
                "nf_mu="<<_nf_mu <<" & "<<
                "nf_gamma="<<_nf_gamma <<
                " QA2END"<<std::endl;
        }
    }
}

void ShearLog::write(std::ostream& os) const
{
    os<<"N Total objects = "<<_ngals<<std::endl;
    os<<"N objects with no input error flags = "<<_ngoodin<<std::endl;
    os<<"N objects with no measurement error flags = "<<_ngood<<std::endl;

    os<<"N_Rejected for range error - distortion = "<<_nf_range1<<std::endl;
    os<<"N_Rejected for range error - psf interp = "<<_nf_range2<<std::endl;

    os<<"Centroid step:\n";
    os<<"  N_Success = "<<_ns_centroid<<std::endl;
    os<<"  N_Fail = "<<_nf_centroid<<std::endl;

    os<<"Native fits:\n";
    os<<"  N_Success = "<<_ns_native<<std::endl;
    os<<"  N_Fail = "<<_nf_native<<std::endl;

    os<<"N_Rejected for being too small = "<<_nf_small<<std::endl;

    os<<"Mu fits:\n";
    os<<"  N_Success = "<<_ns_mu<<std::endl;
    os<<"  N_Fail = "<<_nf_mu<<std::endl;

    os<<"Gamma fits:\n";
    os<<"  N_Success = "<<_ns_gamma<<std::endl;
    os<<"  N_Fail = "<<_nf_gamma<<std::endl;

    os<<"N_Error: TMV Error caught = "<<_nf_tmv_error<<std::endl;
    os<<"N_Error: Other caught = "<<_nf_other_error<<std::endl;
}

PsfLog::PsfLog(
    const ConfigFile& params,
    std::string log_file, std::string fits_file
) :
    Log(params,log_file,fits_file),
    _nstars(0), _ngoodin(0), _ngood(0), _nf_range(0),
    _nf_tmv_error(0), _nf_other_error(0), _ns_psf(0), _nf_psf(0),
    _noutliers(0), _nfit(0), _npc(0), _chisq_fit(0.), _dof_fit(0) {}

PsfLog::~PsfLog() 
{
    this->writeLog();
    this->writeLogToFitsHeader();
}

void PsfLog::writeLogToFitsHeader() const
{
    if ("" == _fits_file) return;
    // Also return if the provided file isn't actually a fits file.
    if (_fits_file.find(".fits") == std::string::npos) return;

    try {
        int hdu = GetHdu(*this->_params,"psf",_fits_file,2);
        CCfits::FITS fits(_fits_file, CCfits::Write, hdu-1);

        CCfits::ExtHDU& table=fits.extension(hdu-1);

        // the enumerated type ExitCode sometimes defaults to long
        table.addKey("mpexit", (int) _exit_code, "Exit code for MeasurePSF");
        table.addKey("mpnstars", _nstars, 
                     "# of PSF star candidates used");
        table.addKey("mpgoodin", _ngoodin,
                     "# of objects with no input error flags");
        table.addKey("mp_ngood", _ngood,
                     "# of measurements with no error flags");
        table.addKey("ns_psf", _ns_psf, 
                     "# of successful PSF decompositions");
        table.addKey("nf_range", _nf_range, 
                     "# PSF failures due to range error");
        table.addKey("nf_mptmv", _nf_tmv_error, 
                     "# PSF failures due to TMV errors");
        table.addKey("nf_mpoth", _nf_other_error, 
                     "# PSF failures due to unclassified errors");
        table.addKey("nf_psf", _nf_psf, 
                     "# PSF failures");
        table.addKey("noutlier", _noutliers,
                     "# outliers rejected during fit");
        table.addKey("n_fit", _nfit,
                     "# final stars used for fit");
        table.addKey("n_pc", _npc,
                     "# principal components used for fit");
        table.addKey("chisqfit", _chisq_fit,
                     "# chisq of fit");
        table.addKey("dof_fit", _dof_fit,
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
    if (_logout) {
        // Emit logging information
        if (_exit_code) {
            *_logout << 
                StatusText1(_exit_code,*this->_params)<<" "<<
                Text(_exit_code)<<" "<<
                _extra_exit_info<<" ";
            if (_fits_file != "")
                *_logout << " (Name="<<_fits_file<<") ";
            *_logout<< StatusText2(_exit_code,*this->_params)<<std::endl;
        } else {
            std::string name = "measurepsf";
            if (_fits_file != "") name = _fits_file;
            *_logout << 
                "QA2BEG "<<
                "Name="<<name <<" & "<<
                "nstars="<<_nstars <<" & "<<
                "ngoodin="<<_ngoodin <<" & "<<
                "ngood="<<_ngood <<" & "<<
                "ns_psf="<<_ns_psf <<" & "<<
                "nf_range="<<_nf_range <<" & "<<
                "nf_tmverror="<<_nf_tmv_error <<" & "<<
                "nf_othererror="<<_nf_other_error <<" & "<<
                "nf_psf="<<_nf_psf <<" & "<<
                "n_outlier="<<_noutliers <<" & "<<
                "n_fit="<<_nfit <<" & "<<
                "n_pc="<<_npc <<" & "<<
                "chisq_fit="<<_chisq_fit <<" & "<<
                "dof_fit="<<_dof_fit <<
                " QA2END"<<std::endl;
        }
    }
}

void PsfLog::write(std::ostream& os) const
{
    os<<"N Total stars = "<<_nstars<<std::endl;
    os<<"N stars with no input error flags = "<<_ngoodin<<std::endl;
    os<<"N stars with no measurement error flags = "<<_ngood<<std::endl;

    os<<"N Rejected for range error - distortion = "<<_nf_range<<std::endl;

    os<<"PSF measurement:\n";
    os<<"  N Success = "<<_ns_psf<<std::endl;
    os<<"  N Fail = "<<_nf_psf<<std::endl;

    os<<"PSF interpolation:\n";
    os<<"  N outliers found = "<<_noutliers<<std::endl;
    os<<"  N stars used in fit = "<<_nfit<<std::endl;
    os<<"  N principal components used for fit = "<<_npc<<std::endl;
    os<<"  chisq of fit = "<<_chisq_fit<<std::endl;
    os<<"  degrees of freedom of fit = "<<_dof_fit<<std::endl;

    os<<"N_Error: TMV Error caught = "<<_nf_tmv_error<<std::endl;
    os<<"N_Error: Other caught = "<<_nf_other_error<<std::endl;
}


FindStarsLog::FindStarsLog(
    const ConfigFile& params,
    std::string log_file, std::string fits_file
) : 
    Log(params,log_file,fits_file),
    _ntot(0), _nr_flag(0), _nr_mag(0), _nsg(0), _nr_size(0),
    _nobj(0), _nallstars(0), _nstars(0) 
{}

FindStarsLog::~FindStarsLog() 
{
    this->writeLog();
    this->writeLogToFitsHeader();
}

void FindStarsLog::writeLogToFitsHeader() const
{
    if ("" == _fits_file) return;
    // Also return if the provided file isn't actually a fits file.
    if (_fits_file.find(".fits") == std::string::npos) return;

    try {
        int hdu = GetHdu(*this->_params,"stars",_fits_file,2);
        CCfits::FITS fits(_fits_file, CCfits::Write, hdu-1);

        CCfits::ExtHDU& table=fits.extension(hdu-1);

        // the enumerated type ExitCode sometimes defaults to long
        table.addKey("fsexit", (int) _exit_code, "Exit code for FindStars");
        table.addKey("fsntot", _ntot, 
                     "# of total objects processed in FindStars");
        table.addKey("fsnrflag", _nr_flag,
                     "# of objects rejected because of input flag");
        table.addKey("fsnrmag", _nr_mag,
                     "# of objects rejected because of mag limits");
        table.addKey("fsnsg", _nsg,
                     "# of objects in inital star-galaxy,mag limits");
        table.addKey("fsnrsize", _nr_size,
                     "# of objects rejected because of size limits");
        table.addKey("fsnobj", _nobj, 
                     "# of objects processed in FindStars");
        table.addKey("fsnall", _nallstars,
                     "# of total stars found by FindStars");
        table.addKey("fsnstars", _nstars, 
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
    if (_logout) {
        // Emit logging information
        if (_exit_code) {
            *_logout << 
                StatusText1(_exit_code,*this->_params)<<" "<<
                Text(_exit_code)<<" "<<
                _extra_exit_info<<" ";
            if (_fits_file != "")
                *_logout << " (Name="<<_fits_file<<") ";
            *_logout<< StatusText2(_exit_code,*this->_params)<<std::endl;
        } else {
            std::string name = "findstars";
            if (_fits_file != "") name = _fits_file;
            *_logout << 
                "QA2BEG "<<
                "Name="<<name <<" & "<<
                "ntot="<<_ntot <<" & "<<
                "nr_flag="<<_nr_flag <<" & "<<
                "nr_mag="<<_nr_mag <<" & "<<
                "nsg="<<_nsg <<" & "<<
                "nr_size="<<_nr_size <<" & "<<
                "nobj="<<_nobj <<" & "<<
                "nallstars="<<_nallstars <<" & "<<
                "nstars="<<_nstars <<
                " QA2END"<<std::endl;
        }
    }
}

void FindStarsLog::write(std::ostream& os) const
{
    os<<"N total input objects = "<<_ntot<<std::endl;
    os<<"N rejected because of input flag = "<<_nr_flag<<std::endl;
    os<<"N rejected because of mag limits = "<<_nr_mag<<std::endl;
    os<<"N rejected because of size limits = "<<_nr_size<<std::endl;
    os<<"N objects processed in FindStars = "<<_nobj<<std::endl;
    os<<"N total stars found by FindStars = "<<_nallstars<<std::endl;
    os<<"N good stars found by FindStars = "<<_nstars<<std::endl;
}


