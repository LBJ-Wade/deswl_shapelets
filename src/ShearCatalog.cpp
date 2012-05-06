
#include <fstream>
#include <iostream>
#include <valarray>
#include <CCfits/CCfits>

#include "dbg.h"
#include "ShearCatalog.h"
#include "ConfigFile.h"
#include "Params.h"
#include "BVec.h"
#include "Ellipse.h"
#include "Pixel.h"
#include "Image.h"
#include "FittedPsf.h"
#include "Log.h"
#include "Params.h"
#include "Form.h"
#include "WlVersion.h"
#include "WriteParam.h"

ShearCatalog::ShearCatalog(
    const InputCatalog& incat, const Transformation& trans,
    const FittedPsf& fitpsf, const ConfigFile& params) :
    _id(incat.getIdList()), _pos(incat.getPosList()), 
    _sky(incat.getSkyList()), _noise(incat.getNoiseList()),
    _flags(incat.getFlagsList()), _skypos(incat.getSkyPosList()),
    _bounds(incat.getBounds()), _skybounds(incat.getSkyBounds()),
    _trans(&trans), _fitpsf(&fitpsf), _params(params)
{
    dbg<<"Create ShearCatalog\n";
    const int ngals = _id.size();

    // Fix flags to only have INPUT_FLAG set.
    // (I don't think this is necessary, but just in case.)
    for (int i=0; i<ngals; ++i) {
        if (_flags[i]) _flags[i] = INPUT_FLAG;
    }

    // Calculate sky positions
    Position skypos_default(DEFVALNEG,DEFVALNEG);
    if (_skypos.size() == 0) {
        _skypos.resize(ngals,skypos_default);
        for(int i=0;i<ngals;++i) {
            try {
                _trans->transform(_pos[i],_skypos[i]);
            } catch (RangeException& e) {
                xdbg<<"distortion range error\n";
                xdbg<<"p = "<<_pos[i]<<", b = "<<e.getBounds()<<std::endl;
            }                       
            if (!_flags[i]) _skybounds += _skypos[i];
        }
    } else {
        Assert(int(_skypos.size()) == ngals);

        xdbg<<"Check transformation:\n";
        double rms_error = 0.;
        int count = 0;
        for(int i=0;i<ngals;++i) {
            try {
                Position temp;
                _trans->transform(_pos[i],temp);
                if (dbgout && XDEBUG) {
                    dbgout->precision(12);
                    xdbg<<_pos[i]<<"  "<<_skypos[i]/3600.<<"  "<<temp/3600.;
                    dbgout->precision(6);
                    xdbg<<"  "<<temp-_skypos[i];
                }
                std::complex<double> diff = temp-_skypos[i];
                // Technically this is illegal since real doesn't have to 
                // return a reference.
                // gcc does, but pgCC doesn't.
                //real(diff) *= cos(imag(diff));
                double diff_ra = real(diff) * cos(imag(diff));
                diff = std::complex<double>(diff_ra,imag(diff));
                double err = std::norm(diff);
                if (err > 0.1) xdbg<<"  XXX";
                xdbg<<std::endl;
                rms_error += err;
                ++count;
            } catch (RangeException& e) {
                xdbg<<"distortion range error\n";
                xdbg<<"p = "<<_pos[i]<<", b = "<<e.getBounds()<<std::endl;
            }
        }
        rms_error /= count;
        rms_error = std::sqrt(rms_error);
        xdbg<<"rms error = "<<rms_error<<" arcsec\n";
        if (_params.read("des_qa",false)) {
            if (rms_error > 0.01) {
                std::cout<<
                    "STATUS3BEG Warning: Positions from WCS transformation "
                    "have rms error of "<<rms_error<<" arcsec relative to "
                    "ra, dec in catalog. STATUS3END"<<std::endl;
            }
        }
    }

    _shear.resize(ngals,std::complex<double>(DEFVALPOS,DEFVALPOS));
    _nu.resize(ngals,DEFVALNEG);

    DSmallMatrix22 cov_default;
    cov_default << DEFVALPOS, 0, 0, DEFVALPOS;
    _cov.resize(ngals,cov_default);

    int galorder = _params.read<int>("shear_gal_order");
    BVec shape_default(galorder,1.);
    shape_default.vec().TMV_setAllTo(DEFVALNEG);
    _meas_galorder.resize(ngals,galorder);
    _shape.resize(ngals,shape_default);

    Assert(int(_id.size()) == size());
    Assert(int(_pos.size()) == size());
    Assert(int(_sky.size()) == size());
    Assert(int(_noise.size()) == size());
    Assert(int(_flags.size()) == size());
    Assert(int(_skypos.size()) == size());
    Assert(int(_shear.size()) == size());
    Assert(int(_nu.size()) == size());
    Assert(int(_cov.size()) == size());
    Assert(int(_meas_galorder.size()) == size());
    Assert(int(_shape.size()) == size());
}

ShearCatalog::ShearCatalog(const ConfigFile& params) : _params(params)
{
    Assert(int(_id.size()) == size());
    Assert(int(_pos.size()) == size());
    Assert(int(_sky.size()) == size());
    Assert(int(_noise.size()) == size());
    Assert(int(_flags.size()) == size());
    Assert(int(_skypos.size()) == size());
    Assert(int(_shear.size()) == size());
    Assert(int(_nu.size()) == size());
    Assert(int(_cov.size()) == size());
    Assert(int(_meas_galorder.size()) == size());
    Assert(int(_shape.size()) == size());
}

void ShearCatalog::writeFits(std::string file) const
{
    Assert(int(_id.size()) == size());
    Assert(int(_pos.size()) == size());
    Assert(int(_sky.size()) == size());
    Assert(int(_noise.size()) == size());
    Assert(int(_flags.size()) == size());
    Assert(int(_shear.size()) == size());
    Assert(int(_nu.size()) == size());
    Assert(int(_cov.size()) == size());
    Assert(int(_meas_galorder.size()) == size());
    Assert(int(_shape.size()) == size());
    const int ngals = size();

    bool output_psf = _params.read("shear_output_psf",false);

    // ! means overwrite existing file
    CCfits::FITS fits("!"+file, CCfits::Write);

    const int nFields= output_psf ? 20 : 17;
    std::vector<string> col_names(nFields);
    std::vector<string> col_fmts(nFields);
    std::vector<string> col_units(nFields);

    std::vector<BVec> interp_psf;
    if (output_psf) {
        Assert(_fitpsf);
        interp_psf.resize(
            ngals,BVec(_fitpsf->getPsfOrder(),_fitpsf->getSigma()));
        for(int i=0;i<ngals;++i) interp_psf[i] = (*_fitpsf)(_pos[i]);
    }

    col_names[0] = _params.get("shear_id_col");
    col_names[1] = _params.get("shear_x_col");
    col_names[2] = _params.get("shear_y_col");
    col_names[3] = _params.get("shear_sky_col");
    col_names[4] = _params.get("shear_noise_col");
    col_names[5] = _params.get("shear_flags_col");
    col_names[6] = _params.get("shear_ra_col");
    col_names[7] = _params.get("shear_dec_col");
    col_names[8] = _params.get("shear_shear1_col");
    col_names[9] = _params.get("shear_shear2_col");
    col_names[10] = _params.get("shear_nu_col");
    col_names[11] = _params.get("shear_cov00_col");
    col_names[12] = _params.get("shear_cov01_col");
    col_names[13] = _params.get("shear_cov11_col");
    col_names[14] = _params.get("shear_order_col");
    col_names[15] = _params.get("shear_sigma_col");
    col_names[16] = _params.get("shear_coeffs_col");
    if (output_psf) {
        col_names[17] = _params.get("shear_psforder_col");
        col_names[18] = _params.get("shear_psfsigma_col");
        col_names[19] = _params.get("shear_psfcoeffs_col");
    }

    col_fmts[0] = "1J"; // id
    col_fmts[1] = "1D"; // x
    col_fmts[2] = "1D"; // y
    col_fmts[3] = "1D"; // sky
    col_fmts[4] = "1D"; // noise
    col_fmts[5] = "1J"; // flags
    col_fmts[6] = "1D"; // ra
    col_fmts[7] = "1D"; // dec
    col_fmts[8] = "1D"; // shear1
    col_fmts[9] = "1D"; // shear2
    col_fmts[10] = "1D"; // nu
    col_fmts[11] = "1D"; // cov00
    col_fmts[12] = "1D"; // cov01
    col_fmts[13] = "1D"; // cov11
    col_fmts[14] = "1J"; // order
    col_fmts[15] = "1D"; // sigma

    int ncoeff = _shape[0].size();
    dbg<<"ncoeff = "<<ncoeff<<std::endl;
    col_fmts[16] = ConvertibleString(ncoeff) + "D"; // shapelet coeffs
    dbg<<"colfmts[16] = "<<col_fmts[16]<<std::endl;

    if (output_psf) {
        col_fmts[17] = "1J"; // psforder
        col_fmts[18] = "1D"; // psfsigma
        int npsf_coeff = interp_psf[0].size();
        dbg<<"npsf_coeff = "<<npsf_coeff<<std::endl;
        col_fmts[19] = ConvertibleString(npsf_coeff) + "D"; // psf coeffs
        dbg<<"col_fmts[19] = "<<col_fmts[19]<<std::endl;
    }

    col_units[0] = "None";   // id
    col_units[1] = "Pixels"; // x
    col_units[2] = "Pixels"; // y
    col_units[3] = "ADU";    // sky
    col_units[4] = "ADU^2";  // noise
    col_units[5] = "None";   // flags
    col_units[6] = "Deg";    // ra
    col_units[7] = "Deg";    // dec
    col_units[8] = "None";   // shear1
    col_units[9] = "None";   // shear2
    col_units[10] = "None";  // nu
    col_units[11] = "None";  // cov00
    col_units[12] = "None";  // cov01
    col_units[13] = "None";  // cov11
    col_units[14] = "None";  // order
    col_units[15] = "Arcsec";// sigma
    col_units[16] = "None";  // coeffs
    if (output_psf) {
        col_units[17] = "None";  // order
        col_units[18] = "Arcsec";// sigma
        col_units[19] = "None";  // coeffs
    }

    dbg<<"Before Create table"<<std::endl;
    CCfits::Table* table;
    table = fits.addTable("shearcat",ngals,col_names,col_fmts,col_units);

    // Header keywords
#ifdef USE_TMV
    std::string tmvVers = tmv::TMV_Version();
#else
    std::string tmvVers = "Eigen";
#endif
    std::string wlVers = GetWlVersion();

    table->addKey("tmvvers", tmvVers, "version of TMV code");
    table->addKey("wlvers", wlVers, "version of weak lensing code");

    std::string str;
    double dbl;
    int intgr;

    // if wlserun= is sent we'll put it in the header.  This allows us to 
    // associate some more, possibly complicated, metadata with this file
    if ( _params.keyExists("wlserun") ) {
        WriteParamToTable(_params, table, "wlserun", str);
    }

    WriteParamToTable(_params, table, "noise_method", str);
    WriteParamToTable(_params, table, "dist_method", str);

    WriteParamToTable(_params, table, "shear_aperture", dbl);
    WriteParamToTable(_params, table, "shear_max_aperture", dbl);
    WriteParamToTable(_params, table, "shear_gal_order", intgr);
    WriteParamToTable(_params, table, "shear_gal_order2", intgr);
    WriteParamToTable(_params, table, "shear_min_gal_size", dbl);
    WriteParamToTable(_params, table, "shear_f_psf", dbl);

    // data
    // make vector copies for writing
    std::vector<double> x(ngals);
    std::vector<double> y(ngals);
    std::vector<double> ra(ngals);
    std::vector<double> dec(ngals);
    std::vector<double> shear1(ngals);
    std::vector<double> shear2(ngals);
    std::vector<double> cov00(ngals);
    std::vector<double> cov01(ngals);
    std::vector<double> cov11(ngals);

    for(int i=0;i<ngals;++i) {
        x[i] = _pos[i].getX();
        y[i] = _pos[i].getY();
        // internally we use arcseconds
        ra[i] = _skypos[i].getX()/3600.;
        dec[i] = _skypos[i].getY()/3600.;
        shear1[i] = real(_shear[i]);
        shear2[i] = imag(_shear[i]);
        cov00[i] = _cov[i](0,0);
        cov01[i] = _cov[i](0,1);
        cov11[i] = _cov[i](1,1);
    }

    int startRow=1;

    table->column(col_names[0]).write(_id,startRow);
    table->column(col_names[1]).write(x,startRow);
    table->column(col_names[2]).write(y,startRow);
    table->column(col_names[3]).write(_sky,startRow);
    table->column(col_names[4]).write(_noise,startRow);
    table->column(col_names[5]).write(_flags,startRow);
    table->column(col_names[6]).write(ra,startRow);
    table->column(col_names[7]).write(dec,startRow);
    table->column(col_names[8]).write(shear1,startRow);
    table->column(col_names[9]).write(shear2,startRow);
    table->column(col_names[10]).write(_nu,startRow);
    table->column(col_names[11]).write(cov00,startRow);
    table->column(col_names[12]).write(cov01,startRow);
    table->column(col_names[13]).write(cov11,startRow);
    table->column(col_names[14]).write(_meas_galorder,startRow);

    for (int i=0; i<ngals; ++i) {
        int row = i+1;

        // MJ: meas_galorder keeps track of the order of the shapelet that
        // was actually measured.  It is <= the full order of the shapelet
        // vector, but the higher order terms are set to zero.
        // Since we can't have different numbers of columns for each row, 
        // we have to write out the full shapelet order for all rows
        // including the zeros.
        //
        //long border = _shape[i].getOrder();
        //table->column(col_names[14]).write(&border,1,row);
        
        double bsigma = _shape[i].getSigma();
        table->column(col_names[15]).write(&bsigma,1,row);

        // There is a bug in the CCfits Column::write function that the 
        // first argument has to be S* rather than const S*, even though
        // it is only read from. 
        // Hence the const_cast.
        double* ptr = const_cast<double*>(TMV_cptr(_shape[i].vec()));
        table->column(col_names[16]).write(ptr, ncoeff, 1, row);
    }

    if (output_psf) {
        for (int i=0; i<ngals; ++i) {
            int row = i+1;
            long psfOrder = interp_psf[i].getOrder();
            double psfSigma = interp_psf[i].getSigma();

            table->column(col_names[17]).write(&psfOrder,1,row);
            table->column(col_names[18]).write(&psfSigma,1,row);
            double* ptr = TMV_ptr(interp_psf[i].vec());
            Assert(interp_psf[i].size() == interp_psf[0].size());
            int npsf_coeff = interp_psf[i].size();
            table->column(col_names[19]).write(ptr, npsf_coeff, 1, row);
        }
    }
}

void ShearCatalog::writeAscii(std::string file, std::string delim) const
{
    Assert(int(_id.size()) == size());
    Assert(int(_pos.size()) == size());
    Assert(int(_sky.size()) == size());
    Assert(int(_noise.size()) == size());
    Assert(int(_flags.size()) == size());
    Assert(int(_skypos.size()) == size());
    Assert(int(_shear.size()) == size());
    Assert(int(_nu.size()) == size());
    Assert(int(_cov.size()) == size());
    Assert(int(_meas_galorder.size()) == size());
    Assert(int(_shape.size()) == size());

    std::ofstream fout(file.c_str());
    if (!fout) {
        throw WriteException("Error opening shear file"+file);
    }

    bool output_psf = _params.read("shear_output_psf",false);

    Form sci8; sci8.sci().trail(0).prec(8);
    Form fix3; fix3.fix().trail(0).prec(3);
    Form fix8; fix8.fix().trail(0).prec(8);
    Form fix6; fix6.fix().trail(0).prec(6);

    const int ngals = size();
    for(int i=0;i<ngals;++i) {
        fout
            << _id[i] << delim
            << fix3(_pos[i].getX()) << delim
            << fix3(_pos[i].getY()) << delim
            << fix3(_sky[i]) << delim
            << sci8(_noise[i]) << delim
            << _flags[i] << delim
            << fix8(_skypos[i].getX()/3600.) << delim 
            << fix8(_skypos[i].getY()/3600.) << delim
            << sci8(real(_shear[i])) << delim
            << sci8(imag(_shear[i])) << delim
            << fix3(_nu[i]) << delim
            << sci8(_cov[i](0,0)) << delim
            << sci8(_cov[i](0,1)) << delim
            << sci8(_cov[i](1,1)) << delim
            //<< _shape[i].getOrder() << delim
            << _meas_galorder[i] << delim
            << fix6(_shape[i].getSigma());
        const int ncoeff = _shape[i].size();
        for(int j=0;j<ncoeff;++j)
            fout << delim << sci8(_shape[i](j));
        if (output_psf) {
            Assert(_fitpsf);
            BVec interp_psf(_fitpsf->getPsfOrder(),_fitpsf->getSigma());
            interp_psf = (*_fitpsf)(_pos[i]);
            fout
                << delim << interp_psf.getOrder()
                << delim << fix6(interp_psf.getSigma());
            const int npsf_coeff = interp_psf.size();
            for(int j=0;j<npsf_coeff;++j)
                fout << delim << sci8(interp_psf(j));
        }
        fout << std::endl;
    }
}

void ShearCatalog::write() const
{
    std::vector<std::string> file_list = MakeMultiName(_params, "shear");  

    const int nfiles = file_list.size();
    for(int i=0; i<nfiles; ++i) {
        const std::string& file = file_list[i];
        dbg<<"Writing shear catalog to file: "<<file<<std::endl;

        bool isFitsIo = false;
        if (_params.keyExists("shear_io")) {
            std::vector<std::string> io_list = _params["shear_io"];
            Assert(io_list.size() == file_list.size());
            isFitsIo = (io_list[i] == "FITS");
        } else if (file.find("fits") != std::string::npos) {
            isFitsIo = true;
        }

        try {
            if (isFitsIo) {
                writeFits(file);
            } else {
                std::string delim = "  ";
                if (_params.keyExists("shear_delim")) {
                    std::vector<std::string> delim_list = _params["shear_delim"];
                    Assert(delim_list.size() == file_list.size());
                    delim = delim_list[i];
                } else if (file.find("csv") != std::string::npos) {
                    delim = ",";
                }
                writeAscii(file,delim);
            }
        } catch (CCfits::FitsException& e) {
            xdbg<<"Caught FitsException: \n"<<e.message()<<std::endl;
            throw WriteException(
                "Error writing to "+file+" -- caught error\n" + e.message());
        } catch (std::exception& e) {
            xdbg<<"Caught std::exception: \n"<<e.what()<<std::endl;
            throw WriteException(
                "Error writing to "+file+" -- caught error\n" + e.what());
        } catch (...) {
            xdbg<<"Caught unknown exception"<<std::endl;
            throw WriteException(
                "Error writing to "+file+" -- caught unknown error");
        }
    }
    dbg<<"Done Write ShearCatalog\n";
}

void ShearCatalog::readFits(std::string file)
{
    int hdu = GetHdu(_params,"shear",file,2);

    dbg<<"Opening ShearCatalog file "<<file<<" at hdu "<<hdu<<std::endl;
    CCfits::FITS fits(file, CCfits::Read);
    if (hdu > 1) fits.read(hdu-1);

    CCfits::ExtHDU& table=fits.extension(hdu-1);

    long nRows=table.rows();

    dbg<<"  nrows = "<<nRows<<std::endl;
    if (nRows <= 0) {
        throw ReadException(
            "ShearCatalog found to have 0 rows.  Must have > 0 rows.");
    }

    std::string id_col=_params.get("shear_id_col");
    std::string x_col=_params.get("shear_x_col");
    std::string y_col=_params.get("shear_y_col");
    std::string sky_col=_params.get("shear_sky_col");
    std::string noise_col=_params.get("shear_noise_col");
    std::string flags_col=_params.get("shear_flags_col");
    std::string ra_col=_params.get("shear_ra_col");
    std::string dec_col=_params.get("shear_dec_col");
    std::string shear1_col=_params.get("shear_shear1_col");
    std::string shear2_col=_params.get("shear_shear2_col");
    std::string nu_col=_params.get("shear_nu_col");
    std::string cov00_col=_params.get("shear_cov00_col");
    std::string cov01_col=_params.get("shear_cov01_col");
    std::string cov11_col=_params.get("shear_cov11_col");
    std::string order_col=_params.get("shear_order_col");
    std::string sigma_col=_params.get("shear_sigma_col");
    std::string coeffs_col=_params.get("shear_coeffs_col");

    // MJ: This is not really optimal.
    // We should get this value from the fits header, since it's stored there.  
    // What is the CCfits command to do this?
    int full_order = _params.get("shear_gal_order");

    long start=1;
    long end=nRows;

    dbg<<"Reading columns"<<std::endl;
    dbg<<"  "<<id_col<<std::endl;
    table.column(id_col).read(_id, start, end);

    dbg<<"  "<<x_col<<"  "<<y_col<<std::endl;
    _pos.resize(nRows);
    std::vector<double> x;
    std::vector<double> y;
    table.column(x_col).read(x, start, end);
    table.column(y_col).read(y, start, end);
    for(long i=0;i<nRows;++i) _pos[i] = Position(x[i],y[i]);

    dbg<<"  "<<sky_col<<std::endl;
    table.column(sky_col).read(_sky, start, end);

    dbg<<"  "<<noise_col<<std::endl;
    table.column(noise_col).read(_noise, start, end);

    dbg<<"  "<<flags_col<<std::endl;
    table.column(flags_col).read(_flags, start, end);

    dbg<<"  "<<ra_col<<"  "<<dec_col<<std::endl;
    _skypos.resize(nRows);
    std::vector<double> ra;
    std::vector<double> dec;
    table.column(ra_col).read(ra, start, end);
    table.column(dec_col).read(dec, start, end);
    for(long i=0;i<nRows;++i) {
        _skypos[i] = Position(ra[i]*3600.,dec[i]*3600.);
    }

    dbg<<"  "<<shear1_col<<"  "<<shear2_col<<std::endl;
    _shear.resize(nRows);
    std::vector<double> shear1;
    std::vector<double> shear2;
    table.column(shear1_col).read(shear1, start, end);
    table.column(shear2_col).read(shear2, start, end);
    for(long i=0;i<nRows;++i) {
        _shear[i] = std::complex<double>(shear1[i],shear2[i]);
    }

    dbg<<"  "<<nu_col<<std::endl;
    table.column(nu_col).read(_nu, start, end);

    dbg<<"  "<<cov00_col<<"  "<<cov01_col<<"  "<<cov11_col<<std::endl;
    _cov.resize(nRows);
    std::vector<double> cov00;
    std::vector<double> cov01;
    std::vector<double> cov11;
    table.column(cov00_col).read(cov00, start, end);
    table.column(cov01_col).read(cov01, start, end);
    table.column(cov11_col).read(cov11, start, end);
    for(long i=0;i<nRows;++i) {
        _cov[i] << cov00[i], cov01[i], cov01[i], cov11[i];
    }

    dbg<<"  "<<sigma_col<<"  "<<order_col<<std::endl;
    std::vector<double> sigma(nRows); // temporary
    table.column(sigma_col).read(sigma, start, end);
    table.column(order_col).read(_meas_galorder, start, end);

    dbg<<"  "<<coeffs_col<<std::endl;
    _shape.reserve(nRows);
    for (int i=0; i<nRows; ++i) {
        int row=i+1;

        _shape.push_back(BVec(full_order,sigma[i]));
        int ncoeff=_shape[i].size();

        std::valarray<double> coeffs(ncoeff);
        table.column(coeffs_col).read(coeffs, row);
        for(int j=0;j<ncoeff;++j) _shape[i](j) = coeffs[j];
        xxdbg<<"shape => "<<_shape[i].vec()<<std::endl;
    }
    dbg<<"Done ReadFits\n";
}

void ShearCatalog::readAscii(std::string file, std::string delim)
{
    std::ifstream fin(file.c_str());
    if (!fin) {
        throw ReadException("Error opening stars file"+file);
    }

    int full_order = _params.get("shear_gal_order");

    _id.clear(); _pos.clear(); _sky.clear(); _noise.clear(); _flags.clear();
    _skypos.clear(); _shear.clear(); _nu.clear(); _cov.clear();
    _meas_galorder.clear(); _shape.clear();

    if (delim == "  ") {
        ConvertibleString flag;
        long id1,border;
        double x,y,sky1,noise1,ra,dec,s1,s2,nu1,c00,c01,c11,bsigma;
        while ( fin >> id1 >> x >> y >> sky1 >> noise1 >> 
                flag >> ra >> dec >> s1 >> s2 >> nu1 >>
                c00 >> c01 >> c11 >> border >> bsigma ) {
            _id.push_back(id1);
            _pos.push_back(Position(x,y));
            _sky.push_back(sky1);
            _noise.push_back(noise1);
            _flags.push_back(flag);
            _skypos.push_back(Position(ra*3600.,dec*3600.));
            _shear.push_back(std::complex<double>(s1,s2));
            _nu.push_back(nu1);
            _cov.push_back(DSmallMatrix22());
            _cov.back() << c00, c01, c01, c11;
            _meas_galorder.push_back(border);
            _shape.push_back(BVec(full_order,bsigma));
            const int ncoeff = _shape.back().size();
            for(int j=0; j<ncoeff; ++j)
                fin >> _shape.back()(j);
        }
    } else {
        if (delim.size() > 1) {
            // getline only works with single character delimiters.
            // Since I don't really expect a multicharacter delimiter to
            // be used ever, I'm just going to throw an exception here 
            // if we do need it, and I can write the workaround then.
            throw ParameterException(
                "ReadAscii delimiter must be a single character");
        }
        char d = delim[0];
        long border;
        double x,y,ra,dec,s1,s2,bsigma,c00,c01,c11;
        ConvertibleString temp;
        while (getline(fin,temp,d)) {
            _id.push_back(temp);
            getline(fin,temp,d); x = temp;
            getline(fin,temp,d); y = temp;
            _pos.push_back(Position(x,y));
            getline(fin,temp,d); _sky.push_back(temp);
            getline(fin,temp,d); _noise.push_back(temp);
            getline(fin,temp,d); _flags.push_back(temp);
            getline(fin,temp,d); ra = temp;
            getline(fin,temp,d); dec = temp;
            _skypos.push_back(Position(ra*3600.,dec*3600.));
            getline(fin,temp,d); s1 = temp;
            getline(fin,temp,d); s2 = temp;
            _shear.push_back(std::complex<double>(s1,s2));
            getline(fin,temp,d); _nu.push_back(temp);
            getline(fin,temp,d); c00 = temp;
            getline(fin,temp,d); c01 = temp;
            getline(fin,temp,d); c11 = temp;
            _cov.push_back(DSmallMatrix22());
            _cov.back() << c00, c01, c01, c11;
            getline(fin,temp,d); border = temp;
            getline(fin,temp,d); bsigma = temp;
            _meas_galorder.push_back(border);
            _shape.push_back(BVec(full_order,bsigma));
            const int ncoeff = _shape.back().size();
            for(int j=0;j<ncoeff-1;++j) {
                getline(fin,temp,d); _shape.back()(j) = temp;
            }
            // Last one doesn't have a following delimiter
            getline(fin,temp); _shape.back()(ncoeff-1) = temp;
        }
    }
}

void ShearCatalog::read()
{
    std::string file = MakeName(_params,"shear",false,true);
    // false,true = input_prefix=false, mustexist=true.
    // It is an input here, but it is in the output_prefix directory.
    dbg<< "Reading Shear cat from file: " << file << std::endl;

    bool isFitsIo = false;
    if (_params.keyExists("shear_io")) 
        isFitsIo = (_params["shear_io"] == "FITS");
    else if (file.find("fits") != std::string::npos) 
        isFitsIo = true;

    if (!DoesFileExist(file)) {
        throw FileNotFoundException(file);
    }
    try {
        if (isFitsIo) {
            readFits(file);
        } else {
            std::string delim = "  ";
            if (_params.keyExists("shear_delim")) {
                delim = _params["shear_delim"];
            } else if (file.find("csv") != std::string::npos) {
                delim = ",";
            }
            readAscii(file,delim);
        }
    } catch (CCfits::FitsException& e) {
        xdbg<<"Caught FitsException: \n"<<e.message()<<std::endl;
        throw ReadException(
            "Error reading from "+file+" -- caught error\n" + e.message());
    } catch (std::exception& e) {
        xdbg<<"Caught std::exception: \n"<<e.what()<<std::endl;
        throw ReadException(
            "Error reading from "+file+" -- caught error\n" + e.what());
    } catch (...) {
        xdbg<<"Caught unknown exception"<<std::endl;
        throw ReadException(
            "Error reading from "+file+" -- caught unknown error");
    }

    // Update the bounds:
    const int ngals = size();
    Assert(int(_pos.size()) == size());
    Assert(int(_skypos.size()) == size());
    for(int i=0;i<ngals;++i) _bounds += _pos[i];

    // Calculate the skybounds using only the good objects:
    Bounds skybounds2; // In degrees, rather than arcsec
    for(int i=0;i<ngals;++i) if (!_flags[i]) {
        _skybounds += _skypos[i];
        Position temp = _skypos[i];
        temp /= 3600.;
        skybounds2 += temp;
    }
    dbg<<"skybounds = "<<_skybounds<<std::endl;
    dbg<<"in degrees: "<<skybounds2<<std::endl;

    dbg<<"Done Read ShearCatalog\n";
}



