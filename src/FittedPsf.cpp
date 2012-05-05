
#include <valarray>
#include <fstream>
#include <CCfits/CCfits>

#include "dbg.h"
#include "FittedPsf.h"
#include "Function2D.h"
#include "Legendre2D.h"
#include "Name.h"
#include "WlVersion.h"
#include "WriteParam.h"

static DVector definePXY(
    int order, double x, double xmin, double xmax)
{
    DVector temp(order+1);
    double newx = (2.*x-xmin-xmax)/(xmax-xmin);
    temp[0] = 1.;
    if(order>0) temp[1] = newx;
    for(int i=2;i<=order;++i) {
        temp[i] = ((2.*i-1.)*newx*temp[i-1] - (i-1.)*temp[i-2])/i;
    }
    return temp;
}

#ifdef USE_TMV
static void setPRow(
    int fitorder, Position pos, const Bounds& bounds, DVectorView prow)
#else
static void setPRow(
    int fitorder, Position pos, const Bounds& bounds, DVector& prow)
#endif
{
    Assert(int(prow.size()) == (fitorder+1)*(fitorder+2)/2);
    DVector px = 
        definePXY(fitorder,pos.getX(),bounds.getXMin(),bounds.getXMax());
    DVector py = 
        definePXY(fitorder,pos.getY(),bounds.getYMin(),bounds.getYMax());
    int pq = 0;
    for(int n=0;n<=fitorder;++n) {
        for(int p=n,q=n-p;q<=n;--p,++q) {
            Assert(pq < int(prow.size()));
            prow(pq) = px[p]*py[q];
            ++pq;
        }
    }
    Assert(pq == int(prow.size()));
}

FittedPsf::FittedPsf(
    PsfCatalog& psfcat, const ConfigFile& params, PsfLog& log) :
    _params(params), _psforder(_params.read<int>("psf_order")),
    _fitorder(_params.read<int>("fitpsf_order")),
    _fitsize((_fitorder+1)*(_fitorder+2)/2)
{
    xdbg<<"FittedPSF constructor\n";
    // Do a polynomial fit of the psf shapelet vectors

    Assert(psfcat.size() > 0);

    // This used to be set to the sigma of psfcat.getPsf(0).
    // This doesn't work if the first object failed. 
    // So now I set it to -1, and then update it when we get to the first
    // object without an error flag (usually n = 0).
    _sigma = -1.;
    const int nstars = psfcat.size();
    for(int n=0;n<nstars;++n) if (!psfcat.getFlags(n)) {
        _sigma = psfcat.getPsf(n).getSigma();
        break;
    }

    std::vector<long>& flags = psfcat.getFlagsList();
    //split_psf_stars(flags);

    calculate(
        psfcat.getPosList(),
        psfcat.getPsfList(),
        psfcat.getNuList(),
        flags,
        log);
}

void FittedPsf::calculate(
    const std::vector<Position>& pos,
    const std::vector<BVec>& psf,
    const std::vector<double>& nu,
    std::vector<long>& flags, 
    PsfLog& log)
{
    const int nstars = pos.size();
    const int psfsize = (_psforder+1)*(_psforder+2)/2;
    const double nsigma_clip = _params.read("fitpsf_nsigma_outlier",3);

    // This is an empirical fit to the chisq level that corresponds to 
    // a 3-sigma outlier for more than 1 dimension.
    // I calculated the values up to n=100 and for n>30, they form a pretty
    // good approximation to a straight line.
    // This is almost certainly wrong for nsigma_clip != 3, so if we start 
    // choosing other values for nsigma_clip, it might be worth doing this 
    // right.
    // That means calculating the 1-d critical value for the given nSigma.
    // e.g. nSigma = 3 -> alpha = P(chisq > 9) = 0.0027.
    // Then calculate the critical value of chisq for that alpha with 
    // the full degrees of freedom = psfsize
    const double chisq_level = 0.14*psfsize + 2.13;

    const double outlier_thresh = nsigma_clip * nsigma_clip * chisq_level;
    dbg<<"outlier_thresh = "<<outlier_thresh<<std::endl;

    int ngood_psf, noutliers, dof;
    double chisq;
    do {

        // Calculate the average psf vector
        _avepsf.reset(new DVector(psfsize));
        Assert(psfsize == int(_avepsf->size()));
        _avepsf->setZero();
        ngood_psf = 0;
        for(int n=0;n<nstars;++n) if ( flags[n]==0 ) {
            Assert(psf[n].getSigma() == _sigma);
            *_avepsf += psf[n].vec();
            ++ngood_psf;
        }
        xdbg<<"ngood_psf = "<<ngood_psf<<std::endl;
        if (ngood_psf == 0) {
            dbg<<"ngoodpsf = 0 in FittedPsf::calculate\n";
            throw ProcessingException("No good stars found for interpolation.");
        }
        *_avepsf /= double(ngood_psf);
        xdbg<<"_avepsf = "<<*_avepsf<<std::endl;

        // Rotate the vectors into their eigen directions.
        // The matrix V is stored to let us get back to the original basis.
        DMatrix mM(ngood_psf,psfsize);
        DDiagMatrix inverseSigma(ngood_psf);
        int i=0;
        for(int n=0;n<nstars;++n) if ( flags[n]==0 ) {
            Assert(int(psf[n].size()) == psfsize);
            Assert(i < ngood_psf);
            mM.row(i) = psf[n].vec() - *_avepsf;
            inverseSigma(i) = nu[n];
            _bounds += pos[n];
            ++i;
        }
        Assert(i == ngood_psf);
        xdbg<<"bounds = "<<_bounds<<std::endl;
        mM = inverseSigma EIGEN_asDiag() * mM;

        int npca_tot = std::min(ngood_psf,psfsize);
        xdbg<<"npca_tot = "<<npca_tot<<std::endl;
        DDiagMatrix mS(npca_tot);
#ifdef USE_TMV
        DMatrixView mU = mM.colRange(0,npca_tot);
        _mV.reset(new tmv::Matrix<double,tmv::RowMajor>(npca_tot,psfsize));
        if (ngood_psf > psfsize) {
            SV_Decompose(mU.view(),mS.view(),_mV->view(),true);
        } else {
            *_mV = mM;
            SV_Decompose(_mV->transpose(),mS.view(),mU.transpose());
        }
        xdbg<<"In FittedPSF: SVD S = "<<mS.diag()<<std::endl;
        //xdbg<<"U = "<<mU<<std::endl;
        //xdbg<<"V = "<<*_mV<<std::endl;
#else
        DMatrix mU(mM.TMV_colsize(),npca_tot);
        _mV_transpose.reset(new DMatrix(psfsize,npca_tot));
        if (ngood_psf > psfsize) {
            Eigen::SVD<DMatrix> svd = TMV_colRange(mM,0,npca_tot).svd();
            mU = svd.matrixU();
            mS = svd.singularValues();
            *_mV_transpose = svd.matrixV();
        } else {
            Eigen::SVD<Eigen::Transpose<DMatrix>::PlainMatrixType > svd = mM.transpose().svd();
            mU = svd.matrixV();
            mS = svd.singularValues();
            *_mV_transpose = svd.matrixU();
        }
        xdbg<<"In FittedPSF: SVD S = "<<EIGEN_Transpose(mS)<<std::endl;
#endif
        if (_params.keyExists("fitpsf_npca")) {
            _npca = _params["fitpsf_npca"];
            dbg<<"npca = "<<_npca<<" from parameter file\n";
        } else {
            double thresh = mS(0);
            if (_params.keyExists("fitpsf_pca_thresh")) 
                thresh *= double(_params["fitpsf_pca_thresh"]);
            else thresh *= std::numeric_limits<double>::epsilon();
            dbg<<"thresh = "<<thresh<<std::endl;
            for(_npca=1;_npca<int(mS.size());++_npca) {
                if (mS(_npca) < thresh) break;
            }
            dbg<<"npca = "<<_npca<<std::endl;
        }
#ifdef USE_TMV
        mU.colRange(0,_npca) *= mS.subDiagMatrix(0,_npca);
#else
        TMV_colRange(mU,0,_npca) *= mS.TMV_subVector(0,_npca).asDiagonal();
#endif
        xdbg<<"After U *= S\n";
        // U S = M(orig) * Vt

        while (ngood_psf <= _fitsize && _fitsize > 1) {
            --_fitorder;
            _fitsize = (_fitorder+1)*(_fitorder+2)/2;
            dbg<<"Too few good stars... reducing order of fit to "<<
                _fitorder<<std::endl;
        }
        DMatrix mP(ngood_psf,_fitsize);
        mP.setZero();
        i=0;
        for(int n=0;n<nstars;++n) if ( flags[n]==0 ) {
            xdbg<<"n = "<<n<<" / "<<nstars<<std::endl;
#ifdef USE_TMV
            setPRow(_fitorder,pos[n],_bounds,mP.row(i));
#else
            DVector mProwi(mP.TMV_rowsize());
            setPRow(_fitorder,pos[n],_bounds,mProwi);
            mP.row(i) = mProwi.transpose();
#endif
            ++i;
        }
        Assert(i == ngood_psf);
        mP = inverseSigma EIGEN_asDiag() * mP;
        xdbg<<"after mP = sigma * mP\n";

#ifdef USE_TMV
        _f.reset(new DMatrix(TMV_colRange(mU,0,_npca)/mP));
#else
        _f.reset(new DMatrix(_fitsize,_npca));
        mP.qr().solve(TMV_colRange(mU,0,_npca),&(*_f));
#endif
        xdbg<<"Done making FittedPSF\n";

        //
        // Remove outliers from the fit using the empirical covariance matrix
        // of the data with respect to the fitted values.
        //
        xdbg<<"Checking for outliers:\n";

        // Calculate the covariance matrix
        DMatrix cov(psfsize,psfsize);
        cov.setZero();
        for(int n=0;n<nstars;++n) if ( flags[n]==0 ) {
            const BVec& data = psf[n];
            DVector fit(psfsize);
            interpolateVector(pos[n],TMV_vview(fit));
            DVector diff = data.vec() - fit;
#ifdef USE_TMV
            cov += diff ^ diff;
#else
            cov += diff * diff.transpose();
#endif
        }

        chisq = (mP * *_f - TMV_colRange(mU,0,_npca)).TMV_normSq();
        dbg<<"chisq calculation #1 = "<<chisq<<std::endl;
        dof = ngood_psf - _fitsize;

        if (dof > 0) { 
            cov /= double(dof);
        }
#ifdef USE_TMV
        cov.divideUsing(tmv::SV);
        cov.saveDiv();
        cov.setDiv();
        dbg<<"cov S = "<<cov.svd().getS().diag()<<std::endl;
#else
        Eigen::SVD<DMatrix> cov_svd = cov.svd();
        dbg<<"cov S = "<<EIGEN_Transpose(cov_svd.singularValues())<<std::endl;

#endif

        // Clip out 3 sigma outliers:
        noutliers = 0;
        chisq = 0;
        for(int n=0;n<nstars;++n) if ( flags[n]==0 ) {
            const BVec& data = psf[n];
            DVector fit(psfsize);
            interpolateVector(pos[n],TMV_vview(fit));
            DVector diff = data.vec() - fit;
#ifdef USE_TMV
            double dev = diff * cov.inverse() * diff;
#else
            DVector temp(psfsize);
            cov_svd.solve(diff,&temp);
            double dev = (diff.transpose() * temp)(0,0);
#endif
            chisq += dev;
            if (dev > outlier_thresh) {
                xdbg<<"n = "<<n<<" is an outlier.\n";
                xdbg<<"data = "<<data.vec()<<std::endl;
                xdbg<<"fit = "<<fit<<std::endl;
                xdbg<<"diff = "<<diff<<std::endl;
#ifdef USE_TMV
                xdbg<<"diff/cov = "<<diff/cov<<std::endl;
#else
                xdbg<<"diff/cov = "<<temp<<std::endl;
#endif
                xdbg<<"dev = "<<dev<<std::endl;
                ++noutliers;
                flags[n] |= PSF_INTERP_OUTLIER;
            }
        }
        dbg<<"ngoodpsf = "<<ngood_psf<<std::endl;
        dbg<<"noutliers = "<<noutliers<<std::endl;
        dbg<<"chisq calculation #2 = "<<chisq<<std::endl;

    } while (noutliers > 0);

    log._noutliers = noutliers;
    log._nfit = ngood_psf;
    log._npc = _npca;
    log._chisq_fit = chisq;
    log._dof_fit = dof;

}

FittedPsf::FittedPsf(const ConfigFile& params) : 
    _params(params), _psforder(_params.read<int>("psf_order")),
    _fitorder(_params.read<int>("fitpsf_order")),
    _fitsize((_fitorder+1)*(_fitorder+2)/2)
{ }

void FittedPsf::writeAscii(std::string file) const
{
    std::ofstream fout(file.c_str());
    fout << _psforder <<"  "<< _sigma <<"  ";
    fout << _fitorder <<"  "<< _npca << std::endl;
    fout << _bounds << std::endl;
    fout << *_avepsf << std::endl;
#ifdef USE_TMV
    fout << _mV->rowRange(0,_npca) <<std::endl;
#else
    fout << TMV_colRange(*_mV_transpose,0,_npca).transpose() <<std::endl;
#endif
    fout << *_f << std::endl;
}

void FittedPsf::readAscii(std::string file)
{
    std::ifstream fin(file.c_str());
    if (!fin) {
        throw ReadException("Error opening fitpsf file"+file);
    }

    xdbg<<"Reading FittedPSF:\n";
    fin >> _psforder >> _sigma >> _fitorder >> _npca >> _bounds;
    xdbg<<"psforder = "<<_psforder<<", sigma_psf = "<<_sigma<<std::endl;
    xdbg<<"fitorder = "<<_fitorder<<", npca = "<<_npca<<std::endl;
    xdbg<<"bounds = "<<_bounds<<std::endl;
    _fitsize = (_fitorder+1)*(_fitorder+2)/2;
    const int psfsize = (_psforder+1)*(_psforder+2)/2;
    _avepsf.reset(new DVector(psfsize));
    _f.reset(new DMatrix(_fitsize,_npca));
#ifdef USE_TMV
    fin >> *_avepsf;
    _mV.reset(new tmv::Matrix<double,tmv::RowMajor>(_npca,psfsize));
    fin >> *_mV;
    fin >> *_f;
#else
    for(int i=0;i<_avepsf->size();++i) fin >> (*_avepsf)(i);
    _mV_transpose.reset(new DMatrix(psfsize,_npca));
    for(int i=0;i<_npca;++i) for(int j=0;j<psfsize;++j)
        fin >> (*_mV_transpose)(i,j);
    for(int i=0;i<_fitsize;++i) for(int j=0;j<_npca;++j)
        fin >> (*_f)(i,j);
#endif
    if (!fin) {
        throw ReadException("Error reading fitpsf file"+file);
    }
}

void FittedPsf::write() const
{
    std::vector<std::string> file_list = MakeMultiName(_params, "fitpsf");  

    const int nfiles = file_list.size();
    for(int i=0; i<nfiles; ++i) {
        const std::string& file = file_list[i];
        dbg<<"Writing fitted psf to file: "<<file<<std::endl;

        bool is_fits_io = false;
        if (_params.keyExists("fitpsf_io")) {
            std::vector<std::string> ios = _params["fitpsf_io"];
            Assert(ios.size() == file_list.size());
            is_fits_io = (ios[i] == "FITS");
        } else if (file.find("fits") != std::string::npos) {
            is_fits_io = true;
        }

        try {
            if (is_fits_io) {
                writeFits(file);
            } else {
                writeAscii(file);
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
            xdbg<<"Caught unknown exception: "<<std::endl;
            throw WriteException(
                "Error writing to "+file+" -- caught unknown error");
        }
    }
    dbg<<"Done Write FittedPSF\n";
}

void FittedPsf::read()
{
    std::string file = MakeName(_params,"fitpsf",false,true);
    read(file);
}

void FittedPsf::read(std::string file)
{
    // false,true = input_prefix=false, mustexist=true.
    // It is an input here, but it is in the output_prefix directory.
    dbg<< "Reading FittedPSF from file: " << file << std::endl;

    bool is_fits_io = false;
    if (_params.keyExists("fitpsf_io")) 
        is_fits_io = (_params["fitpsf_io"] == "FITS");
    else if (file.find("fits") != std::string::npos) 
        is_fits_io = true;

    if (!DoesFileExist(file)) {
        throw FileNotFoundException(file);
    }
    try {
        if (is_fits_io) {
            readFits(file);
        } else {
            readAscii(file);
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
        xdbg<<"Caught unknown exception: "<<std::endl;
        throw ReadException(
            "Error reading from "+file+" -- caught unknown error");
    }
    dbg<<"Done Read FittedPSF\n";
}

void FittedPsf::writeFits(std::string file) const
{
    dbg<<"Start WriteFits"<<std::endl;

    // ! means overwrite existing file
    CCfits::FITS fits("!"+file, CCfits::Write);

    dbg<<"Made fits"<<std::endl;

    // Note the actual coeffs may be less than this the way Mike does things
    // but we will fill the extra with zeros
    int nshapelet_coeff = (_psforder+1)*(_psforder+2)/2;
    int nrot_matrix = _npca*nshapelet_coeff;

    int nfit_coeff = (_fitorder+1)*(_fitorder+2)/2;
    int ninterp_matrix = _npca*nfit_coeff;

    const int nfields=11;

    std::vector<string> col_names(nfields);
    std::vector<string> col_fmts(nfields);
    std::vector<string> col_units(nfields);

    col_names[0] = _params["fitpsf_psf_order_col"];
    col_names[1] = _params["fitpsf_sigma_col"];
    col_names[2] = _params["fitpsf_fit_order_col"];
    col_names[3] = _params["fitpsf_npca_col"];

    col_names[4] = _params["fitpsf_xmin_col"];
    col_names[5] = _params["fitpsf_xmax_col"];
    col_names[6] = _params["fitpsf_ymin_col"];
    col_names[7] = _params["fitpsf_ymax_col"];

    col_names[8] = _params["fitpsf_ave_psf_col"];
    col_names[9] = _params["fitpsf_rot_matrix_col"];
    col_names[10] = _params["fitpsf_interp_matrix_col"];

    col_fmts[0] = "1J"; // psf order
    col_fmts[1] = "1D"; // sigma
    col_fmts[2] = "1J"; // fit order
    col_fmts[3] = "1J"; // npca
    col_fmts[4] = "1E"; // xmin
    col_fmts[5] = "1E"; // xmax
    col_fmts[6] = "1E"; // ymin
    col_fmts[7] = "1E"; // ymax

    std::stringstream fmt;
    fmt << nshapelet_coeff << "D";
    col_fmts[8] = fmt.str(); // average psf shapelets
    fmt.str("");
    fmt << nrot_matrix << "D";
    col_fmts[9] = fmt.str(); // rotation matrix
    fmt.str("");
    fmt << ninterp_matrix << "D";
    col_fmts[10] = fmt.str(); // interp matrix

    col_units[0] = "None";     // psf order
    col_units[1] = "arcsec";   // sigma
    col_units[2] = "None";     // fit order
    col_units[3] = "None";     // npca
    col_units[4] = "pixels";   // xmin
    col_units[5] = "pixels";   // xmax
    col_units[6] = "pixels";   // ymin
    col_units[7] = "pixels";   // ymax
    col_units[8] = "None";     // avg psf shapelets
    col_units[9] = "None";     // rot_matrix
    col_units[10] = "None";    // interp_matrix


    int nrows=1;
    dbg<<"Before Create table"<<std::endl;
    CCfits::Table* table;
    table = fits.addTable("fitpsf",nrows,col_names,col_fmts,col_units);
    dbg<<"After Create table"<<std::endl;


    // Header keywords

    // dimensions information
    dbg<<"Adding dimensional info"<<std::endl;
    std::stringstream tdim10, tdim11;
    tdim10<<"("<<_npca<<","<<nshapelet_coeff<<")";
    tdim11<<"("<<nfit_coeff<<","<<_npca<<")";

    table->addKey("TDIM10",tdim10.str(),"dimensions of rot_matrix");
    table->addKey("TDIM11",tdim11.str(),"dimensions of interp_matrix");


#ifdef USE_TMV
    std::string tmvVersion = tmv::TMV_Version();
#else
    std::string tmvVersion = "Eigen";
#endif
    std::string wlVersion = GetWlVersion();

    table->addKey("tmvvers", tmvVersion, "version of TMV code");
    table->addKey("wlvers", wlVersion, "version of weak lensing code");

    std::string str;
    double dbl;
    int intgr;

    dbg<<"Before Write Par Keys"<<std::endl;

    // if wlserun= is sent we'll put it in the header.  This allows us to 
    // associate some more, possibly complicated, metadata with this file
    if ( _params.keyExists("wlserun") ) {
        WriteParamToTable(_params, table, "wlserun", str);
    }

    WriteParamToTable(_params, table, "noise_method", str);
    WriteParamToTable(_params, table, "dist_method", str);

    WriteParamToTable(_params, table, "fitpsf_order", intgr);
    WriteParamToTable(_params, table, "fitpsf_pca_thresh", dbl);
    dbg<<"After Write Par Keys"<<std::endl;


    // write the data
    int start_row=1;

    // CCfits write() doesn't allow constant scalar arguments, so we have copy
    // out.  Might as well be vectors to simplify the signature

    std::vector<int> psforder_vec(1,_psforder);
    std::vector<double> sigma_vec(1,_sigma);
    std::vector<int> fitorder_vec(1,_fitorder);
    std::vector<int> npca_vec(1,_npca); 
    std::vector<double> xmin(1,_bounds.getXMin());
    std::vector<double> xmax(1,_bounds.getXMax());
    std::vector<double> ymin(1,_bounds.getYMin());
    std::vector<double> ymax(1,_bounds.getYMax());

    table->column(col_names[0]).write(psforder_vec, start_row);
    table->column(col_names[1]).write(sigma_vec, start_row);
    table->column(col_names[2]).write(fitorder_vec, start_row);
    table->column(col_names[3]).write(npca_vec, start_row);
    table->column(col_names[4]).write(xmin, start_row);
    table->column(col_names[5]).write(xmax, start_row);
    table->column(col_names[6]).write(ymin, start_row);
    table->column(col_names[7]).write(ymax, start_row);

    double* ptr;

    ptr = TMV_ptr(*_avepsf);
    table->column(col_names[8]).write(ptr, nshapelet_coeff, nrows, start_row);
#ifdef USE_TMV
    ptr = _mV->ptr();
#else
    ptr = TMV_ptr(*_mV_transpose);
#endif
    table->column(col_names[9]).write(ptr, nrot_matrix, nrows, start_row);
    ptr = TMV_ptr(*_f);
    table->column(col_names[10]).write(ptr, ninterp_matrix, nrows, start_row);
}

void FittedPsf::readFits(std::string file)
{
    // must do this way because of the const thing
    int hdu = GetHdu(_params,"fitpsf",file,2);

    dbg<<"Opening FittedPsf file "<<file<<" at hdu "<<hdu<<std::endl;
    CCfits::FITS fits(file, CCfits::Read);
    if (hdu > 1) fits.read(hdu-1);

    CCfits::ExtHDU& table=fits.extension(hdu-1);

    std::string psforder_col=_params.get("fitpsf_psf_order_col");
    std::string sigma_col=_params.get("fitpsf_sigma_col");
    std::string fitorder_col=_params.get("fitpsf_fit_order_col");
    std::string npca_col=_params.get("fitpsf_npca_col");

    std::string xmin_col=_params.get("fitpsf_xmin_col");
    std::string xmax_col=_params.get("fitpsf_xmax_col");
    std::string ymin_col=_params.get("fitpsf_ymin_col");
    std::string ymax_col=_params.get("fitpsf_ymax_col");

    std::string avepsf_col=_params.get("fitpsf_ave_psf_col");
    std::string rot_matrix_col=_params.get("fitpsf_rot_matrix_col");
    std::string interp_matrix_col=_params.get("fitpsf_interp_matrix_col");

    long start=1;
    long end=1;

    std::vector<int> psforder_vec;
    table.column(psforder_col).read(psforder_vec, start, end);
    _psforder=psforder_vec[0];
    xdbg<<"psforder = "<<_psforder<<std::endl;

    std::vector<double> sigma_vec;
    table.column(sigma_col).read(sigma_vec, start, end);
    _sigma=sigma_vec[0];
    xdbg<<"sigma = "<<_sigma<<std::endl;

    std::vector<int> _fitorder_vec;
    table.column(fitorder_col).read(_fitorder_vec, start, end);
    _fitorder=_fitorder_vec[0];
    xdbg<<"fitorder = "<<_fitorder<<std::endl;

    std::vector<int> npca_vec;
    table.column(npca_col).read(npca_vec, start, end);
    _npca=npca_vec[0];
    xdbg<<"npca = "<<_npca<<std::endl;

    std::vector<double> xmin;
    std::vector<double> xmax;
    std::vector<double> ymin;
    std::vector<double> ymax;

    table.column(xmin_col).read(xmin, start, end);
    xdbg<<"xmin = "<<xmin[0]<<std::endl;
    table.column(xmax_col).read(xmax, start, end);
    xdbg<<"xmax = "<<xmax[0]<<std::endl;
    table.column(ymin_col).read(ymin, start, end);
    xdbg<<"ymin = "<<ymin[0]<<std::endl;
    table.column(ymax_col).read(ymax, start, end);
    xdbg<<"ymax = "<<ymax[0]<<std::endl;

    _bounds.setXMin(xmin[0]);
    _bounds.setXMax(xmax[0]);
    _bounds.setYMin(ymin[0]);
    _bounds.setYMax(ymax[0]);
    xdbg<<"bounds = "<<_bounds<<std::endl;



    // vector columns
    double* dptr=NULL;
    std::valarray<double> dVec;

    const int psfsize = (_psforder+1)*(_psforder+2)/2;
    _avepsf.reset(new DVector(psfsize));

    dVec.resize(0);
    table.column(avepsf_col).read(dVec, 1);
    dptr=TMV_ptr(*_avepsf);
    int dSize = dVec.size();
    Assert(dSize == psfsize);
    for (int j=0; j<dSize; ++j) {
        dptr[j] = dVec[j];
    }

    dVec.resize(0);
    table.column(rot_matrix_col).read(dVec, 1);
#ifdef USE_TMV
    _mV.reset(new tmv::Matrix<double,tmv::RowMajor>(_npca,_avepsf->size()));
    Assert(int(_mV->linearView().size()) == _npca*psfsize);
    dptr=TMV_ptr(*_mV);
#else
    _mV_transpose.reset(new DMatrix(_avepsf->size(),_npca));
    dptr=TMV_ptr(*_mV_transpose);
#endif
    dSize = dVec.size();
    Assert(dSize == _npca*psfsize);
    for (int j=0; j<dSize; ++j) {
        dptr[j] = dVec[j];
    }

    _fitsize = (_fitorder+1)*(_fitorder+2)/2;
    xdbg<<"fitsize = "<<_fitsize<<std::endl;
    _f.reset(new DMatrix(_fitsize,_npca));
#ifdef USE_TMV
    Assert(int(_f->linearView().size()) == _npca*_fitsize);
#endif

    dVec.resize(0);
    table.column(interp_matrix_col).read(dVec, 1);
    dptr=TMV_ptr(*_f);
    dSize = dVec.size();
    Assert(dSize == _npca*_fitsize);
    for (int j=0; j<dSize; ++j) {
        dptr[j] = dVec[j];
    }
}

double FittedPsf::interpolateSingleElement(Position pos, int i) const
{
    DVector P(_fitsize);
#ifdef USE_TMV
    setPRow(_fitorder,pos,_bounds,P.view());
    DVector b1 = P * (*_f);
    double bi = b1 * _mV->col(i,0,_npca);
#else
    setPRow(_fitorder,pos,_bounds,P);
    DVector b1 = _f->transpose() * P;
    double bi = EIGEN_ToScalar(
        EIGEN_Transpose(b1) * _mV_transpose->TMV_rowpart(i,0,_npca));
#endif
    bi += (*_avepsf)(i);
    return bi;
}

void FittedPsf::interpolateVector(Position pos, DVectorView b) const
{
    DVector P(_fitsize);
#ifdef USE_TMV
    setPRow(_fitorder,pos,_bounds,P.view());
    DVector b1 = P * (*_f);
    b = b1 * _mV->rowRange(0,_npca);
#else
    setPRow(_fitorder,pos,_bounds,P);
    DVector b1 = _f->transpose() * P;
    b = TMV_colRange(*_mV_transpose,0,_npca) * b1;
#endif
    b += *_avepsf;
}

BVec FittedPsf::getMean() const
{ return BVec(_psforder,_sigma,*_avepsf); }

