
#include <valarray>
#include <fstream>
#include "TMV.h"
#include <CCfits/CCfits>

#include "FittedPsf.h"
#include "dbg.h"
#include "Function2D.h"
#include "Legendre2D.h"
#include "Name.h"
#include "WlVersion.h"
#include "TMV.h"
#include "WriteParam.h"

static tmv::Vector<double> definePXY(
    int order, double x, double xMin, double xMax)
{
    tmv::Vector<double> temp(order+1);
    double newX = (2.*x-xMin-xMax)/(xMax-xMin);
    temp[0] = 1.;
    if(order>0) temp[1] = newX;
    for(int i=2;i<=order;++i) {
        temp[i] = ((2.*i-1.)*newX*temp[i-1] - (i-1.)*temp[i-2])/i;
    }
    return temp;
}

static void setPRow(
    int fitOrder, Position pos, const Bounds& bounds, 
    const tmv::VectorView<double>& pRow)
{
    Assert(int(pRow.size()) == (fitOrder+1)*(fitOrder+2)/2);
    tmv::Vector<double> px = 
        definePXY(fitOrder,pos.getX(),bounds.getXMin(),bounds.getXMax());
    tmv::Vector<double> py = 
        definePXY(fitOrder,pos.getY(),bounds.getYMin(),bounds.getYMax());
    int pq = 0;
    for(int n=0;n<=fitOrder;++n) {
        for(int p=n,q=n-p;q<=n;--p,++q) {
            Assert(pq < int(pRow.size()));
            pRow(pq) = px[p]*py[q];
            ++pq;
        }
    }
    Assert(pq == int(pRow.size()));
}

FittedPsf::FittedPsf(PsfCatalog& psfCat, const ConfigFile& params) :
    _params(params), _psfOrder(_params.read<int>("psf_order")),
    _fitOrder(_params.read<int>("fitpsf_order")),
    _fitSize((_fitOrder+1)*(_fitOrder+2)/2)
{
    xdbg<<"FittedPSF constructor\n";
    // Do a polynomial fit of the psf shapelet vectors

    Assert(psfCat.size() > 0);

    // This used to be set to the sigma of psfCat.getPsf(0).
    // This doesn't work if the first object failed. 
    // So now I set it to -1, and then update it when we get to the first
    // object without an error flag (usually n = 0).
    _sigma = -1.;
    const int nStars = psfCat.size();
    for(int n=0;n<nStars;++n) if (!psfCat.getFlags(n)) {
        _sigma = psfCat.getPsf(n).getSigma();
        break;
    }

    const int psfSize = (_psfOrder+1)*(_psfOrder+2)/2;
    const double nSigmaClip = _params["fitpsf_nsigma_outlier"];

    // This is an empirical fit to the chisq level that corresponds to 
    // a 3-sigma outlier for more than 1 dimension.
    // I calculated the values up to n=100 and for n>30, they form a pretty
    // good approximation to a straight line.
    // This is almost certainly wrong for nSigmaClip != 3, so if we start 
    // choosing other values for nSigmaClip, it might be worth doing this 
    // right.
    // That means calculating the 1-d critical value for the given nSigma.
    // e.g. nSigma = 3 -> alpha = P(chisq > 9) = 0.0027.
    // Then calculate the critical value of chisq for that alpha with 
    // the full degrees of freedom = psfSize
    const double chisqLevel = 0.14*psfSize + 2.13;

    const double outlierThresh = nSigmaClip * nSigmaClip * chisqLevel;
    dbg<<"outlierThresh = "<<outlierThresh<<std::endl;

    int nOutliers;
    do {

        // Calculate the average psf vector
        _avePsf.reset(new DVector(psfSize));
        Assert(psfSize == int(_avePsf->size()));
        _avePsf->Zero();
        int nGoodPsf = 0;
        for(int n=0;n<nStars;++n) if (!psfCat.getFlags(n)) {
            Assert(psfCat.getPsf(n).getSigma() == _sigma);
            *_avePsf += psfCat.getPsf(n).vec();
            ++nGoodPsf;
        }
        if (nGoodPsf == 0) {
            dbg<<"ngoodpsf = 0 in FittedPsf constructor\n";
            throw ProcessingException("No good stars found for interpolation.");
        }
        *_avePsf /= double(nGoodPsf);

        // Rotate the vectors into their eigen directions.
        // The matrix V is stored to let us get back to the original basis.
        tmv::Matrix<double> mM(nGoodPsf,psfSize);
        tmv::DiagMatrix<double> inverseSigma(nGoodPsf);
        int i=0;
        for(int n=0;n<nStars;++n) if (!psfCat.getFlags(n)) {
            Assert(int(psfCat.getPsf(n).size()) == psfSize);
            Assert(i < nGoodPsf);
            mM.row(i) = psfCat.getPsf(n).vec() - *_avePsf;
            inverseSigma(i) = psfCat.getNu(n);
            _bounds += psfCat.getPos(n);
            ++i;
        }
        Assert(i == nGoodPsf);
        xdbg<<"bounds = "<<_bounds<<std::endl;
        mM = inverseSigma * mM;

        int nPcaTot = std::min(nGoodPsf,psfSize);
        tmv::DiagMatrix<double> mS(nPcaTot);
        tmv::MatrixView<double> mU = mM.Cols(0,nPcaTot);
        _mV.reset(new tmv::Matrix<double,tmv::RowMajor>(nPcaTot,psfSize));
        if (nGoodPsf > psfSize) {
            SV_Decompose(mU.View(),mS.View(),_mV->View(),true);
        } else {
            *_mV = mM;
            SV_Decompose(_mV->Transpose(),mS.View(),mU.Transpose());
        }
        xdbg<<"In FittedPSF: SVD S = "<<mS.diag()<<std::endl;
        if (_params.keyExists("fitpsf_npca")) {
            _nPca = _params["fitpsf_npca"];
            xdbg<<"npca = "<<_nPca<<" from parameter file\n";
        } else {
            double thresh = mS(0);
            if (_params.keyExists("fitpsf_pca_thresh")) 
                thresh *= double(_params["fitpsf_pca_thresh"]);
            else thresh *= tmv::Epsilon<double>();
            xdbg<<"thresh = "<<thresh<<std::endl;
            for(_nPca=1;_nPca<int(mM.rowsize());++_nPca) {
                if (mS(_nPca) < thresh) break;
            }
            xdbg<<"npca = "<<_nPca<<std::endl;
        }
        mU.Cols(0,_nPca) *= mS.SubDiagMatrix(0,_nPca);
        // U S = M(orig) * Vt

        while (nGoodPsf <= _fitSize && _fitSize > 1) {
            --_fitOrder;
            _fitSize = (_fitOrder+1)*(_fitOrder+2)/2;
            dbg<<"Too few good stars... reducing order of fit to "<<
                _fitOrder<<std::endl;
        }
        tmv::Matrix<double> mP(nGoodPsf,_fitSize,0.);
        i=0;
        for(int n=0;n<nStars;++n) if (!psfCat.getFlags(n)) {
            setPRow(_fitOrder,psfCat.getPos(n),_bounds,mP.row(i));
            ++i;
        }
        Assert(i == nGoodPsf);
        mP = inverseSigma * mP;

        _f.reset(new tmv::Matrix<double>(mU.Cols(0,_nPca)/mP));
        xdbg<<"Done making FittedPSF\n";

        //
        // Remove outliers from the fit using the empirical covariance matrix
        // of the data with respect to the fitted values.
        //
        xdbg<<"Checking for outliers:\n";

        // Calculate the covariance matrix
        tmv::Matrix<double> cov(psfSize,psfSize,0.);
        for(int n=0;n<nStars;++n) if (!psfCat.getFlags(n)) {
            const BVec& data = psfCat.getPsf(n);
            DVector fit(psfSize);
            interpolateVector(psfCat.getPos(n),fit.View());
            tmv::Vector<double> diff = data.vec() - fit;
            cov += diff ^ diff;
        }
        if (nGoodPsf > _fitSize) { 
            // only <= if both == 1, so basically always true
            cov /= double(nGoodPsf-_fitSize);
        }
        cov.DivideUsing(tmv::SV);
        cov.SaveDiv();
        cov.SetDiv();
        xdbg<<"cov S = "<<cov.SVD().GetS().diag()<<std::endl;

        // Clip out 3 sigma outliers:
        nOutliers = 0;
        for(int n=0;n<nStars;++n) if (!psfCat.getFlags(n)) {
            const BVec& data = psfCat.getPsf(n);
            DVector fit(psfSize);
            interpolateVector(psfCat.getPos(n),fit.View());
            tmv::Vector<double> diff = data.vec() - fit;
            double dev = diff * cov.Inverse() * diff;
            if (dev > outlierThresh) {
                xdbg<<"n = "<<n<<" is an outlier.\n";
                xdbg<<"data = "<<data.vec()<<std::endl;
                xdbg<<"fit = "<<fit<<std::endl;
                xdbg<<"diff = "<<diff<<std::endl;
                xdbg<<"diff/cov = "<<diff/cov<<std::endl;
                xdbg<<"dev = "<<dev<<std::endl;
                ++nOutliers;
                psfCat.setFlag(n, PSF_INTERP_OUTLIER);
            }
        }
        xdbg<<"ngoodpsf = "<<nGoodPsf<<std::endl;
        xdbg<<"nOutliers = "<<nOutliers<<std::endl;

    } while (nOutliers > 0);

}

FittedPsf::FittedPsf(const ConfigFile& params) : _params(params)
{
    read();
}

// With this one we don't have to use the whole root= thing
FittedPsf::FittedPsf(const ConfigFile& params, std::string file) : 
    _params(params)
{
    read(file);
}



void FittedPsf::writeAscii(std::string file) const
{
    std::ofstream fout(file.c_str());
    fout << _psfOrder <<"  "<< _sigma <<"  ";
    fout << _fitOrder <<"  "<< _nPca << std::endl;
    fout << _bounds << std::endl;
    fout << *_avePsf << std::endl;
    fout << _mV->Rows(0,_nPca) <<std::endl;
    fout << *_f << std::endl;
}

void FittedPsf::readAscii(std::string file)
{
    std::ifstream fin(file.c_str());
    if (!fin) {
        throw ReadException("Error opening fitpsf file"+file);
    }

    xdbg<<"Reading FittedPSF:\n";
    fin >> _psfOrder >> _sigma >> _fitOrder >> _nPca >> _bounds;
    xdbg<<"psforder = "<<_psfOrder<<", sigma_psf = "<<_sigma<<std::endl;
    xdbg<<"fitorder = "<<_fitOrder<<", npca = "<<_nPca<<std::endl;
    xdbg<<"bounds = "<<_bounds<<std::endl;
    _fitSize = (_fitOrder+1)*(_fitOrder+2)/2;
    const int psfSize = (_psfOrder+1)*(_psfOrder+2)/2;
    _avePsf.reset(new DVector(psfSize));
    fin >> *_avePsf;
    _mV.reset(new tmv::Matrix<double,tmv::RowMajor>(_nPca,_avePsf->size()));
    fin >> *_mV;
    _f.reset(new tmv::Matrix<double>(_fitSize,_nPca));
    fin >> *_f;
    if (!fin) {
        throw ReadException("Error reading fitpsf file"+file);
    }
}

void FittedPsf::write() const
{
    std::vector<std::string> fileList = makeMultiName(_params, "fitpsf");  

    const int nFiles = fileList.size();
    for(int i=0; i<nFiles; ++i) {
        const std::string& file = fileList[i];
        dbg<<"Writing fitted psf to file: "<<file<<std::endl;

        bool isFitsIo = false;
        if (_params.keyExists("fitpsf_io")) {
            std::vector<std::string> ios = _params["fitpsf_io"];
            Assert(ios.size() == fileList.size());
            isFitsIo = (ios[i] == "FITS");
        } else if (file.find("fits") != std::string::npos) {
            isFitsIo = true;
        }

        try {
            if (isFitsIo) {
                writeFits(file);
            } else {
                writeAscii(file);
            }
        } catch (CCfits::FitsException& e) {
            throw WriteException(
                "Error writing to "+file+" -- caught error\n" + e.message());
        } catch (std::exception& e) {
            throw WriteException(
                "Error writing to "+file+" -- caught error\n" + e.what());
        } catch (...) {
            throw WriteException(
                "Error writing to "+file+" -- caught unknown error");
        }
    }
    dbg<<"Done Write FittedPSF\n";
}

void FittedPsf::read()
{
    std::string file = makeName(_params,"fitpsf",false,true);
    read(file);
}

void FittedPsf::read(std::string file)
{
    // false,true = input_prefix=false, mustexist=true.
    // It is an input here, but it is in the output_prefix directory.
    dbg<< "Reading FittedPSF from file: " << file << std::endl;

    bool isFitsIo = false;
    if (_params.keyExists("fitpsf_io")) 
        isFitsIo = (_params["fitpsf_io"] == "FITS");
    else if (file.find("fits") != std::string::npos) 
        isFitsIo = true;

    if (!doesFileExist(file)) {
        throw FileNotFoundException(file);
    }
    try {
        if (isFitsIo) {
            readFits(file);
        } else {
            readAscii(file);
        }
    } catch (CCfits::FitsException& e) {
        throw ReadException(
            "Error reading from "+file+" -- caught error\n" + e.message());
    } catch (std::exception& e) {
        throw ReadException(
            "Error reading from "+file+" -- caught error\n" + e.what());
    } catch (...) {
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
    int nShapeletCoeff = (_psfOrder+1)*(_psfOrder+2)/2;
    int nRotMatrix = _nPca*nShapeletCoeff;

    int nFitCoeff = (_fitOrder+1)*(_fitOrder+2)/2;
    int nInterpMatrix = _nPca*nFitCoeff;

    const int nfields=11;

    std::vector<string> colNames(nfields);
    std::vector<string> colFmts(nfields);
    std::vector<string> colUnits(nfields);

    colNames[0] = _params["fitpsf_psf_order_col"];
    colNames[1] = _params["fitpsf_sigma_col"];
    colNames[2] = _params["fitpsf_fit_order_col"];
    colNames[3] = _params["fitpsf_npca_col"];

    colNames[4] = _params["fitpsf_xmin_col"];
    colNames[5] = _params["fitpsf_xmax_col"];
    colNames[6] = _params["fitpsf_ymin_col"];
    colNames[7] = _params["fitpsf_ymax_col"];

    colNames[8] = _params["fitpsf_ave_psf_col"];
    colNames[9] = _params["fitpsf_rot_matrix_col"];
    colNames[10] = _params["fitpsf_interp_matrix_col"];

    colFmts[0] = "1J"; // psf order
    colFmts[1] = "1D"; // sigma
    colFmts[2] = "1J"; // fit order
    colFmts[3] = "1J"; // npca
    colFmts[4] = "1E"; // xmin
    colFmts[5] = "1E"; // xmax
    colFmts[6] = "1E"; // ymin
    colFmts[7] = "1E"; // ymax

    std::stringstream fmt;
    fmt << nShapeletCoeff << "D";
    colFmts[8] = fmt.str(); // average psf shapelets
    fmt.str("");
    fmt << nRotMatrix << "D";
    colFmts[9] = fmt.str(); // rotation matrix
    fmt.str("");
    fmt << nInterpMatrix << "D";
    colFmts[10] = fmt.str(); // interp matrix

    colUnits[0] = "None";     // psf order
    colUnits[1] = "arcsec";   // sigma
    colUnits[2] = "None";     // fit order
    colUnits[3] = "None";     // npca
    colUnits[4] = "pixels";   // xmin
    colUnits[5] = "pixels";   // xmax
    colUnits[6] = "pixels";   // ymin
    colUnits[7] = "pixels";   // ymax
    colUnits[8] = "None";     // avg psf shapelets
    colUnits[9] = "None";     // rot_matrix
    colUnits[10] = "None";    // interp_matrix


    int nRows=1;
    dbg<<"Before Create table"<<std::endl;
    CCfits::Table* table;
    table = fits.addTable("fitpsf",nRows,colNames,colFmts,colUnits);
    dbg<<"After Create table"<<std::endl;


    // Header keywords

    // dimensions information
    dbg<<"Adding dimensional info"<<std::endl;
    std::stringstream tdim10, tdim11;
    tdim10<<"("<<_nPca<<","<<nShapeletCoeff<<")";
    tdim11<<"("<<nFitCoeff<<","<<_nPca<<")";

    table->addKey("TDIM10",tdim10.str(),"dimensions of rot_matrix");
    table->addKey("TDIM11",tdim11.str(),"dimensions of interp_matrix");


    std::string tmvVersion = tmv::TMV_Version();
    std::string wlVersion = getWlVersion();

    table->addKey("tmvvers", tmvVersion, "version of TMV code");
    table->addKey("wlvers", wlVersion, "version of weak lensing code");

    std::string str;
    double dbl;
    int intgr;

    dbg<<"Before Write Par Keys"<<std::endl;

    // if serun= is sent we'll put it in the header.  This allows us to 
    // associate some more, possibly complicated, metadata with this file
    if ( _params.keyExists("serun") ) {
        writeParamToTable(_params, table, "serun", str);
    }

    writeParamToTable(_params, table, "noise_method", str);
    writeParamToTable(_params, table, "dist_method", str);

    writeParamToTable(_params, table, "fitpsf_order", intgr);
    writeParamToTable(_params, table, "fitpsf_pca_thresh", dbl);
    dbg<<"After Write Par Keys"<<std::endl;


    // write the data
    int startRow=1;

    // CCfits write() doesn't allow constant scalar arguments, so we have copy
    // out.  Might as well be vectors to simplify the signature

    std::vector<int> psfOrderV(1,_psfOrder);
    std::vector<double> sigmaV(1,_sigma);
    std::vector<int> fitOrderV(1,_fitOrder);
    std::vector<int> nPcaV(1,_nPca); 
    std::vector<double> xMin(1,_bounds.getXMin());
    std::vector<double> xMax(1,_bounds.getXMax());
    std::vector<double> yMin(1,_bounds.getYMin());
    std::vector<double> yMax(1,_bounds.getYMax());

    table->column(colNames[0]).write(psfOrderV, startRow);
    table->column(colNames[1]).write(sigmaV, startRow);
    table->column(colNames[2]).write(fitOrderV, startRow);
    table->column(colNames[3]).write(nPcaV, startRow);
    table->column(colNames[4]).write(xMin, startRow);
    table->column(colNames[5]).write(xMax, startRow);
    table->column(colNames[6]).write(yMin, startRow);
    table->column(colNames[7]).write(yMax, startRow);

    double* cptr;

    cptr = (double *) _avePsf->cptr();
    table->column(colNames[8]).write(cptr, nShapeletCoeff, nRows, startRow);
    cptr = (double *) _mV->cptr();
    table->column(colNames[9]).write(cptr, nRotMatrix, nRows, startRow);
    cptr = (double *) _f->cptr();
    table->column(colNames[10]).write(cptr, nInterpMatrix, nRows, startRow);
}

void FittedPsf::readFits(std::string file)
{
    // must do this way because of the const thing
    int hdu = getHdu(_params,"fitpsf",file,2);

    dbg<<"Opening FITS file "<<file<<" at hdu "<<hdu<<std::endl;
    CCfits::FITS fits(file, CCfits::Read);
    if (hdu > 1) fits.read(hdu-1);

    CCfits::ExtHDU& table=fits.extension(hdu-1);

    std::string psfOrderCol=_params.get("fitpsf_psf_order_col");
    std::string sigmaCol=_params.get("fitpsf_sigma_col");
    std::string fitOrderCol=_params.get("fitpsf_fit_order_col");
    std::string nPcaCol=_params.get("fitpsf_npca_col");

    std::string xMinCol=_params.get("fitpsf_xmin_col");
    std::string xMaxCol=_params.get("fitpsf_xmax_col");
    std::string yMinCol=_params.get("fitpsf_ymin_col");
    std::string yMaxCol=_params.get("fitpsf_ymax_col");

    std::string avePsfCol=_params.get("fitpsf_ave_psf_col");
    std::string rotMatrixCol=_params.get("fitpsf_rot_matrix_col");
    std::string interpMatrixCol=_params.get("fitpsf_interp_matrix_col");

    long start=1;
    long end=1;

    std::vector<int> psfOrderV;
    table.column(psfOrderCol).read(psfOrderV, start, end);
    _psfOrder=psfOrderV[0];
    xdbg<<"psforder = "<<_psfOrder<<std::endl;

    std::vector<double> sigmaV;
    table.column(sigmaCol).read(sigmaV, start, end);
    _sigma=sigmaV[0];
    xdbg<<"sigma = "<<_sigma<<std::endl;

    std::vector<int> _fitOrderV;
    table.column(fitOrderCol).read(_fitOrderV, start, end);
    _fitOrder=_fitOrderV[0];
    xdbg<<"fitorder = "<<_fitOrder<<std::endl;

    std::vector<int> nPcaV;
    table.column(nPcaCol).read(nPcaV, start, end);
    _nPca=nPcaV[0];
    xdbg<<"npca = "<<_nPca<<std::endl;

    std::vector<float> xMin;
    std::vector<float> xMax;
    std::vector<float> yMin;
    std::vector<float> yMax;

    table.column(xMinCol).read(xMin, start, end);
    xdbg<<"xmin = "<<xMin[0]<<std::endl;
    table.column(xMaxCol).read(xMax, start, end);
    xdbg<<"xmax = "<<xMax[0]<<std::endl;
    table.column(yMinCol).read(yMin, start, end);
    xdbg<<"ymin = "<<yMin[0]<<std::endl;
    table.column(yMaxCol).read(yMax, start, end);
    xdbg<<"ymax = "<<yMax[0]<<std::endl;

    _bounds.setXMin(xMin[0]);
    _bounds.setXMax(xMax[0]);
    _bounds.setYMin(yMin[0]);
    _bounds.setYMax(yMax[0]);
    xdbg<<"bounds = "<<_bounds<<std::endl;



    // vector columns
    double* dptr=NULL;
    std::valarray<double> dVec;

    const int psfSize = (_psfOrder+1)*(_psfOrder+2)/2;
    _avePsf.reset(new DVector(psfSize));
    int nShapeletCoeff = (_psfOrder+1)*(_psfOrder+2)/2;
    Assert(int(_avePsf->size()) == nShapeletCoeff);

    dVec.resize(0);
    table.column(avePsfCol).read(dVec, 1);
    dptr=(double*) _avePsf->cptr();
    int dSize = dVec.size();
    for (int j=0; j<dSize; ++j) {
        dptr[j] = dVec[j];
    }



    int nRotMatrix = _nPca*nShapeletCoeff;
    _mV.reset(new tmv::Matrix<double,tmv::RowMajor>(_nPca,_avePsf->size()));
    Assert(int(_mV->LinearView().size()) == nRotMatrix);

    dVec.resize(0);
    table.column(rotMatrixCol).read(dVec, 1);
    dptr=(double*) _mV->cptr();
    dSize = dVec.size();
    for (int j=0; j<dSize; ++j) {
        dptr[j] = dVec[j];
    }



    _fitSize = (_fitOrder+1)*(_fitOrder+2)/2;
    xdbg<<"fitsize = "<<_fitSize<<std::endl;
    int nInterpMatrix = _nPca*_fitSize;
    _f.reset(new tmv::Matrix<double>(_fitSize,_nPca));
    Assert(int(_f->LinearView().size()) == nInterpMatrix);

    dVec.resize(0);
    table.column(interpMatrixCol).read(dVec, 1);
    dptr=(double*) _f->cptr();
    dSize = dVec.size();
    for (int j=0; j<dSize; ++j) {
        dptr[j] = dVec[j];
    }
}

void FittedPsf::interpolateVector(
    Position pos, const tmv::VectorView<double>& b) const
{
    tmv::Vector<double> P(_fitSize);
    setPRow(_fitOrder,pos,_bounds,P.View());
    tmv::Vector<double> b1 = P * (*_f);
    b = b1 * _mV->Rows(0,_nPca);
    b += *_avePsf;
}

