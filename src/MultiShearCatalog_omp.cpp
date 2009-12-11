
#include "ConfigFile.h"
#include "dbg.h"
#include "Params.h"
#include "Log.h"
#include "TimeVars.h"
#include "MultiShearCatalog.h"
#include "Ellipse.h"
#include "ShearCatalog.h"
#include "ShearCatalogTree.h"
#include "MeasureShearAlgo.h"

int MultiShearCatalog::measureMultiShears(const Bounds& b, ShearLog& log)
{
    dbg<<"Start MeasureMultiShears for b = "<<b<<std::endl;
    int nGals = _skyPos.size();
    dbg<<"ngals = "<<nGals<<std::endl;

    // Read some needed parameters
    double galAperture = _params.read<double>("shear_aperture");
    double maxAperture = _params.read("shear_max_aperture",0.);
    int galOrder = _params.read<int>("shear_gal_order");
    int galOrder2 = _params.read("shear_gal_order2",galOrder);
    double fPsf = _params.read<double>("shear_f_psf");
    double minGalSize = _params.read<double>("shear_min_gal_size");
    bool shouldOutputDots = _params.read("output_dots",false);
    bool isTiming = _params.read("timing",false);

    OverallFitTimes allTimes;

    int nSuccess = 0;

#ifdef ENDAT
    nGals = ENDAT;
#endif

    // Main loop to measure shears
#ifdef _OPENMP
#pragma omp parallel 
    {
        try {
#endif
            OverallFitTimes times; // just for this thread
            ShearLog log1(_params); // just for this thread
            log1.noWriteLog();
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
            for(int i=0;i<nGals;++i) {
                if (!b.includes(_skyPos[i])) continue;
                if (_flags[i]) continue;
#ifdef STARTAT
                if (i < STARTAT) continue;
#endif
#ifdef SINGLEGAL
                if (i < SINGLEGAL) continue;
                if (i > SINGLEGAL) break;
#endif

                if (shouldOutputDots) {
#ifdef _OPENMP
#pragma omp critical (output)
#endif
                    {
                        std::cerr<<"."; std::cerr.flush(); 
                    }
                }

                if (_pixList[i].size() == 0) {
                    dbg<<"no valid single epoch images.\n";
                    _flags[i] = NO_SINGLE_EPOCH_IMAGES;
                    continue;
                }

                // Start with an error code of unknown failure, in case
                // something happens that I don't detect as an error.
                _flags[i] = UNKNOWN_FAILURE;
                long flag1 = 0;

                measureMultiShear(
                    // Input data:
                    _skyPos[i], _pixList[i], _psfList[i],
                    // Parameters:
                    galAperture, maxAperture, galOrder, galOrder2, 
                    fPsf, minGalSize, 
                    // Time stats if desired:
                    isTiming ? &times : 0, 
                    // Log information
                    log1,
                    // Ouput values:
                    _shear[i], _cov[i], _shape[i], _nu[i], flag1);

                _flags[i] = flag1;

                if (!flag1) {
                    dbg<<"Successful shear measurement: "<<
                        _shear[i]<<std::endl;
                    ++nSuccess;
                } else {
                    dbg<<"Unsuccessful shear measurement\n"; 
                }

                if (isTiming) {
                    dbg<<"So far: ns = "<<times._nsGamma;
                    dbg<<",  nf = "<<times._nfNative;
                    dbg<<", "<<times._nfMu<<", "<<times._nfGamma<<std::endl;
                }

            }
#ifdef _OPENMP
#pragma omp critical (add_log)
#endif
            {
                if (isTiming) allTimes += times;
                log += log1;
            }
#ifdef _OPENMP
        } catch (...) {
            // This isn't supposed to happen.
            std::cerr<<"STATUS5BEG Caught some error in parallel region STATUS5END\n";
            exit(1);
        }
    }
#endif

    dbg<<nSuccess<<" successful shear measurements in this pass.\n";
    dbg<<log._nsGamma<<" successful shear measurements so far.\n";

    if (isTiming) {
        dbg<<"From timing structure:\n";
        dbg<<allTimes._nsGamma<<" successful shear measurements, ";
        dbg<<allTimes._nfNative<<" + "<<allTimes._nfMu;
        dbg<<" + "<<allTimes._nfGamma<<" unsuccessful\n";
        std::cerr<<allTimes<<std::endl;
    }
    xdbg<<log<<std::endl;

    return nSuccess;
}

static void getImagePixList(
    std::vector<PixelList>& pixList,
    std::vector<BVec>& psfList,
    std::vector<std::complex<double> >& seShearList,
    std::vector<double>& seSizeList,
    long& inputFlags, int& nImagesFound, int& nImagesGotPix,
    const Position& skyPos,
    const Image<double>& im,
    const Transformation& trans,
    const Transformation& invTrans,
    BVec& psf, const FittedPsf& fitPsf,
    const ShearCatalog& shearCat,
    const ShearCatalogTree& shearCatTree,
    const Image<double>*const weightIm,
    const double noise, const double gain,
    const std::string& skyMethod, const double meanSky, 
    const Image<float>*const skyMap,
    const double galAperture, const double maxAperture)
{
    Assert(psfList.size() == pixList.size());
    Assert(seShearList.size() == pixList.size());
    Assert(seSizeList.size() == pixList.size());

    // Convert ra/dec to x,y in this image

    // First, figure out a good starting point for the nonlinear solver:
    Position pos;
    dbg<<"skypos = "<<skyPos<<std::endl;
    invTrans.transform(skyPos,pos);
    xdbg<<"invtrans(skypos) = "<<pos<<std::endl;

    // Now do the full non-linear solver, which should be pretty fast
    // given the decent initial guess.
    if (!trans.inverseTransform(skyPos, pos) ) {
        dbg << "InverseTransform failed for position "<<skyPos<<".\n";
        dbg << "Initial guess was "<<pos<<".\n";
        inputFlags |= TRANSFORM_EXCEPTION;
    }
    xdbg<<"after exact InverseTransform: pos -> "<<pos<<std::endl;
    if (!(fitPsf.getBounds().includes(pos))) {
        xdbg<<"Reject pos "<<pos<<" not in fitpsf bounds ";
        xdbg<<fitPsf.getBounds()<<std::endl;
        return;
    }

    ++nImagesFound;

    try {
        psf = fitPsf(pos);
    } catch (RangeException& e) {
        xdbg<<"fittedpsf range error: \n";
        xdbg<<"p = "<<pos<<", b = "<<e.getBounds()<<std::endl;
        inputFlags |= FITTEDPSF_EXCEPTION;
        return;
    }

    // Make sure the use of trans in getPixList won't throw:
    try {
        // We don't need to save skyPos.  We just want to catch the range
        // error here, so we don't need to worry about it for dudx, etc.
        Position skyPos1;
        trans.transform(pos,skyPos1);
    } catch (RangeException& e) {
        dbg<<"distortion range error: \n";
        xdbg<<"p = "<<pos<<", b = "<<e.getBounds()<<std::endl;
        inputFlags |= TRANSFORM_EXCEPTION;
        return;
    }

    // Find the nearest object in the shear catalog:
    int nearest = shearCatTree.findNearestTo(pos);

    Assert(psfList.size() == pixList.size());
    Assert(seShearList.size() == pixList.size());
    Assert(seSizeList.size() == pixList.size());

    std::complex<double> seShear = 0.;
    double seSize = 0.;
    double galAp = maxAperture;
    if (std::abs(shearCat.getPos(nearest) - pos) < 1.) { 
        seShear = shearCat.getShear(nearest);
        seSize = shearCat.getShape(nearest).getSigma();
        // If we have a good measurement from the single_epoch image, use
        // that sigma for the size.
        // But expand it by 30% in case we need it.
        galAp = galAperture * seSize * 1.3;
        if (galAp > maxAperture) galAp = maxAperture;
    }

    // Calculate the local sky value.
    long flag = 0;
    double sky=0;
    if (skyMethod == "MEAN") {
        sky = meanSky;
    } else if (skyMethod == "NEAREST") {
        sky = shearCat.getSky(nearest);
    } else {
        Assert(skyMethod == "MAP");
        sky = getLocalSky(*skyMap,pos,trans,galAp,flag);
        // If no pixels in aperture, then
        // a) something is very wrong, but
        // b) use the NEAREST method instead.
        if (flag & BKG_NOPIX) sky = shearCat.getSky(nearest);
        inputFlags |= flag;
    }

    pixList.push_back(PixelList());
#ifdef PIXELLIST_BLOCK
    pixList.back().usePool();
#endif
    xdbg<<"pixlist.size = "<<pixList.size()<<std::endl;

    flag = 0;
    getPixList(
        im,pixList.back(),pos,
        sky,noise,gain,weightIm,trans,galAp,flag);
    xdbg<<"Got pixellist, flag = "<<flag<<std::endl;

    // Make sure not (edge or < 10 pixels) although edge is already
    // checked above
    if (flag == 0) {
        dbg<<"pixlist.size = "<<pixList.size()<<std::endl;
        psfList.push_back(psf);

        // If object in the ShearCatalog is within 1 arcsec of the position
        // then assume it is the correct object.
        // Otherwise put in 0's to indicate we don't have an initial
        // guess for this object (on this image).
        if (std::abs(shearCat.getPos(nearest) - pos) < 1.) { 
            seShearList.push_back(shearCat.getShear(nearest));
            seSizeList.push_back(shearCat.getShape(nearest).getSigma());
        } else {
            seShearList.push_back(0.);
            seSizeList.push_back(0.);
        }
        ++nImagesGotPix;
    } else {
        inputFlags |= flag;
        pixList.pop_back();
    }
}

// Get pixel lists from the file specified in params
void MultiShearCatalog::getImagePixelLists(
    int seIndex, const Bounds& bounds)
{
    dbg<<"Start GetImagePixelLists: se_index = "<<seIndex<<std::endl;

    // If the skybounds for each shear catalog have been saved, then
    // we might be able to skip the ShearCatalog load.
    if (int(_savedSeSkyBounds.size()) > seIndex) {
        Bounds seSkyBounds = _savedSeSkyBounds[seIndex];
        dbg<<"saved bounds for image "<<seIndex<<" = "<<seSkyBounds;
        if (!seSkyBounds.intersects(bounds)) {
            dbg<<"Skipping index "<<seIndex<<
                " because bounds don't intersect\n";
            return;
        }
    }

    // Read the shear catalog
    ShearCatalog shearCat(_params);
    Bounds seSkyBounds = shearCat.getSkyBounds();
    Bounds seBounds = shearCat.getBounds();
    dbg<<"bounds for image "<<seIndex<<" = "<<seSkyBounds;

    // Skip this file if none of the objects in it are in this section of sky.
    if (int(_savedSeSkyBounds.size()) <= seIndex) {
        // Then save the seSkyBounds
        Assert(int(_savedSeSkyBounds.size()) == seIndex);
        _savedSeSkyBounds.push_back(seSkyBounds);
    }
    if (!seSkyBounds.intersects(bounds)) {
        dbg<<"Skipping index "<<seIndex<<" because bounds don't intersect\n";
        return;
    }

    // Read transformation between ra/dec and x/y
    Transformation trans(_params);

    // Read the psf
    FittedPsf fitPsf(_params);

    // Make a tree of the shear catalog to more easily find the nearest
    // single-epoch object to each coadd detection.
    ShearCatalogTree shearCatTree(shearCat);

    // Figure out which method we are going to use to calculate the 
    // local sky values.
    std::string skyMethod = _params.get("multishear_sky_method");
    Assert(skyMethod=="MEAN" || skyMethod=="NEAREST" || skyMethod=="MAP");
    double meanSky=0.;
    std::auto_ptr<Image<float> > skyMap(0);
    if (skyMethod == "MEAN") {
        const int nGals = shearCat.size();
        for(int i=0;i<nGals;++i) meanSky += shearCat.getSky(i);
        meanSky /= shearCat.size();
    }
    if (skyMethod == "MAP") {
        int skyMapHdu = _params.read("skymap_hdu",1);
        skyMap.reset(
            new Image<float>(makeName(_params,"skymap",true,false),skyMapHdu));
    }

    BVec psf(fitPsf.getPsfOrder(), fitPsf.getSigma());

    // Make an inverse transformation that we will use as a starting 
    // point for the more accurate InverseTransform function.
    Transformation invTrans;
    Bounds invBounds = invTrans.makeInverseOf(trans,seBounds,4);
    dbg<<"skybounds = "<<_skyBounds<<std::endl;
    dbg<<"se_skybounds = "<<seSkyBounds<<std::endl;
    dbg<<"se_bounds = "<<seBounds<<std::endl;
    dbg<<"invb = "<<invBounds<<std::endl;

    // We always use the maximum aperture size here, since we don't know
    // how big the galaxy is yet, so we don't know what galAp will be.
    double galAperture = _params.get("shear_aperture");
    double maxAperture = _params.get("shear_max_aperture");

    // Load the image
    // The bounds needed are 
    std::auto_ptr<Image<double> > image;
    std::auto_ptr<Image<double> > weightIm;
    if (_skyBounds.includes(seSkyBounds)) {
        image.reset(new Image<double>(_params,weightIm));
    } else {
        Bounds intersect = invBounds & _skyBounds;
        dbg<<"intersect = invb & skybounds = "<<intersect<<std::endl;

        Bounds subBounds;
        subBounds += invTrans(intersect.get00());
        subBounds += invTrans(intersect.get01());
        subBounds += invTrans(intersect.get10());
        subBounds += invTrans(intersect.get11());

        // Grow bounds by maxAperture
        tmv::SmallMatrix<double,2,2> D;
        trans.getDistortion(subBounds.getCenter(),D);
        double det = std::abs(D.Det());
        double pixScale = sqrt(det); // arcsec/pixel
        subBounds.addBorder(maxAperture / pixScale);

        dbg<<"subb = "<<subBounds<<std::endl;
        image.reset(new Image<double>(_params,weightIm,subBounds));
    }

    // We are using the weight image so the noise and gain are dummy variables
    Assert(weightIm.get());
    double noise = 0.0;
    double gain = 0.0;

    dbg<<"Extracting pixel lists\n";
    // loop over the the objects, if the object falls on the image get
    // the pixel list
    Assert(_skyPos.size() == size());
    Assert(_pixList.size() == size());
    Assert(_psfList.size() == size());
    Assert(_seShearList.size() == size());
    Assert(_seSizeList.size() == size());
    Assert(_flags.size() == size());
    Assert(_inputFlags.size() == size());
    Assert(_nImagesFound.size() == size());
    Assert(_nImagesGotPix.size() == size());

    const int nGals = size();
#ifdef _OPENMP
    //#pragma omp parallel for
#endif
    for (int i=0; i<nGals; ++i) {
        if (_flags[i]) continue;
        if (!bounds.includes(_skyPos[i])) continue;
        if (!invBounds.includes(_skyPos[i])) continue;
        getImagePixList(
            _pixList[i], _psfList[i], _seShearList[i], _seSizeList[i],
            _inputFlags[i], _nImagesFound[i], _nImagesGotPix[i], _skyPos[i], 
            *image, trans, invTrans, psf, fitPsf, shearCat, shearCatTree,
            weightIm.get(), noise, gain,
            skyMethod, meanSky, skyMap.get(),
            galAperture, maxAperture);
    } 
    // Keep track of how much memory we are using.
    // TODO: introduce a parameter max_memory and check to make sure
    // we stay within the allowed memory usage.
    bool shouldOutputDots = _params.read("output_dots",false);
    if (shouldOutputDots) {
        std::cerr<<"Using image# "<<seIndex;
        std::cerr<<"... Memory Usage in MultiShearCatalog = ";
        std::cerr<<calculateMemoryFootprint()<<" MB";
        std::cerr<<"\n";
    }

    dbg<<"Done extracting pixel lists\n";
}

