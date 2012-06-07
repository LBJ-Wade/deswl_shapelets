
#include "dbg.h"
#include "ConfigFile.h"
#include "Params.h"
#include "Log.h"
#include "MultiShearCatalog.h"
#include "Ellipse.h"
#include "ShearCatalog.h"
#include "ShearCatalogTree.h"
#include "MeasureShearAlgo.h"

int MultiShearCatalog::measureMultiShears(const Bounds& b, ShearLog& log)
{
    dbg<<"Start MeasureMultiShears for b = "<<b<<std::endl;
    int ngals = _skypos.size();
    dbg<<"ngals = "<<ngals<<std::endl;

    // Read some needed parameters
    bool output_dots = _params.read("output_dots",false);
#ifdef _OPENMP
    bool des_qa = _params.read("des_qa",false); 
#endif
    bool nativeOnly = _params.read("shear_native_only",false);

    int nSuccess = 0;

#ifdef ENDAT
    ngals = ENDAT;
#endif

    // Main loop to measure shears
#ifdef _OPENMP
#pragma omp parallel 
    {
        try {
#endif
            ShearLog log1(_params); // just for this thread
            log1.noWriteLog();
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
            for(int i=0;i<ngals;++i) {
                if (!b.includes(_skypos[i])) continue;
                if (_flags[i]) continue;
#ifdef STARTAT
                if (i < STARTAT) continue;
#endif
#ifdef SINGLEGAL
                if (i < SINGLEGAL) continue;
                if (i > SINGLEGAL) break;
#endif

                if (output_dots) {
#ifdef _OPENMP
#pragma omp critical (output)
#endif
                    {
                        std::cerr<<"."; std::cerr.flush(); 
                    }
                }
                dbg<<"galaxy "<<i<<":\n";
                dbg<<"id = "<<_id[i]<<std::endl;
                dbg<<"chippos = "<<_chippos[i]<<std::endl;
                dbg<<"skypos = "<<_skypos[i]<<std::endl;

                int nEpoch = _pix_list[i].size();
                if (nEpoch == 0) {
                    dbg<<"no valid single epoch images.\n";
                    dbg<<"FLAG NO_SINGLE_EPOCH_IMAGES\n";
                    _flags[i] = NO_SINGLE_EPOCH_IMAGES;
                    continue;
                }

                dbg<<"Using "<<nEpoch<<" epochs\n";
                Assert(nEpoch == int(_se_num[i].size()));
                Assert(nEpoch == int(_se_pos[i].size()));
                dbg<<"se_image_num   se_image_name   se_pos\n";
                for (int k=0;k<nEpoch;++k) {
                    dbg<<_se_num[i][k]<<"  "<<
                        _image_file_list[_se_num[i][k]]<<"  "<<
                        _se_pos[i][k]<<std::endl;
                }

#if 1
                MeasureSingleShear(
                    // Input data:
                    _pix_list[i], _psf_list[i],
                    // Parameters:
                    _meas_galorder[i], _params,
                    // Log information
                    log1,
                    // Ouput values:
                    _shape[i], _shear[i], _cov[i], _nu[i], _flags[i]);
#else
                _meas_galorder[i] = 0;
                _shear[i] = std::complex<double>(0.1,0.2);
                _cov[i] << 1., 0., 0., 1.;
                _shape[i].vec().setZero();
                _nu[i] = 10.;
#endif

                if (!_flags[i]) {
                    dbg<<"Successful shear measurements: \n";
                    dbg<<"shape = "<<_shape[i]<<std::endl;
                    dbg<<"shear = "<<_shear[i]<<std::endl;
                    dbg<<"cov = "<<_cov[i]<<std::endl;
                    dbg<<"nu = "<<_nu[i]<<std::endl;
                    ++nSuccess;
                } else {
                    dbg<<"Unsuccessful shear measurement\n"; 
                    dbg<<"flag = "<<_flags[i]<<std::endl;
                }
            }
#ifdef _OPENMP
#pragma omp critical (add_log)
#endif
            {
                log += log1;
            }
#ifdef _OPENMP
        } catch (std::exception& e) {
            // This isn't supposed to happen.
            if (des_qa) {
                std::cerr<<"STATUS5BEG Caught error in parallel region STATUS5END\n";
            } 
            std::cerr<<"Caught "<<e.what()<<std::endl;
            std::cerr<<"Caught error in parallel region.  Aborting.\n";
            exit(1);
        } catch (...) {
            if (des_qa) {
                std::cerr<<"STATUS5BEG Caught error in parallel region STATUS5END\n";
            }
            std::cerr<<"Caught error in parallel region.  Aborting.\n";
            exit(1);
        }
    }
#endif

    dbg<<nSuccess<<" successful shear measurements in this pass.\n";
    dbg<<log._ns_gamma<<" successful shear measurements so far.\n";

    // Get value of VmPeak at this point in program.
    double peak_mem = peak_memory_usage(dbgout);
    dbg<<"Peak memory usage so far = "<<peak_mem<<std::endl;
    _params["peak_mem"] = peak_mem; // We'll write it to output file.

    xdbg<<log<<std::endl;

    return nSuccess;
}

static void getImagePixList(
    std::vector<PixelList>& pix_list,
    std::vector<BVec>& psf_list,
    std::vector<int>& se_num, std::vector<Position>& se_pos, int se_index,
    long& input_flags, int& nimages_found, int& nimages_gotpix,
    const Position& skypos,
    const Image<double>& im,
    const Transformation& trans,
    const Transformation& inv_trans,
    BVec& psf, const FittedPsf& fitpsf,
    const ShearCatalog& shearcat,
    const ShearCatalogTree& shearcat_tree,
    const Image<double>*const weight_image,
    const double noise, const double mean_sky, 
    const Image<double>*const skymap,
    double gal_aperture, double max_aperture, const ConfigFile& params)
{
    std::string sky_method = params.get("multishear_sky_method");
    Assert(sky_method=="MEAN" || sky_method=="NEAREST" || sky_method=="MAP");
    bool require_match = params.read("multishear_require_match",false);

    Assert(psf_list.size() == pix_list.size());
    Assert(se_num.size() == pix_list.size());
    Assert(se_pos.size() == pix_list.size());

    // Convert ra/dec to x,y in this image

    // First, figure out a good starting point for the nonlinear solver:
    Position pos;
    xdbg<<"skypos = "<<skypos<<std::endl;
    inv_trans.transform(skypos,pos);
    xdbg<<"invtrans(skypos) = "<<pos<<std::endl;

    // Now do the full non-linear solver, which should be pretty fast
    // given the decent initial guess.
    if (!trans.inverseTransform(skypos, pos) ) {
        dbg << "InverseTransform failed for position "<<skypos<<".\n";
        dbg << "Initial guess was "<<pos<<".\n";
        input_flags |= TRANSFORM_EXCEPTION;
    }
    xdbg<<"after exact InverseTransform: pos -> "<<pos<<std::endl;
    if (!(fitpsf.getBounds().includes(pos))) {
        xdbg<<"Reject pos "<<pos<<" not in fitpsf bounds ";
        xdbg<<fitpsf.getBounds()<<std::endl;
        return;
    }

    ++nimages_found;

    try {
        psf = fitpsf(pos);
    } catch (RangeException& e) {
        xdbg<<"fittedpsf range error: \n";
        xdbg<<"p = "<<pos<<", b = "<<e.getBounds()<<std::endl;
        input_flags |= FITTEDPSF_EXCEPTION;
        return;
    }

    // Make sure the use of trans in GetPixList won't throw:
    try {
        // We don't need to save skypos.  We just want to catch the range
        // error here, so we don't need to worry about it for dudx, etc.
        Position skypos1;
        trans.transform(pos,skypos1);
    } catch (RangeException& e) {
        dbg<<"distortion range error: \n";
        xdbg<<"p = "<<pos<<", b = "<<e.getBounds()<<std::endl;
        input_flags |= TRANSFORM_EXCEPTION;
        return;
    }

    xdbg<<"Start getImagePixList: mem = "<<memory_usage()<<std::endl;

    // Find the nearest object in the shear catalog:
    int nearest = shearcat_tree.findNearestTo(pos);

    double galap = max_aperture;
    if (std::abs(shearcat.getPos(nearest) - pos) < 1.) { 
        double se_size = shearcat.getShape(nearest).getSigma();
        dbg<<"Single epoch id, pos, shear, size = "<<
            shearcat.getId(nearest)<<"  "<<
            shearcat.getPos(nearest)<<"  "<<
            shearcat.getShear(nearest)<<"  "<<
            se_size<<std::endl;
        // If we have a good measurement from the single_epoch image, use
        // that sigma for the size.
        // But expand it by 30% in case we need it.
        galap = gal_aperture * se_size * 1.3;
        if (galap > max_aperture) galap = max_aperture;
    } else if (require_match) {
        xdbg<<"No match found in single epoch catalog\n";
        // This isn't really the same meaning of this flag in the
        // normal flags value, but it doesn't seem worth making a 
        // new name.  Here it means that an image was rejected because
        // there was no valid measurement in the single-epoch run.
        input_flags |= NO_SINGLE_EPOCH_IMAGES;
        return;
    }

    // Calculate the local sky value.
    long flag = 0;
    double sky=0;
    if (sky_method == "MEAN") {
        sky = mean_sky;
    } else if (sky_method == "NEAREST") {
        sky = shearcat.getSky(nearest);
    } else {
        Assert(sky_method == "MAP");
        sky = GetLocalSky(*skymap,pos,trans,galap,params,flag);
        // If no pixels in aperture, then
        // a) something is very wrong, but
        // b) use the NEAREST method instead.
        if (flag & BKG_NOPIX) sky = shearcat.getSky(nearest);
        input_flags |= flag;
    }

    xdbg<<"Before push_back(PixelLise) mem = "<<memory_usage()<<std::endl;
    xdbg<<"pixlist.size = "<<pix_list.size()<<std::endl;
    pix_list.push_back(PixelList());
    xdbg<<"pixlist.size => "<<pix_list.size()<<std::endl;
    xdbg<<"After push_back(PixelLise) mem = "<<memory_usage()<<std::endl;
#ifdef PIXELLIST_BLOCK
    pix_list.back().usePool();
#endif

    flag = 0;
    xdbg<<"Before GetPixList mem = "<<memory_usage()<<std::endl;
    GetPixList(
        im,pix_list.back(),pos,
        sky,noise,weight_image,trans,galap,params,flag);
    xdbg<<"Got pixellist, flag = "<<FlagText(flag)<<std::endl;
    xdbg<<"After GetPixList mem = "<<memory_usage()<<std::endl;
    
    // Make sure not (edge or < 10 pixels) although edge is already
    // checked above
    if (flag == 0) {
        xdbg<<"Before psf_list.push_back(psf): mem = "<<memory_usage()<<std::endl;
        xdbg<<"psflist.size = "<<psf_list.size()<<std::endl;
        psf_list.push_back(psf);
        xdbg<<"psflist.size => "<<psf_list.size()<<std::endl;
        xdbg<<"After psf_list.push_back(psf): mem = "<<memory_usage()<<std::endl;
        se_num.push_back(se_index);
        se_pos.push_back(pos);
        ++nimages_gotpix;
    } else {
        input_flags |= flag;
        xdbg<<"Before pix_list.pop_back: mem = "<<memory_usage()<<std::endl;
        xdbg<<"pixlist.size = "<<pix_list.size()<<std::endl;
        pix_list.pop_back();
        xdbg<<"pixlist.size => "<<pix_list.size()<<std::endl;
        xdbg<<"After pix_list.pop_back: mem = "<<memory_usage()<<std::endl;
    }
    xdbg<<"Done getImagePixList: mem = "<<memory_usage()<<std::endl;
}

// Get pixel lists from the file specified in params
bool MultiShearCatalog::getImagePixelLists(
    int se_index, const Bounds& bounds)
{
    dbg<<"Start GetImagePixelLists: se_index = "<<se_index<<std::endl;

    // If the skybounds for each shear catalog have been saved, then
    // we might be able to skip the ShearCatalog load.
    if (int(_saved_se_skybounds.size()) > se_index) {
        Bounds se_skybounds = _saved_se_skybounds[se_index];
        dbg<<"saved bounds for image "<<se_index<<" = "<<se_skybounds;
        if (!se_skybounds.intersects(bounds)) {
            dbg<<"Skipping index "<<se_index<<
                " because bounds don't intersect\n";
            return true;
        }
    }

    // Read the shear catalog
    ShearCatalog shearcat(_params);
    shearcat.read();
    Bounds se_skybounds = shearcat.getSkyBounds();
    Bounds se_bounds = shearcat.getBounds();
    dbg<<"bounds for image "<<se_index<<" = "<<se_skybounds;

    // Skip this file if none of the objects in it are in this section of sky.
    if (int(_saved_se_skybounds.size()) <= se_index) {
        // Then save the se_skybounds
        Assert(int(_saved_se_skybounds.size()) == se_index);
        _saved_se_skybounds.push_back(se_skybounds);
    }
    if (!se_skybounds.intersects(bounds)) {
        dbg<<"Skipping index "<<se_index<<" because bounds don't intersect\n";
        return true;
    }

    // Read transformation between ra/dec and x/y
    Transformation trans(_params);

    // Read the psf
    FittedPsf fitpsf(_params);
    fitpsf.read();

    // Make a tree of the shear catalog to more easily find the nearest
    // single-epoch object to each coadd detection.
    ShearCatalogTree shearcat_tree(shearcat);

    // Figure out which method we are going to use to calculate the 
    // local sky values.
    std::string sky_method = _params.get("multishear_sky_method");
    Assert(sky_method=="MEAN" || sky_method=="NEAREST" || sky_method=="MAP");
    double mean_sky=0.;
    std::auto_ptr<Image<double> > skymap(0);
    if (sky_method == "MEAN") {
        const int ngals = shearcat.size();
        for(int i=0;i<ngals;++i) mean_sky += shearcat.getSky(i);
        mean_sky /= shearcat.size();
    }
    if (sky_method == "MAP") {
        std::string skymap_name = MakeName(_params,"skymap",true,true);
        int skymap_hdu = GetHdu(_params,"skymap",skymap_name,1);
        skymap.reset(new Image<double>(skymap_name,skymap_hdu));
    }

    BVec psf(fitpsf.getPsfOrder(), fitpsf.getSigma());

    // Make an inverse transformation that we will use as a starting 
    // point for the more accurate InverseTransform function.
    Transformation inv_trans;
    Bounds inv_bounds = inv_trans.makeInverseOf(trans,se_bounds,4);
    dbg<<"skybounds = "<<_skybounds<<std::endl;
    dbg<<"se_skybounds = "<<se_skybounds<<std::endl;
    dbg<<"se_bounds = "<<se_bounds<<std::endl;
    dbg<<"invb = "<<inv_bounds<<std::endl;

    // We always use the maximum aperture size here, since we don't know
    // how big the galaxy is yet, so we don't know what galap will be.
    double gal_aperture = _params.get("shear_aperture");
    double max_aperture = _params.get("shear_max_aperture");

    // Load the image
    // The bounds needed are 
    std::auto_ptr<Image<double> > image;
    std::auto_ptr<Image<double> > weight_image;
    if (_skybounds.includes(se_skybounds)) {
        image.reset(new Image<double>(_params,weight_image));
    } else if (!_skybounds.intersects(inv_bounds)) {
        dbg<<"Skipping index "<<se_index<<" because inv_bounds doesn't intersect\n";
        return true;
    } else {
        Bounds intersect = inv_bounds & _skybounds;
        dbg<<"intersect = invb & skybounds = "<<intersect<<std::endl;

        Bounds sub_bounds;
        sub_bounds += inv_trans(intersect.get00());
        sub_bounds += inv_trans(intersect.get01());
        sub_bounds += inv_trans(intersect.get10());
        sub_bounds += inv_trans(intersect.get11());

        // Grow bounds by max_aperture
        DSmallMatrix22 D;
        trans.getDistortion(sub_bounds.getCenter(),D);
        double det = std::abs(D.TMV_det());
        double pixel_scale = sqrt(det); // arcsec/pixel
        sub_bounds.addBorder(max_aperture / pixel_scale);

        dbg<<"subb = "<<sub_bounds<<std::endl;
        image.reset(new Image<double>(_params,weight_image,sub_bounds));
    }

    // We are using the weight image so the noise and gain are dummy variables
    Assert(weight_image.get());
    double noise = 0.0;
    double gain = 0.0;

    dbg<<"Extracting pixel lists\n";
    // loop over the the objects, if the object falls on the image get
    // the pixel list
    Assert(int(_skypos.size()) == size());
    Assert(int(_pix_list.size()) == size());
    Assert(int(_psf_list.size()) == size());
    Assert(int(_flags.size()) == size());
    Assert(int(_input_flags.size()) == size());
    Assert(int(_nimages_found.size()) == size());
    Assert(int(_nimages_gotpix.size()) == size());

    const int ngals = size();
    double max_mem = _params.read("max_vmem",64)*1024.;
    bool require_match = _params.read("multishear_require_match",false);
    dbg<<"Before getImagePixList loop: memory_usage = "<<memory_usage()<<std::endl;
    for (int i=0; i<ngals; ++i) {
        if (_flags[i]) continue;
        if (!bounds.includes(_skypos[i])) continue;
        if (!inv_bounds.includes(_skypos[i])) continue;
        dbg<<"getImagePixList for galaxy "<<i<<", id = "<<_id[i]<<std::endl;
        getImagePixList(
            _pix_list[i], _psf_list[i],
            _se_num[i], _se_pos[i], se_index,
            _input_flags[i], _nimages_found[i], _nimages_gotpix[i], _skypos[i], 
            *image, trans, inv_trans, psf, fitpsf, shearcat, shearcat_tree,
            weight_image.get(), noise, mean_sky, skymap.get(),
            gal_aperture, max_aperture, _params);
        double mem = memory_usage();
        // memory_usage is in KB
        // max_vmem is in GB
        if (mem > max_mem) {
            dbg<<"VmSize = "<<mem<<" > max_vmem = "<<max_mem<<std::endl;
            return false;
        }
    }
    // Keep track of how much memory we are using.
    dbg<<"Done getImagePixList loop: memory_usage = "<<memory_usage()<<std::endl;
    if (dbgout) PixelList::dumpPool(*dbgout);
    dbg<<"Using image# "<<se_index;
    dbg<<"... Memory Usage in MultiShearCatalog = ";
    dbg<<calculateMemoryFootprint()<<" MB";
    dbg<<"\n";
    dbg<<"VMem = "<<memory_usage()<<std::endl;
    double peak_mem = peak_memory_usage();
    dbg<<"VPeak = "<<peak_mem<<std::endl;
    dbg<<"max_mem allowed = "<<max_mem<<std::endl;
    _params["peak_mem"] = peak_mem;  // Keep track

    dbg<<"Done extracting pixel lists\n";
    return true;
}

template <typename T>
inline long long getMemoryFootprint(const T& x)
{
    return sizeof(x); 
}

template <typename T>
inline long long getMemoryFootprint(T*const x)
{
    long long res=sizeof(x);
    res += getMemoryFootprint(*x);
    return res;
}

template <typename T>
inline long long getMemoryFootprint(const T*const x)
{
    long long res=sizeof(x);
    res += getMemoryFootprint(*x);
    return res;
}

template <typename T>
inline long long getMemoryFootprint(const std::auto_ptr<T>& x)
{
    long long res=sizeof(x);
    res += getMemoryFootprint(*x);
    return res;
}

template <typename T, class Alloc>
inline long long getMemoryFootprint(const std::vector<T,Alloc>& x)
{
    long long res=sizeof(x);
    const int xsize = x.size();
    for(int i=0;i<xsize;++i) res += getMemoryFootprint(x[i]);
    return res;
}

inline long long getMemoryFootprint(const PixelList& x)
{
    long long res=sizeof(x);
    res += x.size() * sizeof(Pixel);
    return res;
}

template <typename T>
inline long long getMemoryFootprint(const TVector(T)& x)
{
    long long res=sizeof(x);
    res += x.size() * sizeof(T);
    return res;
}

inline long long getMemoryFootprint(const BVec& x)
{
    long long res=sizeof(x);
    res += x.size() * sizeof(double);
    return res;
}

template <typename T>
inline long long getMemoryFootprint(const TMatrix(T)& x)
{
    long long res=sizeof(x);
    res += x.colsize() * x.rowsize() * sizeof(T);
    return res;
}

double MultiShearCatalog::calculateMemoryFootprint() const
{
    dbg<<"Memory usage:\n";
    dbg<<"input_flags: "<<
        getMemoryFootprint(_input_flags)/1024./1024.<<" MB\n";
    dbg<<"nimages_found: "<<
        getMemoryFootprint(_nimages_found)/1024./1024.<<" MB\n";
    dbg<<"pixlist: "<<
        getMemoryFootprint(_pix_list)/1024./1024.<<" MB\n";
    dbg<<"nimages_gotpix: "<<
        getMemoryFootprint(_nimages_gotpix)/1024./1024.<<" MB\n";
    dbg<<"id: "<<
        getMemoryFootprint(_id)/1024./1024.<<" MB\n";
    dbg<<"skypos: "<<
        getMemoryFootprint(_skypos)/1024./1024.<<" MB\n";
    dbg<<"flags: "<<
        getMemoryFootprint(_flags)/1024./1024.<<" MB\n";
    dbg<<"shear: "<<
        getMemoryFootprint(_shear)/1024./1024.<<" MB\n";
    dbg<<"nu: "<<
        getMemoryFootprint(_nu)/1024./1024.<<" MB\n";
    dbg<<"cov: "<<
        getMemoryFootprint(_cov)/1024./1024.<<" MB\n";
    dbg<<"shape: "<<
        getMemoryFootprint(_shape)/1024./1024.<<" MB\n";
    dbg<<"psflist: "<<
        getMemoryFootprint(_psf_list)/1024./1024.<<" MB\n";
    dbg<<"image_file_list: "<<
        getMemoryFootprint(_image_file_list)/1024./1024.<<" MB\n";
    dbg<<"shear_file_list: "<<
        getMemoryFootprint(_shear_file_list)/1024./1024.<<" MB\n";
    dbg<<"fitpsf_file_list: "<<
        getMemoryFootprint(_fitpsf_file_list)/1024./1024.<<" MB\n";
    dbg<<"skymap_file_list: "<<
        getMemoryFootprint(_skymap_file_list)/1024./1024.<<" MB\n";
    dbg<<"saved_se_skybounds: "<<
        getMemoryFootprint(_saved_se_skybounds)/1024./1024.<<" MB\n";

    double totmem = 
        getMemoryFootprint(_input_flags) +
        getMemoryFootprint(_nimages_found) +
        getMemoryFootprint(_pix_list) +
        getMemoryFootprint(_nimages_gotpix) +
        getMemoryFootprint(_id) +
        getMemoryFootprint(_skypos) +
        getMemoryFootprint(_flags) +
        getMemoryFootprint(_shear) +
        getMemoryFootprint(_nu) +
        getMemoryFootprint(_cov) +
        getMemoryFootprint(_shape) +
        getMemoryFootprint(_psf_list) +
        getMemoryFootprint(_image_file_list) +
        getMemoryFootprint(_shear_file_list) +
        getMemoryFootprint(_fitpsf_file_list) +
        getMemoryFootprint(_skymap_file_list) +
        getMemoryFootprint(_saved_se_skybounds);
    totmem /= 1024.*1024.;  // B -> MB
    dbg<<"totmem = "<<totmem<<std::endl;

    return totmem;
}

