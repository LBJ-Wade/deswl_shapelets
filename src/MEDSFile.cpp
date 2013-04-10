
#include <valarray>
#include <CCfits/CCfits>

#include "MEDSFile.h"
#include "dbg.h"
#include "Params.h"
#include "Name.h"

MEDSFile::MEDSFile(const ConfigFile& params) : _params(params)
{
    std::string meds_file=_params.get("meds_file");
    _medsPtr = meds_open(meds_file.c_str());

    const struct meds_info_cat* info_cat = _medsPtr->image_info;
    _fitpsf.resize(info_cat->size);
    // Skip i=0, since that's the coadd image.
    for (int i=1;i<info_cat->size;++i) {
        const char* im_name = info_cat->data[i].image_path;
        SetRoot(_params, im_name);
        _fitpsf[i].reset(new FittedPsf(_params));
    }
}

MEDSFile::~MEDSFile() 
{ 
    _medsPtr = meds_free(_medsPtr);
}

std::string MEDSFile::getCoaddCatFile() const
{
    return meds_get_meta(_medsPtr)->coaddcat_file;
}

void MEDSFile::getPixels(std::vector<PixelList>& pix_list, long i, long& flag) const
{
    dbg<<"Getting pixels for i="<<i<<std::endl;

    double noise = 0.;  // dummy variable
    double sky = 0.;  // MEDS images are sky subtracted.
    double aperture = _params.get("shear_max_aperture");

    // meds uses malloc and free, so tell shared_ptr to use free as the deleter:
    long ncutout,nrow,ncol;
    boost::shared_ptr<double> pix(meds_get_mosaicp(_medsPtr,i,&ncutout,&nrow,&ncol),free);
    long npix = nrow * ncol;

    boost::shared_ptr<double> wtpix(meds_get_weight_mosaicp(_medsPtr,i,&ncutout,&nrow,&ncol),free);

    // We'll convert these one at a time to an Image object.
    Image<double> im(ncol,nrow);
    Image<double> wtim(ncol,nrow);

#if 0
    boost::shared_ptr<int> segpix(meds_get_seg_mosaicp(_medsPtr,i,&ncutout,&nrow,&ncol),free);
    Image<int> segim(ncol,nrow);
#endif

    pix_list.resize(ncutout-1);

    double cen_row,cen_col;
    Transformation trans;
    
    for(int j=1;j<ncutout;++j) {
        dbg<<"j = "<<j<<std::endl;
        int id = meds_get_source_file_id(_medsPtr, i, j);
        dbg<<"id = "<<id<<std::endl;
        Assert(id != 0);  // id = 0 should be at j=0.
        if (id == 0) continue; // Skip id = 0, which is the coadd image.

        // Get the images
        im.getM() = tmv::MatrixViewOf(pix.get()+j*npix, ncol, nrow, tmv::ColMajor);
        wtim.getM() = tmv::MatrixViewOf(wtpix.get()+j*npix, ncol, nrow, tmv::ColMajor);
#if 0
        segim.getM() = tmv::MatrixViewOf(segpix.get()+j*npix,ncol,nrow,tmv::ColMajor);
#endif

        // Get the location of the center pixel
        meds_get_cutout_cen(_medsPtr, i, j, &cen_row, &cen_col);
        Position cen(cen_col, cen_row);

        // Convert the distortion into a Transformation
        const struct meds_distort* dist = meds_get_distortion(_medsPtr, i, j);
        trans.setToJacobian(dist->dudcol, dist->dudrow, dist->dvdcol, dist->dvdrow);

        // Finally, build the PixelList.
        GetPixList(im,pix_list[j-1],cen,sky,noise,&wtim,trans,aperture,_params,flag);
    }
}

void MEDSFile::getPSFs(std::vector<BVec>& psf_list, long i, long& flag) const
{
    dbg<<"Getting psfs for i="<<i<<std::endl;
    long ncutout = meds_get_ncutout(_medsPtr, i);
    psf_list.reserve(ncutout-1);

    double cen_row, cen_col;
    for(int j=1;j<ncutout;++j) {
        dbg<<"j = "<<j<<std::endl;
        int id = meds_get_source_file_id(_medsPtr, i, j);
        dbg<<"id = "<<id<<std::endl;
        Assert(id != 0);  // id = 0 should be at j=0.
        if (id == 0) continue; // Skip id = 0, which is the coadd image.

        // Get the location of the center pixel
        meds_get_cutout_cen(_medsPtr, i, j, &cen_row, &cen_col);
        Position cen(cen_col, cen_row);

        // Get the interpolated BVec
        BVec psf = (*_fitpsf[meds_get_source_file_id(_medsPtr, i, j)])(cen);
    }
}

