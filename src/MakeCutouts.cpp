/*
   Author: Erin Sheldon, BNL

   usage
   ---------
   make-cutouts coadd_file= cat_file= coadd_srclist= cutout_file=

   required keywords
   -----------------
   coadd_file
       The coadd image file.  Holds the image and weight image.
   cat_file
       The meds input file.  This is an ascii file with the following
       columns.  Below are the types the columns are read into
           id ra dec row col box_size
           i8 f8 f8  f8  f8  i8
   coadd_srclist
       An ascii file with a list of the coadd input images.  it has the
       following columns.  Below are the types these are read into.

           image_id flags image  background_image segmentation_image zeropoint (wcs_file)
           i8       i8    string string           string             f8        (string)

       The image is the single-epoch "red" image and background is the "bkg"
       file.  The segmentation image holds the sextractor segmentation map.
       The image_id field is arbitrary; for DES it is generally set to the
       image database id.

       Note the wcs_file is an optional separate file to hold the wcs solution.
       Send use_alt_wcs=1 to specify that column will exist, and se_wcs_hdu
       to specify the HDU.

   cutout_file
       The name of the output file.

   More keywords can be given, some of which are used in the code to over-ride
   defaults.  These are the se_hdu, se_wt_hdu, se_badpix_hdu, sky_hdu,
   coadd_hdu, and coadd_wt_hdu.  Any keywords given will be written to the
   metadata table.
   
   algorithm
   ---------

   Make a first pass through the files using only the wcs to determine the
   cutouts.  Then pass through writing the cutouts one at a time into the
   pre-made mosaic image.  This limits the memory usage significantly.
   Memory usage is typically around 100-200Mb.

   The cutouts are sky subtracted.  The box sizes are forced to be even.

   outputs
   -------

   The output file name is determined by the cutout_file= option.

   The file has the following extensions.  They are all named as follows;
   the order in the file may change, so use the names.
       object_data: 
           A binary table with one row for each coadd object.  This table
           aligns exactly row-for-row with the input catalog file.
       image_info:
           A binary table holding paths to the images used as well as the
           background files, as well as the zeropoints and applied scale
           factors.
       metadata: 
           A binary table holding all parameters used in making the files, as
           well as the value of all keywords sent to the program, whether used
           or not.  This is good for adding other information such as code
           versions.  Also stored is the reference magnitude zero point for
           rescaling the SE images.
       image_cutouts:
           An image extension holding the image cutouts.  The images are sky
           subtracted. All SE images are put on the same zeropoint relative to
           the median zeropoint of single-epoch images.  This reference
           zeropoint is also stored in the metadata table. The coadd is not
           rescaled.
       weight_cutouts:
           An image extension holding the weight image cutouts.  Note the
           same scaling for single-epoch images is applied as for the cutouts,
           i.e. the weight is scaled by 1/scale^2
       seg_cutouts:
           An image extension holding the segmentation cutouts.

   The "object_data" extension is a table with an entry for each object in the
   coadd.  This table lines up row-by-row with the input catalog.

     id                 i8       The id column from the input file
     number             i8       SExtractor "NUMBER" column from catalog
     ncutout            i8       number of cutouts for this object
     box_size           i8       box size for each cutout
     file_id            i8[NMAX] zero-offset id into the file names in the 
                                 second extension
     start_row          i8[NMAX] zero-offset, points to start of each cutout.
     orig_row           f8[NMAX] zero-offset position in original image
     orig_col           f8[NMAX] zero-offset position in original image
     orig_start_row     i8[NMAX] zero-offset start corner in original image
     orig_start_col     i8[NMAX] zero-offset start corner in original image
     cutout_row         f8[NMAX] zero-offset position in cutout imag
     cutout_col         f8[NMAX] zero-offset position in cutout image
     dudrow             f8[NMAX] jacobian of transformation 
     dudcol             f8[NMAX] row,col->ra,dec tangent plane (u,v)
     dvdrow             f8[NMAX] Note Source Extractor uses y for row
     dvdcol             f8[NMAX]

   The array fields are constant size NMAX, where NMAX is the max number of
   cutouts of any object in the list.  When a value is not used, it is set to
   -9999 except for box_size, which is always the box_size as input from
   the cat_file, and ncutout which is the number of cutouts and is zero when
   there are none.

   The first cutout entry for each object is always the coadd cutout.
   
   Note in cfitsio you will need to convert start_row to 1-offset.

   The extension "image_info", contains info about each image.  Currently this
   is just the file paths.  The file_id column above point into this structure.

     image_id             i8       an id for the image
     image_flags          i8       flag bitmask for the image
     image_path           SXX      full path to source images
     wcs_path             SXX      full path to wcs file
     sky_path             SYY      full path to sky images
     seg_path             SZZ      full path to segmentation images
     magzp                f4       magnitude zero point
     scale                f4       scale factor used to place on common system

   The extension "image_cutouts", contains a mosaic of all the image cutouts.
   This mosaic is a single-dimensional image rather than two-dimensional in
   order to support variable cutout sizes.  The zero-offset starting row in
   this image is given in the start_row field in the object_data table.  There
   is an entry for each cutout.  The cutouts for a single object are always
   contiguous, and the cutouts are ordered in the file as objects are ordered
   in the table, which makes sequential access more efficient. The total number
   of pixels for all cutouts associated with an object is
   box_size*box_size*ncutout.

   The extension "weight_cutouts", is the same format as the "image_cutouts"
   extension, but contains the cutouts of the weight images.

   The extension "seg_cutouts", is a similar format as the "image_cutouts"
   extension, but contains the cutouts of the seg images.

   TODO
       - incorporate Daniel's feature map.
*/

// CCfits has some unused variables, so use this to suppress warnings that would otherwise 
// show up.

// icpc pretends to be GNUC, since it thinks it's compliant, but it's not.
// It doesn't understand "pragma GCC"
#ifndef __INTEL_COMPILER

// I think this starts being necessary at version 4.5
#if defined(__GNUC__) && __GNUC__ >= 4 && (__GNUC__ >= 5 || __GNUC_MINOR__ >= 5)
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif

#endif

#include <valarray>
#include <sys/time.h>
#include <iostream>
#include <fstream>
#include <CCfits/CCfits>
#include <stdint.h>
#include <algorithm>

#include "Image.h"
#include "Transformation.h"
#include "BasicSetup.h"

using std::auto_ptr;
using std::cerr;
using std::string;
using std::vector;

typedef float cutout_t;
typedef int badpix_t;
typedef int seg_t;

const double ARCSEC_PER_RAD = 206264.806247;
const int DEFVAL=-9999;

// defaults
const int SE_HDU        = 2;
const int SE_BADPIX_HDU = 3;
const int SE_WT_HDU     = 4;

const int SKY_HDU       = 2;
//const int SEG_HDU       = 1;
// now second hdu for fz files
const int SEG_HDU       = 2;

const int COADD_HDU     = 2;
const int COADD_WT_HDU  = 3;
//const int COADD_SEG_HDU = 1;
const int COADD_SEG_HDU = 2;

enum cutout_type {
    CUTOUT_IMAGE,
    CUTOUT_WEIGHT,
    CUTOUT_BADPIX,
    CUTOUT_SKY,
    CUTOUT_SEG
};

class CutoutMaker
{
    public:
        CutoutMaker(const ConfigFile *input_params);

        ~CutoutMaker();

        void set_skybounds();
        void read_coadd_ascii();
        void append_cutout_info(const Position *pos,
                                const Transformation *trans,
                                long index);
        void get_file_cutout_info(int ifile, const Bounds *bounds);
        void get_cutout_info();

        void offset_positions();
        void load_coadd_trans();

        void set_start_rows_and_npix();

        template <typename T>
            void write_cutout(const Image<T> *image,
                              long iobj, long icut,
                              cutout_t scale);

        void write_fake_coadd_seg();
        void write_cutouts_from_coadd(enum cutout_type cut_type);
        void write_cutouts_from_file(long ifile);
        void write_weight_cutouts_from_file(int ifile);
        void write_seg_cutouts_from_file(int ifile);
        void write_mosaic(enum cutout_type cut_type);

        void open_fits(); // after catalogs written
        void close_fits();

        void push_image_file(string fname);
        void push_wcs_file(string fname);
        void push_sky_file(string fname);
        void push_seg_file(string fname);
        void load_data();
        void make_scale_factors();

        string setup_se_file(long ifile,
                             enum cutout_type cut_type,
                             int *hdu);
        bool load_file_cutouts(int ifile, const Bounds *bounds);

        void write_metadata(CCfits::FITS *fits);
        void write_image_info(CCfits::FITS *fits);
        void write_catalog();

        void print_area_stats();

        //int count_noflags();
        long get_max_cutouts();

        void print_file_progress(long ifile, const char *type, string fname);

    private:

        vector<long> id;
        vector<Position> pos;
        vector<Position> skypos;
        //vector<long> flags;

        vector<long> cutout_count;
        vector<long> box_size;

        long total_pixels;

        vector<vector<long> > start_row;

        vector<vector<Position> > cutout_pos; // positions in the cutout

        vector<vector<long> > orig_file_id; // id of SE file
        vector<vector<Position> > orig_pos; // position in SE image
        vector<vector<long> > orig_start_col; // box corner
        vector<vector<long> > orig_start_row; // box corner

        vector<vector<double> > dudrow;
        vector<vector<double> > dudcol;
        vector<vector<double> > dvdrow;
        vector<vector<double> > dvdcol;

        Bounds skybounds; 
        string imlist_path;

        vector<long> image_id_list;
        vector<long> image_flags_list;
        vector<string> image_file_list;
        vector<string> wcs_file_list; // only used if wcs is in a separate file
        vector<string> sky_file_list;
        vector<string> seg_file_list;
        vector<cutout_t> scale_list;
        vector<cutout_t> magzp_list;
        cutout_t magzp_ref;

        long max_image_filename_len;
        long max_wcs_filename_len;
        long max_sky_filename_len;
        long max_seg_filename_len;

        Transformation coadd_trans;

        ConfigFile params;

        fitsfile *fits;

        long nobj;
        long ncutout;

        string cutout_filename;
};


CutoutMaker::CutoutMaker(const ConfigFile *input_params) :

        params(*input_params),
        fits(NULL),
        ncutout(0)

{
    this->read_coadd_ascii();

    this->offset_positions();
    this->load_coadd_trans();

    this->imlist_path=this->params["coadd_srclist"];
    this->cutout_filename=this->params["cutout_file"];

    this->load_data();
    this->make_scale_factors();

    this->cutout_count.resize(this->nobj,0);
    this->cutout_pos.resize(this->nobj);

    this->orig_file_id.resize(this->nobj);
    this->orig_pos.resize(this->nobj);
    this->orig_start_col.resize(this->nobj);
    this->orig_start_row.resize(this->nobj);

    this->dudrow.resize(this->nobj);
    this->dudcol.resize(this->nobj);
    this->dvdrow.resize(this->nobj);
    this->dvdcol.resize(this->nobj);
}

CutoutMaker::~CutoutMaker()
{
    this->close_fits();
}

void CutoutMaker::set_skybounds()
{
    long nrows=this->skypos.size();

    for(long i=0; i<nrows; ++i) {
        //if (!this->flags[i]) {
        //    this->skybounds += this->skypos[i];
        //}
        this->skybounds += this->skypos[i];
    }
}

static long make_box_size_even(long box_size)
{
    if ((box_size % 2) != 0) {
        box_size += 1;
    }
    return box_size;
}
void CutoutMaker::read_coadd_ascii()
{
    string fname=this->params["cat_file"];
    if (!DoesFileExist(fname)) {
        cerr<<"File does not exist: "<<fname<<"\n";
        exit(1);
    }

    std::ifstream cat(fname.c_str(), std::ios::in);

    long id;
    double ra,dec,row,col;
    long bsize;
    this->nobj=0;

    while (cat >> id >> ra >> dec >> row >> col >> bsize) {

        bsize = make_box_size_even(bsize);

        this->id.push_back(id);

        Position sky_pos(ra,dec);
        Position image_pos(col,row);

        sky_pos *= 3600.; // deg->arcsec

        this->pos.push_back(image_pos);
        this->skypos.push_back(sky_pos);
        this->box_size.push_back(bsize);

        this->nobj++;
    }

    this->set_skybounds();
}

void CutoutMaker::offset_positions()
{
    std::complex<double> off1(1,1);
    for (int i=0; i<this->nobj; i++) {
        this->pos[i] -= off1;
    }
}

void CutoutMaker::load_coadd_trans()
{
    string fname=this->params["coadd_wcs_file"];
    int hdu = this->params["coadd_wcs_hdu"];
    if (fname != this->params["coadd_file"]) {
        cerr<<"Reading coadd wcs: "<<fname<<" "<<hdu<<"\n";
    }
    this->coadd_trans.readWCS(fname,hdu);
}

static fitsfile *_close_fits(fitsfile *fits)
{
    if (!fits) {
        return NULL;
    }
    int fitserr=0;
    fits_close_file(fits, &fitserr);
    if (fitserr != 0) {
        fits_report_error(stderr,fitserr);
        throw WriteException("Error closing fits file");
    }
    return NULL;
}
//READWRITE
//READONLY
static fitsfile *_open_fits(string filename, int mode)
{
    fitsfile *fits=NULL;
    int fitserr=0;

    fits_open_file(&fits,filename.c_str(),mode,&fitserr);
    if (fitserr != 0) {
        fits_report_error(stderr,fitserr);
        throw WriteException(
                "Error opening fits file " + filename);
    }

    return fits;
}

static void get_image_shape(string filename, int hdu, long *ncol, long *nrow)
{
    int fitserr=0;

    fitsfile *fits=_open_fits(filename,READONLY);

    fits_movabs_hdu(fits,hdu,0,&fitserr);
    if (fitserr != 0) {
        fits_report_error(stderr,fitserr);
        std::stringstream ss;
        ss<<"Error moving to hdu "<<hdu<<" in file "<<filename;
        throw ReadException(ss.str());
    }

    int ndim=0;
    fits_get_img_dim(fits, &ndim, &fitserr);
    if (fitserr != 0) {
        fits_report_error(stderr,fitserr);
        throw ReadException("Error getting ndims from file "+filename);
    }

    long dims[ndim];
    fits_get_img_size(fits, ndim, dims, &fitserr);
    if (fitserr != 0) {
        fits_report_error(stderr,fitserr);
        throw ReadException("Error getting dims from file "+filename);
    }

    fits=_close_fits(fits);

    *ncol=dims[0];
    *nrow=dims[1];
}

// not efficient! only good if you just need a single keyword
static double get_dbl_keyword(string filename, int hdu, const char *keyname)
{
    int status=0;

    fitsfile *fits=_open_fits(filename,READONLY);

    fits_movabs_hdu(fits,hdu,0,&status);
    if (status != 0) {
        fits_report_error(stderr,status);
        std::stringstream ss;
        ss<<"Error moving to hdu "<<hdu<<" in file "<<filename;
        throw ReadException(ss.str());
    }

    double val=0;
    if (fits_read_key_dbl(fits,(char*)keyname,&val,NULL, &status)) {
        fits_report_error(stderr,status);
        std::stringstream ss;
        ss<<"Error reading header keyword "<<keyname<<" in file "<<filename;
        throw ReadException(ss.str());
    }

    fits=_close_fits(fits);

    return val;
}





#if 0
static fitsfile *_open_new_fits(string filename)
{
    fitsfile *fits=NULL;
    int fitserr=0;

    fits_create_file(&fits,("!"+filename).c_str(),&fitserr);
    if (fitserr != 0) {
        fits_report_error(stderr,fitserr);
        throw WriteException(
                "Error creating fits file " + filename);
    }

    return fits;
}
#endif

void CutoutMaker::close_fits()
{
    if (this->fits) {
        this->fits=_close_fits(this->fits);
    }
}

// at least size 2 because CCfits requires
// length 2 for array columns
long CutoutMaker::get_max_cutouts()
{
    long max_cutouts=2;
    for (long i=0; i<this->nobj; i++) {
        long count=this->orig_file_id[i].size();
        if (count > max_cutouts) {
            max_cutouts=count;
        }
    }
    return max_cutouts;
}
template <typename T>
void to_valarray_vec(const vector<vector<T> > *vvec,
                     vector<std::valarray<T> > *vval,
                     long arrsize)
{
    long n=vvec->size();
    vval->resize(n);
    for (long i=0; i<n; i++) {
        long s=vvec->at(i).size();
        if (s > arrsize) {
            std::stringstream ss;
            ss<<"number of cutouts "<<s<<" exceeded maximum "<<arrsize;
            throw std::runtime_error(ss.str());
        }

        vval->at(i).resize(arrsize,DEFVAL);
        for (long j=0; j<s; j++) {
            vval->at(i)[j] = vvec->at(i)[j];
        }
    }
}


static void vposition_to_xy(const vector<vector<Position> > *pvec,
                            vector<std::valarray<double> > *xvec,
                            vector<std::valarray<double> > *yvec,
                            long arrsize)
{
    long n=pvec->size();
    xvec->resize(n);
    yvec->resize(n);

    for (long i=0; i<n; i++) {
        long s=pvec->at(i).size();
        if (s > arrsize) {
            std::stringstream ss;
            ss<<"number of cutouts "<<s<<" exceeded maximum "<<arrsize;
            throw std::runtime_error(ss.str());
        }

        xvec->at(i).resize(arrsize,DEFVAL);
        yvec->at(i).resize(arrsize,DEFVAL);
        for (long j=0; j<s; j++) {
            xvec->at(i)[j] = pvec->at(i)[j].getX();
            yvec->at(i)[j] = pvec->at(i)[j].getY();
        }
    }
}
string get_desdata() {
    string desdata("undefined");
    const char *desdata_c=getenv("DESDATA");
    if (desdata_c) {
        desdata=desdata_c;
        if (desdata.size() == 0) {
            desdata="undefined";
        }
    }

    return desdata;
}

// write all params to an extension                                 
void CutoutMaker::write_metadata(CCfits::FITS *fits)
{
    std::map<std::string,ConvertibleString>::const_iterator iter;
    std::stringstream ss;
    string fmt;

    vector<string> col_names_req;
    vector<string> col_fmts_req;
    vector<string> col_names_add;
    vector<string> col_fmts_add;

    vector<string> col_names;
    vector<string> col_fmts;

    // 8 byte integer
    const char *kfmt="K";

    col_names_req.push_back("magzp_ref");
    col_fmts_req.push_back("E");

    col_names_req.push_back("se_hdu");
    col_fmts_req.push_back(kfmt);
    col_names_req.push_back("se_wt_hdu");
    col_fmts_req.push_back(kfmt);
    col_names_req.push_back("se_badpix_hdu");
    col_fmts_req.push_back(kfmt);
    col_names_req.push_back("sky_hdu");
    col_fmts_req.push_back(kfmt);
    col_names_req.push_back("seg_hdu");
    col_fmts_req.push_back(kfmt);
    col_names_req.push_back("coadd_hdu");
    col_fmts_req.push_back(kfmt);
    col_names_req.push_back("coadd_wt_hdu");
    col_fmts_req.push_back(kfmt);
    col_names_req.push_back("coadd_seg_hdu");
    col_fmts_req.push_back(kfmt);
    col_names_req.push_back("fake_coadd_seg");
    col_fmts_req.push_back(kfmt);


    ss.str("");

    vector<string> desdata_vec(1);
    desdata_vec[0]=get_desdata();
    ss<<desdata_vec[0].size()<<"A";
    fmt=ss.str();
    col_names_req.push_back("DESDATA");
    col_fmts_req.push_back(fmt);

    col_names=col_names_req;
    col_fmts = col_fmts_req;

    // add non-required keywords
    for (iter=this->params.begin(); iter!= this->params.end(); iter++) {
        string key=iter->first;
        if ( std::find(col_names.begin(), col_names.end(), key)==col_names.end()) {
            string val=iter->second;
            ss.str("");
            ss<<val.size()<<"A";
            fmt=ss.str();

            col_names.push_back(key);
            col_fmts.push_back(fmt);
        }
    }

    vector<string> col_units(col_fmts.size());

    long nrows=1;
    CCfits::Table* table = fits->addTable("metadata",nrows,
                                          col_names,col_fmts,col_units);

    long firstrow=1;

    vector<double> dvec(1);
    vector<long> ivec(1);

    dvec[0] = this->magzp_ref;
    table->column("magzp_ref").write(dvec,firstrow);

    ivec[0] = this->params["se_hdu"];
    table->column("se_hdu").write(ivec,firstrow);
    ivec[0] = this->params["se_wt_hdu"];
    table->column("se_wt_hdu").write(ivec,firstrow);
    ivec[0] = this->params["se_badpix_hdu"];
    table->column("se_badpix_hdu").write(ivec,firstrow);
    ivec[0] = this->params["sky_hdu"];
    table->column("sky_hdu").write(ivec,firstrow);
    ivec[0] = this->params["seg_hdu"];
    table->column("seg_hdu").write(ivec,firstrow);
    ivec[0] = this->params["coadd_hdu"];
    table->column("coadd_hdu").write(ivec,firstrow);
    ivec[0] = this->params["coadd_wt_hdu"];
    table->column("coadd_wt_hdu").write(ivec,firstrow);
    ivec[0] = this->params["fake_coadd_seg"];
    table->column("fake_coadd_seg").write(ivec,firstrow);

    table->column("DESDATA").write(desdata_vec,firstrow);

    // add not-required keywords
    vector<string> vval(1);  // CCfits really wants you to use a vector
    for (iter=this->params.begin(); iter!= this->params.end(); iter++) {
        string key=iter->first;
        if ( std::find(col_names_req.begin(), col_names_req.end(), key)==col_names_req.end()) {
            cerr<<"        writing key: "<<key<<"\n";
            vval[0]=iter->second;
            table->column(key).write(vval,firstrow);
        }
    }

}

void CutoutMaker::write_image_info(CCfits::FITS *fits)
{
    long nfiles=this->image_file_list.size();

    long ncols=8;
    vector<string> col_names(ncols);
    vector<string> col_fmts(ncols);
    vector<string> col_units(ncols);

    col_names[0] = "image_id";
    col_fmts[0] = "K"; // 8-byte integer

    col_names[1] = "image_flags";
    col_fmts[1] = "K"; // 8-byte integer

    std::stringstream ss;
    ss<<this->max_image_filename_len<<"A";
    col_names[2] = "image_path";
    col_fmts[2] = ss.str();

    ss.str("");
    ss<<this->max_wcs_filename_len<<"A";
    col_names[3] = "wcs_path";
    col_fmts[3] = ss.str();

    ss.str("");
    ss<<this->max_sky_filename_len<<"A";
    col_names[4] = "sky_path";
    col_fmts[4] = ss.str();

    ss.str("");
    ss<<this->max_seg_filename_len<<"A";
    col_names[5] = "seg_path";
    col_fmts[5] = ss.str();

    col_names[6] = "magzp";
    col_fmts[6] = "E";

    col_names[7] = "scale";
    col_fmts[7] = "E";

    //long max_cutouts=this->get_max_cutouts();

    CCfits::Table* table = fits->addTable("image_info",nfiles,
                                          col_names,col_fmts,col_units);


    long firstrow=1;
    table->column("image_id").write(this->image_id_list,firstrow);
    table->column("image_flags").write(this->image_flags_list,firstrow);
    table->column("image_path").write(this->image_file_list,firstrow);
    table->column("wcs_path").write(this->wcs_file_list,firstrow);
    table->column("sky_path").write(this->sky_file_list,firstrow);
    table->column("seg_path").write(this->seg_file_list,firstrow);
    table->column("magzp").write(this->magzp_list,firstrow);
    table->column("scale").write(this->scale_list,firstrow);

}

void make_number_vec(vector<long> *vec, long n)
{
    vec->resize(n);
    for (long i=0; i<n; i++) {
        vec->at(i) = i+1;
    }
}
// write catalog, file table, and metadata table
void CutoutMaker::write_catalog()
{
    cerr<<"opening output file: "<<this->cutout_filename<<"\n";
    CCfits::FITS fits("!"+this->cutout_filename, CCfits::Write);

    cerr<<"  writing catalog\n";
    long max_cutouts=this->get_max_cutouts();

    vector<long> number_vec;
    make_number_vec(&number_vec, this->box_size.size());

    std::stringstream ss;
    ss<<max_cutouts<<"K";
    string kfmt = ss.str();
    ss.str("");
    ss<<max_cutouts<<"D";
    string dfmt = ss.str();


    vector<string> col_names;
    vector<string> col_fmts;

    col_names.push_back("id");
    col_fmts.push_back("1K"); // 8-byte integer

    col_names.push_back("number");
    col_fmts.push_back("1K"); // 8-byte integer

    col_names.push_back("ncutout");
    col_fmts.push_back("1K");

    col_names.push_back("box_size");
    col_fmts.push_back("1K");

    col_names.push_back("file_id");
    col_fmts.push_back(kfmt);

    col_names.push_back("start_row");
    col_fmts.push_back(kfmt);

    col_names.push_back("orig_row");
    col_fmts.push_back(dfmt);
    col_names.push_back("orig_col");
    col_fmts.push_back(dfmt);


    col_names.push_back("orig_start_row");
    col_fmts.push_back(kfmt);
    col_names.push_back("orig_start_col");
    col_fmts.push_back(kfmt);

    col_names.push_back("cutout_row");
    col_fmts.push_back(dfmt);
    col_names.push_back("cutout_col");
    col_fmts.push_back(dfmt);

    col_names.push_back("dudrow");
    col_fmts.push_back(dfmt);
    col_names.push_back("dudcol");
    col_fmts.push_back(dfmt);
    col_names.push_back("dvdrow");
    col_fmts.push_back(dfmt);
    col_names.push_back("dvdcol");
    col_fmts.push_back(dfmt);


    // dummy
    vector<string> col_units(col_fmts.size());

    CCfits::Table* table = fits.addTable("object_data",
                                         this->nobj,
                                         col_names,
                                         col_fmts,
                                         col_units);


    //this->write_image_info(&fits);
    //this->write_metadata(&fits);

    long firstrow=1;
    table->column("id").write(this->id,firstrow);
    table->column("number").write(number_vec,firstrow);
    table->column("ncutout").write(this->cutout_count,firstrow);
    table->column("box_size").write(this->box_size,firstrow);

    vector<std::valarray<long> > orig_file_id;
    to_valarray_vec(&this->orig_file_id, &orig_file_id, max_cutouts);
    table->column("file_id").writeArrays(orig_file_id,firstrow);

    vector<std::valarray<long> > start_row;
    to_valarray_vec(&this->start_row, &start_row, max_cutouts);
    table->column("start_row").writeArrays(start_row,firstrow);

    vector<std::valarray<double> > x;
    vector<std::valarray<double> > y;

    vposition_to_xy(&this->orig_pos, &x, &y, max_cutouts);
    table->column("orig_row").writeArrays(y,firstrow);
    table->column("orig_col").writeArrays(x,firstrow);

    // re-using start_row storage
    to_valarray_vec(&this->orig_start_row, &start_row, max_cutouts);
    table->column("orig_start_row").writeArrays(start_row,firstrow);

    to_valarray_vec(&this->orig_start_col, &start_row, max_cutouts);
    table->column("orig_start_col").writeArrays(start_row,firstrow);

    vposition_to_xy(&this->cutout_pos, &x, &y, max_cutouts);
    table->column("cutout_row").writeArrays(y,firstrow);
    table->column("cutout_col").writeArrays(x,firstrow);


    vector<std::valarray<double> > deriv;

    to_valarray_vec(&this->dudrow, &deriv, max_cutouts);
    table->column("dudrow").writeArrays(deriv,firstrow);

    to_valarray_vec(&this->dudcol, &deriv, max_cutouts);
    table->column("dudcol").writeArrays(deriv,firstrow);

    to_valarray_vec(&this->dvdrow, &deriv, max_cutouts);
    table->column("dvdrow").writeArrays(deriv,firstrow);

    to_valarray_vec(&this->dvdcol, &deriv, max_cutouts);
    table->column("dvdcol").writeArrays(deriv,firstrow);

    cerr<<"  writing image info\n";
    this->write_image_info(&fits);
    cerr<<"  writing meta data\n";
    this->write_metadata(&fits);

    // RAII will write and flush the data
}


template <typename T> 
inline int getBitPix() { return 0; }
template <> 
inline int getBitPix<double>() { return DOUBLE_IMG; }
template <> 
inline int getBitPix<float>() { return FLOAT_IMG; }
template <> 
inline int getBitPix<int>() { return LONG_IMG; }


template <typename T>
void create_mosaic(fitsfile* fits, long total_pixels,
                   string extname)
{
    int fitserr=0;

    int n_axes = 1;
    long dims[1] = {total_pixels};
    cerr<<"        mosaic size: "<<dims[0]<<"\n";

    Assert(getBitPix<T>());
    int bitpix = getBitPix<T>();
    fits_create_img(fits, bitpix, n_axes, dims, &fitserr);
    if (fitserr != 0) {
        fits_report_error(stderr,fitserr);
        throw WriteException("Error creating mosaic");
    }

    fits_write_key_str(fits, 
                       "extname", 
                       (char*)extname.c_str(),
                       NULL,
                       &fitserr);
    if (fitserr != 0) {
        fits_report_error(stderr,fitserr);
        throw WriteException("Error writing extname");
    }
}

void CutoutMaker::open_fits()
{
    this->fits=_open_fits(this->cutout_filename,READWRITE);
}


void CutoutMaker::print_area_stats()
{
    cerr<<"Total bounds are "<<this->skybounds<<std::endl;

    double area = this->skybounds.getArea();
    area /= 3600.;
    cerr<<"Raw area = "<<area<<" square arcmin\n";

    double dec = this->skybounds.getCenter().getY();
    double cosdec = cos(dec / ARCSEC_PER_RAD);
    area *= cosdec;
    cerr<<"Corrected area = "<<area<<" square arcmin";

    area /= 3600.;
    cerr<<" = "<<area<<" square degrees\n";
}


void CutoutMaker::push_image_file(string fname)
{
    this->image_file_list.push_back(fname);

    if (int(fname.size()) > this->max_image_filename_len) {
        this->max_image_filename_len=fname.size();
    }
}
void CutoutMaker::push_wcs_file(string fname)
{
    this->wcs_file_list.push_back(fname);

    if (int(fname.size()) > this->max_wcs_filename_len) {
        this->max_wcs_filename_len=fname.size();
    }
}

void CutoutMaker::push_sky_file(string fname)
{
    this->sky_file_list.push_back(fname);

    if (int(fname.size()) > this->max_sky_filename_len) {
        this->max_sky_filename_len=fname.size();
    }
}
void CutoutMaker::push_seg_file(string fname)
{
    this->seg_file_list.push_back(fname);

    if (int(fname.size()) > this->max_seg_filename_len) {
        this->max_seg_filename_len=fname.size();
    }
}



template <typename T> T get_median(vector<T> &vec)
{
    vector<T> vcpy = vec;

    std::sort( vcpy.begin(), vcpy.end() );

    long i=vcpy.size()/2;
    return vcpy[i];
}

// make scales to put all the images on the same zero point
void CutoutMaker::make_scale_factors()
{
    //this->scale_list.push_back(1.0);  // the coadd is not scaled
    //this->magzp_ref = get_median(this->magzp_list);

    this->magzp_ref = this->params["magzp_ref"];
    cerr<<"ref zero point: "<<this->magzp_ref<<"\n";
    for (long i=0; i<long(this->magzp_list.size()); i++) {
        cutout_t magzp = this->magzp_list[i];
        cutout_t scale = pow(10.0, 0.4*(this->magzp_ref-magzp) );
        this->scale_list.push_back(scale);
    }

}
void CutoutMaker::load_data()
{
    if (!DoesFileExist(this->imlist_path)) {
        cerr<<"file not found: "<<this->imlist_path<<"\n";
        throw FileNotFoundException(this->imlist_path);
    }

    dbg<<"Opening coadd srclist\n";
    std::ifstream sedatafile(this->imlist_path.c_str(), std::ios::in);
    if (!sedatafile) {
        string err="Unable to open source list file " + this->imlist_path;
        cerr<<err<<std::endl;
        throw ReadException(err);
    }


    this->max_image_filename_len=0;
    this->max_wcs_filename_len=0;
    this->max_sky_filename_len=0;
    this->max_seg_filename_len=0;

    // first push the coadd, always first
    // note sky file for coadd is itself


    double coadd_magzp = get_dbl_keyword(this->params["coadd_file"],
                                         this->params["coadd_hdu"],
                                         "SEXMGZPT");

    this->image_id_list.push_back(this->params["coadd_image_id"]);
    this->image_flags_list.push_back(0);
    this->magzp_list.push_back(coadd_magzp);
    this->push_image_file(this->params["coadd_file"]);
    this->push_sky_file(this->params["coadd_file"]);
    this->push_seg_file(this->params["coaddseg_file"]);

    // only for bookkeeping right now
    this->push_wcs_file(this->params["coadd_wcs_file"]);

    string image_path, sky_path, seg_path, wcs_path;
    long image_id, image_flags;
    double magzp;

    if (this->params["use_alt_wcs"]) {
        cerr<<"    reading with alt wcs\n";
        while (sedatafile >> image_id >> image_flags >> image_path >> sky_path >> seg_path >> magzp >> wcs_path) {
            this->image_id_list.push_back(image_id);
            this->image_flags_list.push_back(image_flags);
            this->push_image_file(image_path);
            this->push_sky_file(sky_path);
            this->push_seg_file(seg_path);
            this->magzp_list.push_back(magzp);
            this->push_wcs_file(wcs_path);
        }

    } else {
        while (sedatafile >> image_id >> image_flags >> image_path >> sky_path >> seg_path >> magzp) {
            this->image_id_list.push_back(image_id);
            this->image_flags_list.push_back(image_flags);
            this->push_image_file(image_path);
            this->push_sky_file(sky_path);
            this->push_seg_file(seg_path);
            this->magzp_list.push_back(magzp);
            this->push_wcs_file(image_path);
        }
    }
    cerr<<"read "
        <<this->image_file_list.size()
        <<" se images"<<std::endl;
}

// use the rough transform as input to the nonlinear finder
// If it doesn't work, use the rough one
static int get_image_pos(const Transformation *trans,
                         const Transformation *inv_trans,
                         const Position *skypos,
                         Position *pos)
{
    int good=1;
    inv_trans->transform(*skypos,*pos);

    Position pos_tmp = (*pos);
    if (!trans->inverseTransform(*skypos, pos_tmp) ) {
        dbg << "InverseTransform failed for position "<<*skypos<<".\n";
        dbg << "Initial guess was "<<*pos<<".\n";
        good=0;
    } else {
        (*pos) = pos_tmp;
    }

    return good;
}

// position is position in the original image
void CutoutMaker::append_cutout_info(const Position *pos,
                                     const Transformation *trans,
                                     long index)
{
    double dudx=0, dudy=0, dvdx=0, dvdy=0;
    trans->getDistortion(*pos, dudx, dudy, dvdx, dvdy); // by ref

    double px=int(pos->getX());
    double py=int(pos->getY());

    // note the box size is always even
    int offset = this->box_size[index]/2;
    int startx = px-offset;
    int starty = py-offset;

    Position cpos=(*pos);
    cpos -= Position(startx,starty);

    this->orig_pos[index].push_back((*pos));
    this->cutout_pos[index].push_back(cpos);

    this->orig_start_col[index].push_back(startx);
    this->orig_start_row[index].push_back(starty);

    this->dudcol[index].push_back(dudx);
    this->dudrow[index].push_back(dudy);
    this->dvdcol[index].push_back(dvdx);
    this->dvdrow[index].push_back(dvdy);

    this->cutout_count[index]++;
    this->ncutout++;
}

void CutoutMaker::print_file_progress(long ifile, const char *type, string fname)
{
    cerr<<"    "
        <<(ifile+1)
        <<"/"
        <<this->image_file_list.size()
        <<" "
        <<type
        <<" "
        <<fname<<"\n";
}

string CutoutMaker::setup_se_file(long ifile,
                                  enum cutout_type cut_type,
                                  int *hdu)
{
    string fname;
    
    if (cut_type==CUTOUT_SKY) {
        fname = this->sky_file_list[ifile];
        (*hdu) = this->params["sky_hdu"];
    } else if (cut_type==CUTOUT_SEG) {
        fname = this->seg_file_list[ifile];
        (*hdu) = this->params["seg_hdu"];
    } else if (cut_type==CUTOUT_WEIGHT) {
        fname = this->image_file_list[ifile];
        (*hdu) = this->params["se_wt_hdu"];
    } else if (cut_type==CUTOUT_BADPIX) {
        fname = this->image_file_list[ifile];
        (*hdu) = this->params["se_badpix_hdu"];
    } else {
        fname = this->image_file_list[ifile];
        (*hdu) = this->params["se_hdu"];
    }

    return fname;
}




void CutoutMaker::set_start_rows_and_npix()
{
    long current_row=0;
    this->start_row.resize(this->nobj);
    this->total_pixels=0;

    for (long i=0; i<this->nobj; i++) {
        long ncut=this->cutout_pos[i].size();
        if (ncut > 0) {

            long bsize=this->box_size[i];
            long npix_per_cutout = bsize*bsize;

            this->start_row[i].resize(ncut);

            for (long j=0; j<ncut; j++) {

                this->start_row[i][j] = current_row;

                current_row += npix_per_cutout;
                this->total_pixels += npix_per_cutout;
            }
        }
    }
}



/*
   Scale is only for the cutouts and weight images. Must be applied after sky
   subtraction.  Send as negative to avoid application.

   badpix can be null;  can send with weight image to set zero when
   there is a bad pixel
*/

template <typename T>
void CutoutMaker::write_cutout(const Image<T> *image,
                               long iobj, long icut,
                               cutout_t scale)
{
    long xmax=image->getXMax();
    long ymax=image->getYMax();

    long startx=this->orig_start_col[iobj][icut];
    long starty=this->orig_start_row[iobj][icut];

    long box_size=this->box_size[iobj];

    cutout_t val=0;
    Image<cutout_t> tim(box_size,box_size);
    for (long ix=0, x=startx; ix<box_size; ix++, x++) {
        if (x < 0 || x >= xmax ) continue;

        for (long iy=0, y=starty; iy<box_size; iy++, y++) {
            if (y < 0 || y >= ymax ) continue;

            val = (*image)(x,y);
            
            if (scale > 0) {
                val *= scale;
            }
            tim(ix,iy) = val;
        }
    }

    long start_row=1+this->start_row[iobj][icut];
    tim.write_sub_flat(this->fits, start_row);
}

void CutoutMaker::write_fake_coadd_seg() 
{
    
    //Image<cutout_t> *sky_null=NULL;
    //Image<badpix_t> *badpix_null=NULL;
    cutout_t scale=-9999; // negative means don't apply it

    long ncol=0, nrow=0;
    get_image_shape(this->params["coadd_file"],
                    this->params["coadd_hdu"],
                    &ncol,
                    &nrow);

    Image<seg_t> seg_image(ncol,nrow);
    for (long i=0; i<ncol; i++) {
        for (long j=0; j<ncol; j++) {
            seg_image(i,j) = DEFVAL;
        }
    }

    for (long i=0; i<this->nobj; ++i) { 
        long ncut=this->cutout_pos[i].size();
        for (long j=0; j<ncut; j++) {
            if (this->orig_file_id[i][j] == 0) { 
                this->write_cutout(&seg_image, i, j, scale);
            }
        }
    }

}


// here we use the fact that the coadd is always in position ifile==0
// note either image or seg_image is written, depending on cut_type
void CutoutMaker::write_cutouts_from_coadd(enum cutout_type cut_type)
{
    Image<cutout_t> image;
    Image<seg_t> seg_image;
    cutout_t scale=-9999; // negative means don't apply it
    int hdu=0;

    const char *print_type=NULL;
    string fname;
    if (cut_type==CUTOUT_SEG && this->params["fake_coadd_seg"]) {
        cerr<<"    faking coadd seg image\n";
        this->write_fake_coadd_seg();
        return;
    }

    if (cut_type==CUTOUT_SEG) {
        cerr<<"    loading coadd seg image\n";
        hdu = this->params["coadd_seg_hdu"];
        print_type="seg";
        fname=this->params["coaddseg_file"];
        seg_image.load(fname, hdu);
    } else if (cut_type==CUTOUT_WEIGHT) {
        cerr<<"    loading coadd weightimage\n";
        hdu = this->params["coadd_wt_hdu"];
        print_type="weight";
        fname=this->params["coadd_file"];
        image.load(fname, hdu);
        scale = this->scale_list[0];
        scale = 1.0/(scale*scale);
    } else {
        cerr<<"    loading coadd image\n";
        hdu = this->params["coadd_hdu"];
        print_type="image";
        fname=this->params["coadd_file"];
        image.load(fname, hdu);
        scale = this->scale_list[0];
    }

    this->print_file_progress(0,print_type,fname);

    for (long i=0; i<this->nobj; ++i) {
        long ncut=this->cutout_pos[i].size();
        for (long j=0; j<ncut; j++) {
            if (this->orig_file_id[i][j] == 0) {
                if (cut_type==CUTOUT_SEG) {
                    this->write_cutout(&seg_image, i, j, scale);
                } else {
                    this->write_cutout(&image, i, j, scale);
                }
            }
        }
    }
}

// subtract the sky, setting to zero for bad pixels
void subtract_sky_image_with_bpm(Image<cutout_t> *image,
                                 const Image<cutout_t> *sky,
                                 const Image<badpix_t> *bpm)
{
    long xmax=image->getXMax();
    long ymax=image->getYMax();

    for (long ix=0; ix<xmax; ix++) {
        for (long iy=0; iy<ymax; iy++) {
            if ( (*bpm)(ix,iy) != 0) {
                (*image)(ix,iy) = 0.0;
            } else {
                (*image)(ix,iy) -= (*sky)(ix,iy);
            }
        }
    }
}

// these are the SE image cutouts
void CutoutMaker::write_cutouts_from_file(long ifile)
{

    // ifile==0 means coadd
    if (ifile==0) {
        return;
    }

    int hdu=0, sky_hdu=0, badpix_hdu=0;
    string filename=this->setup_se_file(ifile, CUTOUT_IMAGE, &hdu);
    string sky_filename=this->setup_se_file(ifile, CUTOUT_SKY, &sky_hdu);
    string badpix_filename=this->setup_se_file(ifile, 
            CUTOUT_BADPIX, &badpix_hdu);

    this->print_file_progress(ifile,"image",filename);

    Image<cutout_t> image(filename, hdu);
    Image<cutout_t> sky_image(sky_filename, sky_hdu);
    Image<badpix_t> bpm(badpix_filename, badpix_hdu);

    subtract_sky_image_with_bpm(&image, &sky_image, &bpm);

    cutout_t scale = this->scale_list[ifile];

    for (long i=0; i<this->nobj; ++i) {
        long ncut=this->cutout_pos[i].size();
        for (long j=0; j<ncut; j++) {
            if (this->orig_file_id[i][j] == ifile) {
                this->write_cutout(&image, i, j, scale);
            }
        }
    }
}

void set_weight_zero_at_bad_pix(Image<cutout_t> *wt,
                                Image<badpix_t> *bpm)
{
    long xmax=wt->getXMax();
    long ymax=wt->getYMax();

    for (long ix=0; ix<xmax; ix++) {
        for (long iy=0; iy<ymax; iy++) {
            if ( (*bpm)(ix,iy) != 0) {
                (*wt)(ix,iy) = 0.0;
            }
        }
    }
}
// For weight image we scale by 1/scale^2
void CutoutMaker::write_weight_cutouts_from_file(int ifile)
{
    // ifile==0 means coadd
    if (ifile==0) {
        return;
    }

    int hdu=0, badpix_hdu=0;
    string filename=this->setup_se_file(ifile, CUTOUT_WEIGHT, &hdu);
    string badpix_filename=this->setup_se_file(ifile, 
            CUTOUT_BADPIX, &badpix_hdu);

    this->print_file_progress(ifile,"weight",filename);

    Image<cutout_t> wt(filename, hdu);
    Image<badpix_t> bpm(badpix_filename, badpix_hdu);

    set_weight_zero_at_bad_pix(&wt, &bpm);

    cutout_t scale = this->scale_list[ifile];
    cutout_t scale_inv2 = 1.0/(scale*scale);

    for (long i=0; i<this->nobj; ++i) {
        long ncut=this->cutout_pos[i].size();
        for (long j=0; j<ncut; j++) {
            if (this->orig_file_id[i][j] == ifile) {
                this->write_cutout(&wt, i, j, scale_inv2);
            }
        }
    }
}

void CutoutMaker::write_seg_cutouts_from_file(int ifile)
{
    // ifile==0 means coadd
    if (ifile==0) {
        return;
    }
    cutout_t scale=-9999; // negative means don't apply it

    int hdu=0;
    //int sky_hdu=0;
    string filename=this->setup_se_file(ifile, CUTOUT_SEG, &hdu);

    this->print_file_progress(ifile,"seg",filename);

    Image<seg_t> seg_image(filename, hdu);

    for (long i=0; i<this->nobj; ++i) {
        long ncut=this->cutout_pos[i].size();
        for (long j=0; j<ncut; j++) {
            if (this->orig_file_id[i][j] == ifile) {
                this->write_cutout(&seg_image, i, j, scale);
            }
        }
    }
}

// Write the idividual cutouts into the big mosaic layed out on disk.
void CutoutMaker::write_mosaic(enum cutout_type cut_type)
{
    this->open_fits();
    string extname;

    if (cut_type==CUTOUT_WEIGHT) {
        extname="weight_cutouts";
    } else if (cut_type==CUTOUT_SEG) {
        extname="seg_cutouts";
    } else {
        extname="image_cutouts";
    }

    if (cut_type==CUTOUT_SEG) {
        create_mosaic<seg_t>(this->fits, this->total_pixels,
                           extname);
    } else {
        create_mosaic<cutout_t>(this->fits, this->total_pixels,
                extname);
    }

    this->write_cutouts_from_coadd(cut_type);

    const long nfiles = this->image_file_list.size();
    for (long ifile=1; ifile<nfiles; ifile++) {
        if (cut_type==CUTOUT_IMAGE) {
            this->write_cutouts_from_file(ifile);
        } else if (cut_type==CUTOUT_WEIGHT) {
            this->write_weight_cutouts_from_file(ifile);
        } else if (cut_type==CUTOUT_SEG) {
            this->write_seg_cutouts_from_file(ifile);
        }
    }

    this->close_fits();
}



// the input bounds are on the sky
// can avoid reading the image by just grabbing the header keywords
void CutoutMaker::get_file_cutout_info(int ifile, const Bounds *bounds)
{

    int hdu=0;
    string filename=this->setup_se_file(ifile, CUTOUT_IMAGE, &hdu);
    this->print_file_progress(ifile,"info",filename);

    long ncol=0, nrow=0;
    get_image_shape(filename, hdu, &ncol, &nrow);

    string wcs_filename=this->wcs_file_list[ifile];
    Transformation trans;
    hdu = this->params["se_wcs_hdu"];
    if (wcs_filename != filename) {
        cerr<<"                wcs: "<<wcs_filename<<" "<<hdu<<"\n";
    }
    trans.readWCS(wcs_filename, hdu);
    Transformation inv_trans;

    // rough inverse trans, only valid on the bounds of this image
    // wcs uses 1-offset
    Bounds se_bounds(1, ncol, 1, nrow);
    Bounds inv_bounds = inv_trans.makeInverseOf(trans,se_bounds,4);

    if (!this->skybounds.intersects(inv_bounds)) {
        cerr<<"Skipping "<<filename<<" because inv_bounds doesn't intersect\n";
        return;
    }

    long nkeep=0, ngood=0, ntry=0;
    std::complex<double> off1(1,1);
    Position pos(0,0);

    for (long i=0; i<this->nobj; ++i) {
        //if (this->flags[i]) continue;
        if (!bounds->includes(this->skypos[i])) continue;
        if (!inv_bounds.includes(this->skypos[i])) continue;

        ntry++;
        ngood+=get_image_pos(&trans, &inv_trans, &this->skypos[i], &pos);

        if (!se_bounds.includes(pos)) continue;

        pos -= off1; // zero-offset

        // always push the coadd image on first
        if (this->orig_file_id[i].size()==0) {
            this->append_cutout_info(&this->pos[i],&this->coadd_trans,i);
            this->orig_file_id[i].push_back(0); // coadd is always at zero
        }
        this->append_cutout_info(&pos, &trans, i);
        this->orig_file_id[i].push_back(ifile);
        nkeep++;
    }

    cerr<<"        kept "<<nkeep<<"/"<<this->nobj
        <<" good trans "<<ngood<<"/"<<ntry<<"\n";
}


// loop through and find the cutout boundaries, etc.  note ifile=0 would be the
// coadd, but we deal with that specially
void CutoutMaker::get_cutout_info()
{

    const long nfiles = this->image_file_list.size();
    for (long ifile=1; ifile<nfiles; ifile++) {
        this->get_file_cutout_info(ifile,&this->skybounds);
    }
    cerr<<"Found "<<this->ncutout<<" cutouts\n";

    // determine where the objects will go in the mosaic
    this->set_start_rows_and_npix();
}

static void make_cutouts(ConfigFile *params) 
{

    CutoutMaker maker(params);

    maker.print_area_stats();
    maker.get_cutout_info();

    maker.write_catalog();
    cerr<<"  writing cutout mosaic\n";
    maker.write_mosaic(CUTOUT_IMAGE);
    cerr<<"  writing weight mosaic\n";
    maker.write_mosaic(CUTOUT_WEIGHT);
    cerr<<"  writing seg mosaic\n";
    maker.write_mosaic(CUTOUT_SEG);

    cerr<<"output is in : "<<params->get("cutout_file")<<"\n";
}


static void usage_and_exit()
{
    cerr<<"make-cutouts coadd_file= cat_file= coadd_srclist= cutout_file=\n";
    exit(1);
}

// these can be over-written in load_params
void set_default_params(ConfigFile *params)
{
    params->set("se_hdu", SE_HDU);
    params->set("se_wt_hdu", SE_WT_HDU);
    params->set("se_badpix_hdu", SE_BADPIX_HDU);
    params->set("sky_hdu", SKY_HDU);
    params->set("seg_hdu", SEG_HDU);
    params->set("coadd_hdu", COADD_HDU);
    params->set("coadd_wt_hdu", COADD_WT_HDU);
    params->set("coadd_seg_hdu", COADD_SEG_HDU);

    params->set("fake_coadd_seg", 0);

    // by default we grab the wcs from the same file and hdu as the image
    params->set("coadd_wcs_hdu",COADD_HDU);
    params->set("se_wcs_hdu",SE_HDU);

    params->set("use_alt_wcs",0);

    params->set("coadd_image_id",-9999);

}
// check params.  also set some default params
bool check_params(const ConfigFile *params)
{
    bool ok=true;

    if (!params->keyExists("cat_file")) {
        cerr<<"send cat_file=\n";
        ok=false;
    }
    if (!params->keyExists("coadd_file")) {
        cerr<<"send coadd_file=\n";
        ok=false;
    }
    if (!params->keyExists("coaddseg_file")) {
        cerr<<"send coaddseg_file=\n";
        ok=false;
    }
    if (!params->keyExists("coadd_srclist")) {
        cerr<<"send coadd_srclist=\n";
        ok=false;
    }
    if (!params->keyExists("cutout_file")) {
        cerr<<"send cutout_file=\n";
        ok=false;
    }
    if (!params->keyExists("magzp_ref")) {
        cerr<<"send magzp_ref=\n";
        ok=false;
    }


    return ok;
}

// values set in set_default_params can be over-written
void load_params(int argc, char **argv, ConfigFile *params)
{
    for(int k=1;k<argc;k++) {
        params->append(argv[k]);
    }

    // if coadd_wcs_file is not sent, default to the coadd file
    // note coadd_wcs_hdu already has a default set
    if (!params->keyExists("coadd_wcs_file")) {
        params->set("coadd_wcs_file", params->get("coadd_file"));
    }

}

int main(int argc, char **argv)
{
    ConfigFile params;
    set_default_params(&params);
    load_params(argc, argv, &params);

    if (!check_params(&params)) {
        usage_and_exit();
    }

    make_cutouts(&params);
    return 0;
}
