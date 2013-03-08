/*
   Author: Erin Sheldon, BNL.  This began as a hack on the multishear code.

   notes

       - memory usage is moderate, ~1G on the coadd I've tested.  When adding
       the weight images and seg images this will increase to <~3G.  Maybe the
       max for any tile would be 10G.  
       
       Could probably remove all the memory management machinery for typical
       coadds, which would simplify the code. On the other hand, we may want to
       run this on the supernova fields and eventually lsst.  Also, this code
       is reasonably fast as is...

   TODO
       - decide if using 1 or 0 offset, both in terms of locations in
         images and row numbers for our cutouts
       - add sky values in table and sky subtract the images
       - add wcs info to the second extension along with the filenames.
       - when seg images are available, add them
       - create matching output file names
*/
#include <valarray>
#include <sys/time.h>
#include <iostream>
#include <fstream>
#include <CCfits/CCfits>

#include "Image.h"
#include "Transformation.h"
#include "CoaddCatalog.h"
#include "BasicSetup.h"

const double ARCSEC_PER_RAD = 206264.806247;
const int MAX_HDU=999;
const int DEFVAL=-9999;

enum cutout_type {
    CUTOUT_IMAGE,
    CUTOUT_WEIGHT
};

// define just work on the first few, good for testing
#define MAXIMAGES 20

using std::auto_ptr;
using std::cerr;
using std::string;

class CutoutMaker
{
    public:
        // no need to construct coadd outside, should
        // do this internally
        CutoutMaker(const CoaddCatalog *coaddCat, 
                    const ConfigFile *input_params);

        ~CutoutMaker();

        void append_cutout_info(const Position *pos, int index);
        void get_file_cutout_info(int ifile, const Bounds *bounds);
        void get_cutout_info();

        void load_coadd_image();
        void set_box_size();

        void set_start_rows();
        void write_cutout(const Image<double> *image, int iobj, int icut);
        void write_mosaics_from_coadd(enum cutout_type cut_type);
        void write_mosaics_from_file(int ifile, enum cutout_type cut_type);
        void write_mosaic(enum cutout_type cut_type);

        void open_fits(); // after catalogs written
        void close_fits();

        void push_image_file(string fname);
        void load_image_list();

        string setup_se_file(int ifile);
        bool load_file_cutouts(int ifile, const Bounds *bounds);

        void write_catalog_filenames(CCfits::FITS *fits);
        void write_catalog();

        void print_area_stats();

        int get_box_size(int index);

        int count_noflags();
        int get_max_cutouts();

        void print_file_progress(int ifile);

    private:

        Image<double> coadd_image;
        //auto_ptr<Image<double> > coadd_weight_image;
        Image<double> coadd_weight_image;

        // for input

        std::vector<long> id;
        std::vector<Position> pos;
        std::vector<Position> skypos;
        std::vector<long> flags;

        std::vector<int> cutout_count;
        std::vector<int> box_sizes; // constant for now

        // these are zero-offset to be consistent with other
        // indices in this file.  Need to convert to 1-offset
        // when using cfitsio
        std::vector<std::vector<int> > start_row;

        std::vector<std::vector<Position> > cutout_pos; // positions in the cutout

        std::vector<std::vector<int> > orig_file_id; // id of SE file
        std::vector<std::vector<Position> > orig_pos; // position in SE image
        std::vector<std::vector<int> > orig_start_col; // box corner
        std::vector<std::vector<int> > orig_start_row; // box corner

        Bounds skybounds; 
        string imlist_path;

        std::vector<string> image_file_list;
        int max_image_filename_len;

        ConfigFile params;

        fitsfile *fits;

        int nobj;
        int box_size;
        int ncutout;

        string cutouts_filename;
};


CutoutMaker::CutoutMaker(const CoaddCatalog *coaddCat,
                         const ConfigFile *input_params) :

        id(coaddCat->getIdList()), 
        pos(coaddCat->getPosList()),
        skypos(coaddCat->getSkyPosList()),
        flags(coaddCat->getFlagsList()), 
        skybounds(coaddCat->getSkyBounds()),
        params(*input_params),
        fits(NULL),
        ncutout(0)

{
    this->load_coadd_image();

    this->imlist_path=this->params["coadd_srclist"];

    this->nobj = coaddCat->size();

    this->load_image_list();
    this->set_box_size();

    this->cutout_count.resize(this->nobj,0);
    this->cutout_pos.resize(this->nobj);

    // constant for now
    this->box_sizes.resize(this->nobj, this->box_size);

    this->orig_file_id.resize(this->nobj);
    this->orig_pos.resize(this->nobj);
    this->orig_start_col.resize(this->nobj);
    this->orig_start_row.resize(this->nobj);

    this->cutouts_filename=this->params["cutouts_file"];

}

CutoutMaker::~CutoutMaker()
{
    this->close_fits();
}

void CutoutMaker::load_coadd_image()
{
    // I'm testing on DC6b, different structure so
    // must put in hdu by hand and not use weight image
    int hdu=2;
    int weight_hdu=3;
    this->coadd_image.load(this->params["coaddimage_file"],
                           hdu);
    this->coadd_weight_image.load(this->params["coaddimage_file"],
                                  weight_hdu);

    cerr<<"  image dims "
        <<this->coadd_image.getYMax()<<" "
        <<this->coadd_image.getXMax()<<"\n";
    cerr<<"  weight dims "
        <<this->coadd_weight_image.getYMax()<<" "
        <<this->coadd_weight_image.getXMax()<<"\n";

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

static fitsfile *_open_fits(string filename)
{
    fitsfile *fits=NULL;
    int fitserr=0;

    fits_open_file(&fits,filename.c_str(),READWRITE,&fitserr);
    if (fitserr != 0) {
        fits_report_error(stderr,fitserr);
        throw WriteException(
                "Error opening fits file " + filename);
    }

    return fits;
}


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

void CutoutMaker::close_fits()
{
    if (this->fits) {
        this->fits=_close_fits(this->fits);
    }
}

// at least size 2 because CCfits requires
// length 2 for array columns
int CutoutMaker::get_max_cutouts()
{
    int max_cutouts=2;
    for (int i=0; i<this->nobj; i++) {
        int count=this->orig_file_id[i].size();
        if (count > max_cutouts) {
            max_cutouts=count;
        }
    }
    return max_cutouts;
}
template <typename T>
void to_valarray_vec(const std::vector<std::vector<T> > *vvec,
                     std::vector<std::valarray<T> > *vval,
                     int arrsize)
{
    int n=vvec->size();
    vval->resize(n);
    for (int i=0; i<n; i++) {
        int s=vvec->at(i).size();
        if (s > arrsize) {
            std::stringstream ss;
            ss<<"number of cutouts "<<s<<" exceeded maximum "<<arrsize;
            throw std::runtime_error(ss.str());
        }

        vval->at(i).resize(arrsize,DEFVAL);
        for (int j=0; j<s; j++) {
            vval->at(i)[j] = vvec->at(i)[j];
        }
    }
}


static void vposition_to_xy(const std::vector<std::vector<Position> > *pvec,
                            std::vector<std::valarray<double> > *xvec,
                            std::vector<std::valarray<double> > *yvec,
                            int arrsize)
{
    int n=pvec->size();
    xvec->resize(n);
    yvec->resize(n);

    for (int i=0; i<n; i++) {
        int s=pvec->at(i).size();
        if (s > arrsize) {
            std::stringstream ss;
            ss<<"number of cutouts "<<s<<" exceeded maximum "<<arrsize;
            throw std::runtime_error(ss.str());
        }

        xvec->at(i).resize(arrsize,DEFVAL);
        yvec->at(i).resize(arrsize,DEFVAL);
        for (int j=0; j<s; j++) {
            xvec->at(i)[j] = pvec->at(i)[j].getX();
            yvec->at(i)[j] = pvec->at(i)[j].getY();
        }
    }

}
                                 

void CutoutMaker::write_catalog_filenames(CCfits::FITS *fits)
{
    int nfiles=this->image_file_list.size();

    std::vector<string> col_names(1);
    std::vector<string> col_fmts(1);
    std::vector<string> col_units(1);

    std::stringstream ss;
    ss<<this->max_image_filename_len<<"A";
    col_names[0] = "filename";
    col_fmts[0] = ss.str();

    CCfits::Table* table = fits->addTable("images_info",nfiles,
                                          col_names,col_fmts,col_units);

    int firstrow=1;
    table->column("filename").write(this->image_file_list,firstrow);
}
void CutoutMaker::write_catalog()
{


    cerr<<"opening output file: "<<this->cutouts_filename<<"\n";
    CCfits::FITS fits("!"+this->cutouts_filename, CCfits::Write);

    cerr<<"  writing catalog\n";
    int max_cutouts=this->get_max_cutouts();

    std::stringstream ss;
    ss<<max_cutouts<<"J";
    string jfmt = ss.str();
    ss.str("");
    ss<<max_cutouts<<"D";
    string dfmt = ss.str();


    std::vector<string> col_names;
    std::vector<string> col_fmts;

    col_names.push_back("id");
    col_fmts.push_back("1J");

    //col_names.push_back("file_id");
    //col_fmts.push_back("1J");

    //col_names.push_back("ext");
    //col_fmts.push_back("1J");

    col_names.push_back("ncutout");
    col_fmts.push_back("1J");

    col_names.push_back("box_size");
    col_fmts.push_back("1J");
    //col_names.push_back("start_row");
    //col_fmts.push_back("1J");

    col_names.push_back("start_row");
    col_fmts.push_back(jfmt);

    col_names.push_back("file_id");
    col_fmts.push_back(jfmt);
    col_names.push_back("col_orig");
    col_fmts.push_back(dfmt);
    col_names.push_back("row_orig");
    col_fmts.push_back(dfmt);

    col_names.push_back("x_cutout");
    col_fmts.push_back(dfmt);
    col_names.push_back("y_cutout");
    col_fmts.push_back(dfmt);

    // dummy
    std::vector<string> col_units(col_fmts.size());

    CCfits::Table* table = fits.addTable("object_data",this->id.size(),
                                         col_names,col_fmts,col_units);

    int firstrow=1;
    table->column("id").write(this->id,firstrow);
    table->column("ncutout").write(this->cutout_count,firstrow);
    table->column("box_size").write(this->box_sizes,firstrow);

    std::vector<std::valarray<int> > orig_file_id;
    to_valarray_vec(&this->orig_file_id, &orig_file_id, max_cutouts);
    table->column("file_id").writeArrays(orig_file_id,firstrow);

    std::vector<std::valarray<int> > start_row;
    to_valarray_vec(&this->start_row, &start_row, max_cutouts);
    table->column("start_row").writeArrays(start_row,firstrow);

    std::vector<std::valarray<double> > x;
    std::vector<std::valarray<double> > y;

    vposition_to_xy(&this->orig_pos, &x, &y, max_cutouts);
    table->column("col_orig").writeArrays(x,firstrow);
    table->column("row_orig").writeArrays(y,firstrow);

    vposition_to_xy(&this->cutout_pos, &x, &y, max_cutouts);
    table->column("x_cutout").writeArrays(x,firstrow);
    table->column("y_cutout").writeArrays(y,firstrow);

    this->write_catalog_filenames(&fits);

    // RAII will write and flush the data
}


// Image has no resize method, so we create a new image and send it
// back
//
// assume all are the same shape
template <typename T>
Image<T> *make_mosaic(const std::vector<Image<T> > *vec)
{
    int nx=vec->at(0).getXMax();
    int ny=vec->at(0).getYMax();

    int nim=vec->size();
    int nxtot = nx*nim;

    Image<T> *mosaic=new Image<T>(nxtot, ny);

    for (int i=0; i<nim; i++) {
        // sub-image extraction is [min,max)
        const Image<T> *cutout = &vec->at(i);

        int startx= i*nx;
        int endx  = (i+1)*nx;
        Image<T> subim=mosaic->subImage(startx,endx,0,ny); // shared storage
        subim = (*cutout);
    }
    return mosaic;
}


template <typename T> 
inline int getBitPix() { return 0; }
template <> 
inline int getBitPix<double>() { return DOUBLE_IMG; }
template <> 
inline int getBitPix<float>() { return FLOAT_IMG; }

template <typename T> 
inline int getDataType() { return 0; }
template <> 
inline int getDataType<double>() { return TDOUBLE; }
template <> 
inline int getDataType<float>() { return TFLOAT; }


// box_size*ncutout rows by box_size columns
template <typename T>
void create_mosaic(fitsfile* fits, int box_size, int ncutout)
{
    int fitserr=0;
    int n_axes = 2;

    int ncols=box_size;
    int nrows=box_size*ncutout;

    long dims[2] = { ncols, nrows };

    Assert(getBitPix<T>());
    int bitpix = getBitPix<T>();
    fits_create_img(fits, bitpix, n_axes, dims, &fitserr);
    if (fitserr != 0) {
        fits_report_error(stderr,fitserr);
        throw WriteException("Error creating mosaic");
    }


}
// we assume the catalog and file metadata have already been written, just
// reopen and add the new image hdu(s)

void CutoutMaker::open_fits()
{
    this->fits=_open_fits(this->cutouts_filename);
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

    if (fname.size() > this->max_image_filename_len) {
        this->max_image_filename_len=fname.size();
    }
}
void CutoutMaker::load_image_list()
{
    if (!DoesFileExist(this->imlist_path)) {
        cerr<<"file not found: "<<this->imlist_path<<"\n";
        throw FileNotFoundException(this->imlist_path);
    }

    dbg<<"Opening coadd srclist\n";
    std::ifstream flist(this->imlist_path.c_str(), std::ios::in);
    if (!flist) {
        string err="Unable to open source list file " + this->imlist_path;
        cerr<<err<<std::endl;
        throw ReadException(err);
    }

    this->max_image_filename_len=0;

    // first push the coadd, always first
    this->push_image_file(this->params["coaddimage_file"]);
    string fname;
    while (flist >> fname) {
        dbg<<"  "<<fname<<std::endl;
        this->push_image_file(fname);
    }

#ifdef MAXIMAGES 
    this->image_file_list.resize(MAXIMAGES);
#endif

    cerr<<"read "
        <<this->image_file_list.size()
        <<" se images"<<std::endl;
}


int CutoutMaker::count_noflags()
{
    int ngood = 0;
    for(int i=0;i<this->nobj;i++) {
        if (this->flags[i]== 0) {
            ngood++;
        }
    }
    return ngood;
}



static double get_pixscale(const Transformation *trans)
{
    DSmallMatrix22 D;
    Position pos(1000,1000);
    trans->getDistortion(pos,D);
    double det = std::abs(D.TMV_det());
    double pixel_scale = sqrt(det); // arcsec/pixel
    return pixel_scale;
}


/* 
   use the rough transform as input to the nonlinear finder
   If it doesn't work, use the rough one
*/
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

// we may make this be different for different galaxies
// constant for now
// always force odd
void CutoutMaker::set_box_size()
{
    this->box_size = this->params["multishear_cutout_size"];
    if ( (this->box_size % 2) == 0) {
        this->box_size += 1;
    }
}

// position is position in the original image
void CutoutMaker::append_cutout_info(const Position *pos, int index)
{
    double px=int(pos->getX());
    double py=int(pos->getY());

    // note the box size is always odd
    int offset = (this->box_size-1)/2;
    int startx = px-offset;
    int starty = py-offset;

    Position cpos=(*pos);
    cpos -= Position(startx,starty);

    this->orig_pos[index].push_back((*pos));
    this->cutout_pos[index].push_back(cpos);

    this->orig_start_col[index].push_back(startx);
    this->orig_start_row[index].push_back(starty);

    this->ncutout++;

}



string CutoutMaker::setup_se_file(int ifile)
{
    string image_file = this->image_file_list[ifile];
    this->print_file_progress(ifile);

    // for image reading and transform reading
    this->params["image_file"] = image_file;
    SetRoot(this->params,image_file);

    return image_file;
}



static void get_image_shape(string filename, int hdu, long *ncol, long *nrow)
{

    int fitserr=0;

    fitsfile *fits=_open_fits(filename);

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
    fits_get_img_size(fits, ndim,  dims, &fitserr);
    if (fitserr != 0) {
        fits_report_error(stderr,fitserr);
        throw ReadException("Error getting dims from file "+filename);
    }

    fits=_close_fits(fits);

    *ncol=dims[0];
    *nrow=dims[1];
}



void CutoutMaker::set_start_rows()
{
    long current_row=0;
    this->start_row.resize(this->nobj);

    for (int i=0; i<this->nobj; i++) {
        int ncut=this->cutout_pos[i].size();
        if (ncut > 0) {
            this->start_row[i].resize(ncut);
            for (int j=0; j<ncut; j++) {
                this->start_row[i][j] = current_row;
                current_row += this->box_size;
            }
        }
    }
}

void CutoutMaker::write_cutout(const Image<double> *image, int iobj, int icut)
{
    int xmax=image->getXMax();
    int ymax=image->getYMax();

    int startx=this->orig_start_col[iobj][icut];
    int starty=this->orig_start_row[iobj][icut];

    Image<double> tim(this->box_size,this->box_size);
    for (int ix=0, x=startx; ix<this->box_size; ix++, x++) {
        if (x < 0 || x >= xmax ) continue;

        for (int iy=0, y=starty; iy<this->box_size; iy++, y++) {
            if (y < 0 || y >= ymax ) continue;

            tim(ix,iy) = (*image)(x,y);
        }
    }

    long col=1;
    long row=1+this->start_row[iobj][icut];
    tim.write_sub(this->fits, col, row);

}

void CutoutMaker::print_file_progress(int ifile)
{
    cerr<<"    "
        <<(ifile+1)
        <<"/"
        <<this->image_file_list.size()
        <<" "
        <<this->image_file_list[ifile]<<"\n";
}

void CutoutMaker::write_mosaics_from_coadd(enum cutout_type cut_type)
{
    this->print_file_progress(0);
    for (int i=0; i<this->nobj; ++i) {
        int ncut=this->cutout_pos[i].size();
        for (int j=0; j<ncut; j++) {
            if (this->orig_file_id[i][j] == 0) {
                if (cut_type==CUTOUT_WEIGHT) {
                    //this->write_cutout(this->coadd_weight_image.get(), i, j);
                    this->write_cutout(&this->coadd_weight_image, i, j);
                } else {
                    this->write_cutout(&this->coadd_image, i, j);
                }
            }
        }
    }
}
void CutoutMaker::write_mosaics_from_file(int ifile, enum cutout_type cut_type)
{
    if (ifile==0) {
        return;
    }

    string filename=this->setup_se_file(ifile);

    auto_ptr<Image<double> > weight_image;
    Image<double> image(this->params,weight_image);

    for (int i=0; i<this->nobj; ++i) {
        int ncut=this->cutout_pos[i].size();
        for (int j=0; j<ncut; j++) {
            if (this->orig_file_id[i][j] == ifile) {
                if (cut_type==CUTOUT_WEIGHT) {
                    this->write_cutout(weight_image.get(), i, j);
                } else {
                    this->write_cutout(&image, i, j);
                }
            }
        }
    }

}

// write the idividual cutouts into the big mosaic
void CutoutMaker::write_mosaic(enum cutout_type cut_type)
{
    this->open_fits();

    create_mosaic<double>(this->fits, this->box_size, this->ncutout);

    this->write_mosaics_from_coadd(cut_type);

    const int nfiles = this->image_file_list.size();
    for (int ifile=1; ifile<nfiles; ++ifile) {
        this->write_mosaics_from_file(ifile, cut_type);
    }

    this->close_fits();
}



// the input bounds are on the sky
// can avoid reading the image by just grabbing the header keywords
void CutoutMaker::get_file_cutout_info(int ifile, const Bounds *bounds)
{

    string filename=this->setup_se_file(ifile);
    long ncol=0, nrow=0;
    int hdu=GetHdu(this->params,"image",filename,1);
    get_image_shape(filename, hdu, &ncol, &nrow);

    // Read transformation (x,y) -> (ra,dec)
    Transformation trans(this->params);
    Transformation inv_trans;

    // rough inverse trans, only valid on the bounds of this image
    Bounds se_bounds(0, ncol, 0, nrow);
    Bounds inv_bounds = inv_trans.makeInverseOf(trans,se_bounds,4);

    if (!this->skybounds.intersects(inv_bounds)) {
        cerr<<"Skipping "<<filename<<" because inv_bounds doesn't intersect\n";
        return;
    }

    long nkeep=0, ngood=0, ntry=0;
    Position pos(0,0);

    for (int i=0; i<this->nobj; ++i) {
        if (this->flags[i]) continue;
        if (!bounds->includes(this->skypos[i])) continue;
        if (!inv_bounds.includes(this->skypos[i])) continue;

        ntry++;
        ngood+=get_image_pos(&trans, &inv_trans, &this->skypos[i], &pos);

        if (!se_bounds.includes(pos)) continue;

        // always push the coadd image on first
        if (this->orig_file_id[i].size()==0) {
            this->append_cutout_info(&this->pos[i],i);
            this->orig_file_id[i].push_back(0); // coadd is always at zero
        }
        this->append_cutout_info(&pos, i);
        this->orig_file_id[i].push_back(ifile);
        nkeep++;
    }

    cerr<<"        kept "<<nkeep<<"/"<<this->nobj
        <<" good trans "<<ngood<<"/"<<ntry<<"\n";
}


// loop through and find the cutout boundaries, etc.
// in each file
void CutoutMaker::get_cutout_info()
{

    const int nfiles = this->image_file_list.size();
    for (int ifile=1; ifile<nfiles; ++ifile) {
        this->get_file_cutout_info(ifile,&this->skybounds);
    }
    cerr<<"Found "<<this->ncutout<<" cutouts\n";

    // determine where the objects will go in the mosaic
    this->set_start_rows();
}

static void make_cutouts(ConfigFile *params) 
{
    CoaddCatalog coaddcat(*params);
    coaddcat.read();

    CutoutMaker maker(&coaddcat, params);
    cerr<<"number without flags: "<<maker.count_noflags();

    maker.print_area_stats();
    maker.get_cutout_info();

    maker.write_catalog();
    cerr<<"  writing cutout mosaic\n";
    maker.write_mosaic(CUTOUT_IMAGE);
    cerr<<"  writing weight mosaic\n";
    maker.write_mosaic(CUTOUT_WEIGHT);
}


static void usage_and_exit()
{
    cerr<<"make-cutouts config_file coaddimage_file= coaddcat_file= coadd_srclist= cutouts_file=\n";
    exit(1);
}
bool check_params(const ConfigFile *params)
{

    if (!params->keyExists("coaddimage_file")
            || !params->keyExists("coaddcat_file")
            || !params->keyExists("coadd_srclist")
            || !params->keyExists("cutouts_file") ) {
        return false;
    }

    return true;
}

int main(int argc, char **argv)
{
    if (argc < 2)
        usage_and_exit();

    ConfigFile params;
    if (BasicSetup(argc,argv,params,"make_cutouts")) return EXIT_FAILURE;

    if (!check_params(&params))
        usage_and_exit();

    make_cutouts(&params);
}
