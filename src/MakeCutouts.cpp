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
       - add the coadd image in the first entry of the mosaic
       - add wcs info to the second extension along with the filenames.
       - add weight images
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

        void load_coadd_image();
        string get_cutout_filename(int filenum);
        void start_new_fits();

        template <typename T> void write_mosaic(const Image<T> *mosaic);

        void write_mosaics();

        void get_cutouts_bound(const Bounds& bounds);
        void clear_cutout_list();

        void push_image_file(string fname);
        void load_image_list();

        std::vector<Bounds> split_bounds();

        void append_cutout_and_pos(const Image<double> *image,
                                   const Position *pos,
                                   int index);

        string setup_se_file(int ifile);
        bool load_file_cutouts(int ifile,
                               const Bounds *bounds);

        bool make_section_cutouts(const Bounds *bounds);
        void make_cutouts();
        void write_catalog_filenames(CCfits::FITS *fits);
        void write_catalog();

        int get_nobj_with_pixels() const;

        void print_area_stats();

        void bombout_memory();
        void set_memory_usage();

        bool vmem_ok();

        int get_box_size(int index);

        int count_noflags();
        int get_max_cutouts();

        void close_files();

    private:

        Image<double> coadd_image;
        auto_ptr<Image<double> > coadd_weight_image;

        // for input

        std::vector<long> id;
        std::vector<Position> pos;
        std::vector<Position> skypos;
        std::vector<long> flags;

        // these tell us how to find the cutout
        // for an object
        std::vector<int> cutout_file_id;
        std::vector<int> cutout_ext; // 0 offset
        std::vector<int> cutout_count;

        std::vector<std::vector<int> > file_id_orig; // id of SE file
        std::vector<std::vector<Position> > orig_pos; // position in SE image
        std::vector<std::vector<Position> > cutout_pos; // positions in the cutout
        std::vector<std::vector<Image<double> > > cutout_list;

        Bounds skybounds; 
        string imlist_path;


        std::vector<string> image_file_list;
        int max_image_filename_len;

        ConfigFile params;

        fitsfile *fits;
        
        string current_filename;
        int current_file_id; // the file number
        int current_hdu; // where we are in the current file

        double max_mem;
        int nobj;
};

static void bombout_resplit(const ConfigFile *params)
{
    cerr<<"No more allowed resplits.\n";
    cerr<<"Either reduce multishear_section_size, or\n";
    cerr<<"increase max_vmem or multishear_max_resplits\n";
    cerr<<"Current values = "<<
        (*params)["multishear_section_size"]<<" , "<<
        (*params)["max_vmem"]<<" , "<<
        (*params)["multishear_max_resplits"]<<std::endl;
    throw ProcessingException(
            "Memory exceeded maximum allowed, and unable to "
            "recover by resplitting bounds");
}

// Split this section and push it on the end of the Bounds vector
static void resplit(std::vector<Bounds> *section_bounds, const Bounds* bounds)
{
    std::vector<Bounds> split = bounds->quarter();
    section_bounds->insert( section_bounds->end(),
                            split.begin(),split.end());

    dbg<<"Split bounds into four new bounds:\n";
    dbg<<"b1 = "<<split[0];
    dbg<<"b2 = "<<split[1];
    dbg<<"b3 = "<<split[2];
    dbg<<"b4 = "<<split[3];

}





CutoutMaker::CutoutMaker(const CoaddCatalog *coaddCat,
                         const ConfigFile *input_params) :

        id(coaddCat->getIdList()), 
        pos(coaddCat->getPosList()),
        skypos(coaddCat->getSkyPosList()),
        flags(coaddCat->getFlagsList()), 
        skybounds(coaddCat->getSkyBounds()),
        params(*input_params),
        current_file_id(-1),
        current_hdu(99999),
        fits(NULL)

{
    this->load_coadd_image();

    this->imlist_path=this->params["coadd_srclist"];

    this->nobj = coaddCat->size();
    this->cutout_list.resize(this->nobj);
    this->cutout_file_id.resize(this->nobj,-9999);
    this->cutout_ext.resize(this->nobj,-9999);
    this->cutout_count.resize(this->nobj,0);

    this->file_id_orig.resize(this->nobj);
    this->orig_pos.resize(this->nobj);
    this->cutout_pos.resize(this->nobj);

    this->max_mem = this->params.read("max_vmem",64)*1024.;

    this->load_image_list();
}

CutoutMaker::~CutoutMaker()
{
    this->close_files();
}

void CutoutMaker::load_coadd_image()
{
    // I'm testing on DC6b, different structure so
    // must put in hdu by hand and not use weight image
    int hdu=2;
    this->coadd_image.load(this->params["coaddimage_file"],
                           hdu);
    cerr<<"  image dims ";
    cerr<<this->coadd_image.getYMax()<<" "
        <<this->coadd_image.getXMax()<<"\n";
}


static fitsfile *close_fits(fitsfile *fits)
{
    int fitsErr=0;
    fits_close_file(fits, &fitsErr);
    if (fitsErr != 0) {
        fits_report_error(stderr,fitsErr);
        throw WriteException("Error closing fits file");
    }
    return NULL;
}
static fitsfile *open_new_fits(string filename)
{
    fitsfile *fits=NULL;
    int fitsErr=0;

    fits_create_file(&fits,("!"+filename).c_str(),&fitsErr);
    if (fitsErr != 0) {
        fits_report_error(stderr,fitsErr);
        throw WriteException(
                "Error creating fits file " + filename);
    }

    return fits;
}



void CutoutMaker::close_files()
{
    if (this->fits) {
        this->fits=close_fits(this->fits);
    }
}

// at least size 2 because CCfits requires
// length 2 for array columns
int CutoutMaker::get_max_cutouts()
{
    int max_cutouts=2;
    for (int i=0; i<this->nobj; i++) {
        int count=this->file_id_orig[i].size();
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

        vval->at(i).resize(arrsize,-9999);
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

        xvec->at(i).resize(arrsize,-9999);
        yvec->at(i).resize(arrsize,-9999);
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

    int start_row=1;
    table->column("filename").write(this->image_file_list,start_row);
}
void CutoutMaker::write_catalog()
{

    string catname="cutouts-cat.fits";
    cerr<<"writing catalog: "<<catname<<"\n";
    CCfits::FITS fits("!"+catname, CCfits::Write);

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

    col_names.push_back("file_id");
    col_fmts.push_back("1J");

    col_names.push_back("ext");
    col_fmts.push_back("1J");

    col_names.push_back("ncutout");
    col_fmts.push_back("1J");

    col_names.push_back("file_id_orig");
    col_fmts.push_back(jfmt);
    col_names.push_back("x_orig");
    col_fmts.push_back(dfmt);
    col_names.push_back("y_orig");
    col_fmts.push_back(dfmt);

    col_names.push_back("x_cutout");
    col_fmts.push_back(dfmt);
    col_names.push_back("y_cutout");
    col_fmts.push_back(dfmt);

    // dummy
    std::vector<string> col_units(col_fmts.size());

    CCfits::Table* table = fits.addTable("cutouts",this->id.size(),
                                         col_names,col_fmts,col_units);

    int start_row=1;
    table->column("id").write(this->id,start_row);
    table->column("file_id").write(this->cutout_file_id,start_row);
    table->column("ext").write(this->cutout_ext,start_row);
    table->column("ncutout").write(this->cutout_count,start_row);

    std::vector<std::valarray<int> > file_id_orig;
    to_valarray_vec(&this->file_id_orig, &file_id_orig, max_cutouts);
    table->column("file_id_orig").writeArrays(file_id_orig,start_row);

    std::vector<std::valarray<double> > x;
    std::vector<std::valarray<double> > y;

    vposition_to_xy(&this->orig_pos, &x, &y, max_cutouts);
    table->column("x_orig").writeArrays(x,start_row);
    table->column("y_orig").writeArrays(y,start_row);

    vposition_to_xy(&this->cutout_pos, &x, &y, max_cutouts);
    table->column("x_cutout").writeArrays(x,start_row);
    table->column("y_cutout").writeArrays(y,start_row);

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


// need to base this on the input coadd image
string CutoutMaker::get_cutout_filename(int filenum)
{
    std::stringstream ss;
    char numstr[10];

    std::sprintf(numstr,"%04d", filenum);

    ss << "cutouts-" << numstr << ".fits";
    return ss.str();
}

// upon entry, current_file_id should be for the last write, or -1
void CutoutMaker::start_new_fits()
{
    std::stringstream ss;
    this->current_file_id++;

    this->current_filename = this->get_cutout_filename(this->current_file_id);

    cerr<<"starting new cutout file: "<<this->current_filename<<"\n";

    if (this->fits) {
        this->fits=close_fits(this->fits);
    }
    this->fits = open_new_fits(this->current_filename);
}

// upon entry, curren_hdu and current_file_id should be that of the last write
//
// upon exit, current_hdu and current_file_id are updated to the file and hdu
// just written (current_file_id is modified by start_new_fits)

template <typename T>
void CutoutMaker::write_mosaic(const Image<T> *mosaic)
{

    // up on entry, current hdu and fileid are from the last write
    this->current_hdu++;
    if (this->current_hdu > MAX_HDU) {
        this->start_new_fits();
        this->current_hdu=1;
    }

    mosaic->write(this->fits);
}

// append a new extension for each object that has images
// if working on a particular section, many will not 
// have cutouts
void CutoutMaker::write_mosaics()
{
    for (int i=0; i<this->nobj; i++) {
        int ncutout=this->cutout_list[i].size();
        if (ncutout > 0) {
            const std::vector<Image<double> > *cutouts=
                &this->cutout_list[i];
            Image<double> *mosaic=make_mosaic(cutouts);

            mosaic->set_compression(RICE_1);
            this->write_mosaic(mosaic);
            delete mosaic;

            this->cutout_file_id[i] = this->current_file_id;
            this->cutout_ext[i] = this->current_hdu-1;
            this->cutout_count[i] = ncutout;
        }
    }
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
    this->push_image_file(this->params["coaddcat_file"]);
    string fname;
    while (flist >> fname) {
        dbg<<"  "<<fname<<std::endl;
        this->push_image_file(fname);
    }

    cerr<<"read "
        <<this->image_file_list.size()
        <<" se images"<<std::endl;
}
void CutoutMaker::clear_cutout_list()
{
    dbg<<"Start getPixels: memory_usage = "<<memory_usage()<<std::endl;
    const int n = this->cutout_list.size();
    for (int i=0; i<n; ++i) {
        this->cutout_list[i].clear();
    }
    dbg<<"After clear: memory_usage = "<<memory_usage()<<std::endl;
}

int CutoutMaker::get_nobj_with_pixels() const
{
    int nobj_withpix=0;
    for (int i=0; i<this->nobj; ++i) {
        if (this->cutout_list[i].size() > 0) {
            nobj_withpix++;
        }
    }
    return nobj_withpix;
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



std::vector<Bounds> CutoutMaker::split_bounds()
{
    double sectionSize = this->params.get("multishear_section_size");
    xdbg<<"sectionSize = "<<sectionSize<<std::endl;
    sectionSize *= 60.; // arcmin -> arcsec

    double xRange = 
        this->skybounds.getXMax() - this->skybounds.getXMin();

    double yRange = 
        this->skybounds.getYMax() - this->skybounds.getYMin();

    double dec = this->skybounds.getCenter().getY();
    double cosdec = cos(dec / ARCSEC_PER_RAD);
    xRange *= cosdec;

    int nx = int(floor(xRange/sectionSize))+1;
    int ny = int(floor(yRange/sectionSize))+1;
    return this->skybounds.divide(nx,ny);
}


void CutoutMaker::bombout_memory()
{
    bool des_qa = this->params.read("des_qa",false); 

    dbg<<"Caught bad_alloc\n";
    double mem = memory_usage(dbgout);
    double peak_mem = peak_memory_usage();
    dbg<<"memory usage = "<<mem<<" MB\n";
    dbg<<"peak memory usage = "<<peak_mem<<" MB\n";
    if (des_qa) cerr<<"STATUS5BEG ";
    cerr
        << "Memory exhausted in MultShearCatalog.\n"
        << "Actual Virtual Memory Usage = "
        << mem<<" MB \n"
        << "Try reducing multishear_section_size or "
        << "reducing mam_vmem\n"
        << "(Current values = "
        << this->params["multishear_section_size"]
        << " , "<< this->params["max_vmem"]<<")";
    if (des_qa) cerr<<" STATUS5END";
    cerr<<std::endl;
    exit(1);
}

void CutoutMaker::set_memory_usage()
{
    double mem = memory_usage(dbgout);
    double peak_mem = peak_memory_usage();
    double max_mem = double(this->params["max_vmem"])*1024.;
    this->params["peak_mem"] = peak_mem; // Keep track...
    dbg<<"Actual memory usage = "<<mem<<" MB\n";
    dbg<<"Peak memory usage = "<<peak_mem<<" MB\n";
    dbg<<"Max allowed memory usage = "<<max_mem<<" MB\n";
}

bool CutoutMaker::vmem_ok()
{
    bool status=true;
    double mem = memory_usage();
    if (mem > this->max_mem) {
        dbg<<"VmSize = "<<mem<<" > max_vmem = "<<max_mem<<std::endl;
        status = false;
    }
    return status;
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
int CutoutMaker::get_box_size(int index)
{
    int box_size = this->params["multishear_cutout_size"];
    if ( (box_size % 2) == 0) {
        box_size += 1;
    }

    return box_size;
}

// append a new image for this object, data copied from the SE image.  Zeros
// are used outside the bounds of the SE image
//
// also append the position within the SE and cutout images
//
// Return the position in the cutout

void CutoutMaker::append_cutout_and_pos(const Image<double> *image,
                                        const Position *pos,
                                        int index)
{
    int box_size = this->get_box_size(index);

    double px=int(pos->getX());
    double py=int(pos->getY());

    int xmax=image->getXMax();
    int ymax=image->getYMax();

    // note the box size is always odd
    int offset = (box_size-1)/2;
    int startx = px-offset;
    int starty = py-offset;

    Image<double> tim(box_size,box_size);
    for (int ix=0, x=startx; ix<box_size; ix++, x++) {
        if (x < 0 || x >= xmax ) continue;

        for (int iy=0, y=starty; iy<box_size; iy++, y++) {
            if (y < 0 || y >= ymax ) continue;

            tim(ix,iy) = (*image)(x,y);
        }
    }

    Position cpos=(*pos);
    cpos -= Position(startx,starty);

    this->cutout_list[index].push_back(tim);
    this->orig_pos[index].push_back((*pos));
    this->cutout_pos[index].push_back(cpos);
}

string CutoutMaker::setup_se_file(int ifile)
{
    string image_file = this->image_file_list[ifile];
    cerr<<"    file: "<<image_file<<" "
        <<(ifile+1)<<"/"<<this->image_file_list.size()<<"\n";

    // for image reading and transform reading
    this->params["image_file"] = image_file;
    SetRoot(this->params,image_file);

    return image_file;
}

// the input bounds are on the sky
bool CutoutMaker::load_file_cutouts(int ifile, const Bounds *bounds)
{

    string filename=this->setup_se_file(ifile);

    auto_ptr<Image<double> > weight_image;
    Image<double> image(this->params,weight_image);

    // Read transformation (x,y) -> (ra,dec)
    Transformation trans(this->params);
    Transformation inv_trans;

    // rough inverse trans, only valid on the bounds of this image
    Bounds se_bounds = image.getBounds();
    Bounds inv_bounds = inv_trans.makeInverseOf(trans,se_bounds,4);

    if (!this->skybounds.intersects(inv_bounds)) {
        cerr<<"Skipping "<<filename<<" because inv_bounds doesn't intersect\n";
        return true;
    }

    long nkeep=0;
    Position pos(0,0);

    for (int i=0; i<this->nobj; ++i) {
        if (this->flags[i]) continue;
        if (!bounds->includes(this->skypos[i])) continue;
        if (!inv_bounds.includes(this->skypos[i])) continue;

        get_image_pos(&trans, &inv_trans, &this->skypos[i], &pos);

        if (!se_bounds.includes(pos)) continue;

        // always push the coadd image on first
        if (this->cutout_list[i].size()==0) {
            this->append_cutout_and_pos(&this->coadd_image,&this->pos[i],i);
            this->file_id_orig[i].push_back(0); // coadd is always at zero
        }
        this->append_cutout_and_pos(&image, &pos, i);
        this->file_id_orig[i].push_back(ifile);
        nkeep++;

        if (!this->vmem_ok()) return false;
    }

    cerr<<"        kept "<<nkeep<<"/"<<this->nobj
        <<"  memory_usage: "<<memory_usage()/1024.<<" Gb\n";;
    return true;
}

// note ifile loop starts at 1 because the coadd is in
// position zero
bool CutoutMaker::make_section_cutouts(const Bounds *bounds)
{

    this->clear_cutout_list();

    try {
        dbg<<"getting cutouts for b = "<<bounds<<std::endl;
        memory_usage(dbgout);

        //const int nfiles = this->image_file_list.size();
        const int nfiles = 10;

        for (int ifile=1; ifile<nfiles; ++ifile) {
            if (!this->load_file_cutouts(ifile,bounds)) {
                this->clear_cutout_list();
                return false;
            }
        }
    } catch (std::bad_alloc) {
        this->bombout_memory();
    }

    dbg <<"Done getting cutouts\n";
    this->set_memory_usage();
    return true;
}



//   Break the area up into sections and load each section separately.
//   It is most efficient to load everything at once, but it can take 
//   a lot of memory.  So this bit keeps the memory requirement for each
//   section manageable.  So long as each section does a significant amount
//   of work, the extra I/O time won't be much of an issue. 


void CutoutMaker::make_cutouts()
{

    std::vector<Bounds> section_bounds = this->split_bounds();
    int nsec = section_bounds.size();
    int nresplit = this->params.read("multishear_max_resplits",1);

    for(int i=0;i<nsec;++i) {

        cerr<<"Starting section ";
        cerr<<(i+1)<<"/"<<nsec<<std::endl;

        // Load the pixel information for each galaxy in the section.
        if (!this->make_section_cutouts(&section_bounds[i])) {

            cerr<<"Section exceeded maximum memory usage.\n";
            cerr<<"Will try splitting it up and continuing.\n";
            if (nresplit > 0) {
                resplit(&section_bounds, &section_bounds[i]);

                nresplit--;
                cerr<<nresplit<<" more resplits allowed.\n";

                continue;
            } else {
                bombout_resplit(&this->params);
            }
        }
        cerr<<this->get_nobj_with_pixels()<<" obj in section.\n";

        cerr<<"writing mosaics\n";
        this->write_mosaics();
    }


}

static void make_cutouts(ConfigFile *params) 
{
    CoaddCatalog coaddcat(*params);
    coaddcat.read();

    CutoutMaker maker(&coaddcat, params);
    cerr<<"number without flags: "<<maker.count_noflags();

    maker.print_area_stats();
    maker.make_cutouts();
    maker.write_catalog();
}

int main(int argc, char **argv)
{
    ConfigFile params;
    if (BasicSetup(argc,argv,params,"make_cutouts")) return EXIT_FAILURE;

    make_cutouts(&params);
}
