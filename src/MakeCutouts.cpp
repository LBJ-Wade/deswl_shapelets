/*
   Author: Erin Sheldon, BNL.  This began as a hack on the multishear code.
*/
#include <valarray>
#include <sys/time.h>
#include <iostream>
#include <fstream>

#include "Image.h"
#include "FittedPsf.h"
#include "CoaddCatalog.h"
#include "MultiShearCatalog.h"
#include "BasicSetup.h"

const double ARCSEC_PER_RAD = 206264.806247;

// If desired stop after this many sections.
//#define ENDAT 2

class CutoutMaker
{
    public:
        CutoutMaker(const CoaddCatalog *coaddCat, 
                    std::string imlist,
                    const ConfigFile *input_params);

        ~CutoutMaker() {}


        void get_cutouts_bound(const Bounds& bounds);
        void clear_cutout_list(void);
        void load_image_list(void);

        std::vector<Bounds> split_bounds();

        void append_image_cutout(const Image<double> *image,
                                 int index,
                                  const Position *pos);
        bool load_file_cutouts(std::string filename, 
                               const Bounds *bounds);
        bool make_section_cutouts(const Bounds *bounds);
        void make_cutouts(void);

        int get_ngals_with_pixels() const;

        void print_stats(void);

        double calc_memory_footprint(void);

        void bombout_memory(void);
        void set_memory_usage(void);

        int get_box_size(int index);

    private:
        std::vector<long> id;
        std::vector<Position> chippos;
        std::vector<Position> skypos;
        std::vector<long> flags;
        Bounds skybounds; 
        std::string imlist_path;

        std::vector<std::vector<Image<double> > > cutout_list;

        std::vector<std::string> image_file_list;

        ConfigFile params;
};

template <typename T>
inline long long getMemoryFootprint(const T& x)
{
    return sizeof(x); 
}

template <typename T, class Alloc>
inline long long getMemoryFootprint(const std::vector<T,Alloc>& x)
{
    long long res=sizeof(x);
    const int xsize = x.size();
    for(int i=0;i<xsize;++i) res += getMemoryFootprint(x[i]);
    return res;
}

// assume the sub vector contains fixed-size images
template <typename T>
inline long long getMemoryFootprint(
        const std::vector<std::vector<Image<T> > > &x)
{
    long long res=sizeof(x);
    long long xsize=x.size();
    long long nimage=0;
    long long nobj_used=0;
    long long elsize=sizeof(T);

    for(int i=0; i<xsize; ++i) {
        long long nsub = x[i].size();
        nimage+= nsub;
        if (nsub > 0) {
            nobj_used+=1;
            long long imsize=x[i][0].getXMax()*x[i][0].getYMax();
            nsub *= imsize*elsize;

            res += nsub;
        }
    }
    dbg<<"  cutout memory usage "<<res/1024./1024.<<" Mb in "<<nimage<<" cutouts\n";
    if (nobj_used > 0) {
        dbg<<"  average cutouts per object: "
            <<nimage<<"/"<<nobj_used
            <<"="<<nimage/float(nobj_used)<<"\n";
    }
    return res;
}


static void bombout_resplit(const ConfigFile *params)
{
    std::cerr<<"No more allowed resplits.\n";
    std::cerr<<"Either reduce multishear_section_size, or\n";
    std::cerr<<"increase max_vmem or multishear_max_resplits\n";
    std::cerr<<"Current values = "<<
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
                         std::string imlist,
                         const ConfigFile *input_params) :

        id(coaddCat->getIdList()), 
        chippos(coaddCat->getPosList()),
        skypos(coaddCat->getSkyPosList()),
        flags(coaddCat->getFlagsList()), 
        skybounds(coaddCat->getSkyBounds()),
        imlist_path(imlist),
        params(*input_params)
{

    const int n = coaddCat->size();
    this->cutout_list.resize(n);

    this->load_image_list();

}


double CutoutMaker::calc_memory_footprint(void)
{
    double mem=0;
    
    mem += getMemoryFootprint(this->skypos);
    mem += getMemoryFootprint(this->chippos);
    mem += getMemoryFootprint(this->flags);
    mem += getMemoryFootprint(this->skybounds);
    mem += getMemoryFootprint(this->cutout_list);

    mem /= 1024.*1024.;
    return mem;
}

void CutoutMaker::print_stats(void)
{
    std::cerr<<"Total bounds are "<<this->skybounds<<std::endl;

    double area = this->skybounds.getArea();
    area /= 3600.;
    std::cerr<<"Raw area = "<<area<<" square arcmin\n";

    double dec = this->skybounds.getCenter().getY();
    double cosdec = cos(dec / ARCSEC_PER_RAD);
    area *= cosdec;
    std::cerr<<"Corrected area = "<<area<<" square arcmin";

    area /= 3600.;
    std::cerr<<" = "<<area<<" square degrees\n";
}


void CutoutMaker::load_image_list(void)
{
    if (!DoesFileExist(this->imlist_path)) {
        std::cerr<<"file not found: "<<this->imlist_path<<"\n";
        throw FileNotFoundException(this->imlist_path);
    }

    dbg<<"Opening coadd srclist\n";
    std::ifstream flist(this->imlist_path.c_str(), std::ios::in);
    if (!flist) {
        std::string err="Unable to open source list file " + this->imlist_path;
        std::cerr<<err<<std::endl;
        throw ReadException(err);
    }

    std::string fname;
    while (flist >> fname) {
        dbg<<"  "<<fname<<std::endl;
        this->image_file_list.push_back(fname);
    }

    std::cerr<<"read "
        <<this->image_file_list.size()
        <<" se images"<<std::endl;
}
void CutoutMaker::clear_cutout_list(void)
{
    dbg<<"Start getPixels: memory_usage = "<<memory_usage()<<std::endl;
    const int n = this->cutout_list.size();
    for (int i=0; i<n; ++i) {
        this->cutout_list[i].clear();
    }
    dbg<<"After clear: memory_usage = "<<memory_usage()<<std::endl;
}

int CutoutMaker::get_ngals_with_pixels() const
{
    const int ngals = this->cutout_list.size();
    int ngals_withpix=0;
    for (int i=0; i<ngals; ++i) {
        if (this->cutout_list[i].size() > 0) {
            ngals_withpix++;
        }
    }
    return ngals_withpix;
}



/*
static int get_ngood(const MultiShearCatalog *shearcat)
{
    int ngood = 0;
    int nflags = shearcat->getFlagsList().size();

    for(int i=0;i<nflags;++i) {
        if (shearcat->getFlagsList()[i] == 0) {
            ++ngood;
        }
    }

    return ngood;
}
*/



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


void CutoutMaker::bombout_memory(void)
{
    bool des_qa = this->params.read("des_qa",false); 

    dbg<<"Caught bad_alloc\n";
    double mem = memory_usage(dbgout);
    double peak_mem = peak_memory_usage();
    dbg<<"memory usage = "<<mem<<" MB\n";
    dbg<<"peak memory usage = "<<peak_mem<<" MB\n";
    if (des_qa) std::cerr<<"STATUS5BEG ";
    std::cerr
        << "Memory exhausted in MultShearCatalog.\n"
        << "Memory Usage in MultiShearCatalog = "
        << this->calc_memory_footprint()<<" MB \n"
        << "Actual Virtual Memory Usage = "
        << mem<<" MB \n"
        << "Try reducing multishear_section_size or "
        << "reducing mam_vmem\n"
        << "(Current values = "
        << this->params["multishear_section_size"]
        << " , "<< this->params["max_vmem"]<<")";
    if (des_qa) std::cerr<<" STATUS5END";
    std::cerr<<std::endl;
    exit(1);
}

void CutoutMaker::set_memory_usage(void)
{
    dbg << "Memory Usage in getting cutouts = "
        << this->calc_memory_footprint()<<" MB \n";
    double mem = memory_usage(dbgout);
    double peak_mem = peak_memory_usage();
    double max_mem = double(this->params["max_vmem"])*1024.;
    this->params["peak_mem"] = peak_mem; // Keep track...
    dbg<<"Actual memory usage = "<<mem<<" MB\n";
    dbg<<"Peak memory usage = "<<peak_mem<<" MB\n";
    dbg<<"Max allowed memory usage = "<<max_mem<<" MB\n";

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


/* use the rough transform as input to the nonlinear finder
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

void CutoutMaker::append_image_cutout(const Image<double> *image,
                                      int index,
                                      const Position *pos)
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

    this->cutout_list[index].push_back(tim);
}

// the input bounds are on the sky
bool CutoutMaker::load_file_cutouts(std::string filename, const Bounds *bounds)
{

    std::auto_ptr<Image<double> > weight_image;
    Image<double> image(this->params,weight_image);

    // Read transformation (x,y) -> (ra,dec)
    Transformation trans(this->params);
    Transformation inv_trans;

    // rough inverse trans, only valid on the bounds of this image
    Bounds se_bounds = image.getBounds();
    Bounds inv_bounds = inv_trans.makeInverseOf(trans,se_bounds,4);

    if (!this->skybounds.intersects(inv_bounds)) {
        std::cerr<<"Skipping "<<filename<<" because inv_bounds doesn't intersect\n";
        return true;
    }

    double max_mem = this->params.read("max_vmem",64)*1024.;
    long ngals=this->cutout_list.size();
    long nkeep=0;
    long ngood=0, ntry=0;

    Position pos(0,0);

    for (int i=0; i<ngals; ++i) {
        if (this->flags[i]) continue;
        if (!bounds->includes(this->skypos[i])) continue;
        if (!inv_bounds.includes(this->skypos[i])) continue;

        ntry++;
        ngood+=get_image_pos(&trans, &inv_trans, &this->skypos[i], &pos);

        if (!se_bounds.includes(pos)) continue;

        this->append_image_cutout(&image, i, &pos);

        double mem = memory_usage();
        if (mem > max_mem) {
            dbg<<"VmSize = "<<mem<<" > max_vmem = "<<max_mem<<std::endl;
            return false;
        }
        nkeep++;
    }

    long long cutmem=getMemoryFootprint(this->cutout_list);
    std::cerr<<"        kept "<<nkeep<<"/"<<ngals
        <<"  memory_usage: "<<memory_usage()/1024.<<" Gb\n";;
    return true;
}
bool CutoutMaker::make_section_cutouts(const Bounds *bounds)
{

    this->clear_cutout_list();

    try {
        dbg<<"getting cutouts for b = "<<bounds<<std::endl;
        memory_usage(dbgout);

        const int nfiles = this->image_file_list.size();

        for (int ifile=0; ifile<nfiles; ++ifile) {
            std::string image_file = this->image_file_list[ifile];
            std::cerr<<"    file: "<<image_file<<" "
                <<(ifile+1)<<"/"<<nfiles<<"\n";
            
            // for image reading and transform reading
            this->params["image_file"] = image_file;
            SetRoot(this->params,image_file);

            if (!this->load_file_cutouts(image_file,bounds)) {
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


void CutoutMaker::make_cutouts(void)
{

    std::vector<Bounds> section_bounds = this->split_bounds();
    int nsec = section_bounds.size();
    int nresplit = this->params.read("multishear_max_resplits",1);

    for(int i=0;i<nsec;++i) {

#ifdef ENDAT
        if (i == ENDAT) break;
#endif

        std::cerr<<"Starting section ";
        std::cerr<<(i+1)<<"/"<<nsec<<std::endl;

        // Load the pixel information for each galaxy in the section.
        if (!this->make_section_cutouts(&section_bounds[i])) {

            std::cerr<<"Section exceeded maximum memory usage.\n";
            std::cerr<<"Will try splitting it up and continuing.\n";
            if (nresplit > 0) {
                resplit(&section_bounds, &section_bounds[i]);

                --nresplit;
                std::cerr<<nresplit<<" more resplits allowed.\n";

                continue;
            } else {
                bombout_resplit(&this->params);
            }
        }
        std::cerr<<this->get_ngals_with_pixels()<<
                " galaxies in this section.\n";

    }


}

static void make_cutouts(ConfigFile *params) 
{
    // Read the coadd catalog.
    CoaddCatalog coaddcat(*params);
    coaddcat.read();
    dbg<<"Constructed coaddcat\n";

    std::string imlist_file=(*params)["coadd_srclist"];
    CutoutMaker maker(&coaddcat, imlist_file, params);

    maker.print_stats();
    maker.make_cutouts();

}

int main(int argc, char **argv)
{
    ConfigFile params;
    if (BasicSetup(argc,argv,params,"make_cutouts")) return EXIT_FAILURE;

    make_cutouts(&params);
}
