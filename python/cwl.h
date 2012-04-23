#include "../src/ConfigFile.h"
#include "../src/Image.h"
#include "../src/InputCatalog.h"
#include "../src/Transformation.h"
#include "../src/StarCatalog.h"
#include "../src/PsfCatalog.h"
#include "../src/FittedPsf.h"
#include "../src/ShearCatalog.h"
#include "../src/Pixel.h"
#include "../src/MeasureShearAlgo.h"
#include "errors.h"

#include <Python.h>
#include "numpy/arrayobject.h" 

// a wrapper class for the wl classes and functions.  This is designed to
// be simple to call from python and supports less general functionality.
// Another principle is encapsulation to simplify calls

// Also, python will abort if an exception is thrown and not caught.  for now
// I'm converting all exceptions to const char* since these are easy to work
// with in SWIG, but I plan to configure the SWIG wrapper to catch other types
// later.

using std::string;

// A class to make it easy to run on an input image and PSF.
class WLQuick {
    public:

        WLQuick(int order_psf, int order_shear_gal) 
                : psf_order(order_psf), shear_gal_order(order_shear_gal),
                  psf(NULL), image(NULL), 
                  psf_shapelets(order_psf,1.0),
                  gal_shapelets(order_shear_gal,1.0) {
            import_array();

            // set all the parametes that are not entered
            // these are defaults ala wl.config
            this->shear_aperture = 4.0; //sigma?
            // this will get set below
            this->shear_max_aperture = 44.0; // pixels?  We use 12'' so 44 pix
            this->shear_gal_order2 = 20;
            this->shear_min_gal_order = 4;
            this->shear_min_gal_size = 0.5; // in multiples of psf size
            this->shear_base_order_on_nu = false;


            this->minfpsf = 2.0;
            this->maxfpsf = 2.0;
            this->shear_fix_centroid=false;

            this->shear_do_force_sigma=false;
            this->shear_force_sigma_val=0;

            this->native_only=false;

            // gain=0 means ignore gain and don't add im/gain to variance
            this->gain=0;
            this->xoffset=0;
            this->yoffset=0;

            this->psf_sigma=-9999;
        };
        ~WLQuick() {
            if (this->psf) delete this->psf;
            if (this->image) delete this->image;

            this->psf=NULL;
            this->image=NULL;
        }

        void set_image(PyObject* image_obj, 
                       double rcen, 
                       double ccen, 
                       double sky, 
                       double skyvar,
                       double aperture)  // this is max aperture pixels
                       throw (const char*) {
            check_numpy_image(image_obj);
            if (this->image) { delete this->image; this->image=NULL; }
            this->image = copy_numpy_image(image_obj);

            std::complex<double> tcen(rcen,ccen);
            this->image_cen = tcen;

            this->image_sky=sky;
            this->image_skyvar=skyvar;

            this->shear_max_aperture = aperture;

        }
        void set_psf(PyObject* psf_obj, 
                     double rcen, 
                     double ccen,
                     double sky,
                     double aperture)  // max aperture
                     throw (const char*) {
            check_numpy_image(psf_obj);
            if (this->psf) { delete this->psf; this->psf=NULL; }

            this->psf = copy_numpy_image(psf_obj);

            std::complex<double> tcen(rcen,ccen);
            this->psf_cen = tcen;
            this->psf_sky = sky;
            this->psf_max_aperture=aperture;
        }


        long calculate_psf_sigma(double guess) throw (const char*) {
            if (!this->psf) {
                throw "psf is not set";
            }
            this->psf_sigma=guess;

            Image<double>* weightIm=NULL;
            Transformation trans; // defaults to identity

            bool shouldUseShapeletSigma=0;
            double skyvar=1;
            long flags=0;
            calculateSigma(
                    this->psf_sigma,
                    *this->psf, this->psf_cen, this->psf_sky, 
                    skyvar, this->gain, weightIm,
                    trans, this->psf_max_aperture,
                    this->xoffset, this->yoffset, flags, 
                    shouldUseShapeletSigma);
            return flags;
        }
        long calculate_sigma0(double guess) throw (const char*) {
            if (!this->image) {
                throw "image is not set";
            }
            this->sigma0=guess;

            Image<double>* weightIm=NULL;
            Transformation trans; // defaults to identity

            bool shouldUseShapeletSigma=0;
            double skyvar=1;
            long flags=0;
            calculateSigma(
                    this->sigma0,
                    *this->image, this->image_cen, this->image_sky, 
                    this->image_skyvar, this->gain, weightIm,
                    trans, this->shear_max_aperture,
                    this->xoffset, this->yoffset, flags, 
                    shouldUseShapeletSigma);
            return flags;
        }


        long calculate_psf_shapelets() throw (const char*) {
            if (!this->psf) {
                throw "psf is not set";
            }
            if (this->psf_sigma <= 0) {
                throw "run calculate_psf_sigma first";
            }

            Image<double>* weightIm=NULL;
            Transformation trans; // defaults to identity
            double skyvar=1;
            long flags=0;
            bool fix_cen=false;
            //ConfigFile fake_config;
            //PsfLog log(fake_config,"","");
            PsfLog log;
            int maxm=psf_order;

            measureSinglePsf( 
                    // Input data:
                    this->psf_cen, *this->psf, this->psf_sky, 
                    trans, 
                    // Noise values:
                    skyvar, this->gain, weightIm,
                    // Parameters:
                    this->psf_sigma, this->shear_max_aperture, 
                    this->psf_order, maxm,
                    fix_cen, this->xoffset, this->yoffset,
                    // Log information
                    log,
                    // Ouput value:
                    this->psf_shapelets, this->psf_nu, flags);
            return flags;
        }

        long calculate_shear() throw (const char*) {
            if (!this->image || !this->psf) {
                throw "set both psf and image first";
            }

            Image<double>* weightIm=NULL;
            Transformation trans; // defaults to identity

            std::vector<PixelList> pix(1);
            //ConfigFile fake_config;
            //ShearLog log(fake_config); 
            ShearLog log;
            // we are forced to make a copy of our psf shapelets
            // because the function takes a vector
            std::vector<BVec> vpsf(1, this->psf_shapelets);
            long flags=0;

            // on return will hold order we used
            this->meas_gal_order = this->shear_gal_order;
            getPixList(*this->image,pix[0],this->image_cen, 
                       this->image_sky,this->image_skyvar,
                       this->gain,weightIm, trans,
                       this->shear_max_aperture,
                       xoffset,yoffset,flags);

            if (flags != 0) {
                return flags;
            }

            int maxm=shear_gal_order;

            measureSingleShear(
                    // Input data:
                    pix, vpsf,
                    // Parameters:
                    this->shear_aperture, this->shear_max_aperture,
                    meas_gal_order, shear_gal_order2, maxm,
                    shear_min_gal_order, shear_base_order_on_nu,
                    minfpsf, maxfpsf, shear_min_gal_size, shear_fix_centroid,
                    shear_do_force_sigma,shear_force_sigma_val, 
                    native_only,
                    // Log information
                    log,
                    // Ouput values:
                    this->gal_shapelets, 
                    this->shear, this->shear_cov, this->shear_nu, flags);
            return flags;
        }


        double get_image_val(int row, int col) throw (const char*) {
            if (!this->image) { throw "image is NULL"; }
            return (*this->image)(row,col);
        }
        double get_psf_val(int row, int col) throw (const char*) {
            if (!this->psf) { throw "psf is NULL"; }
            return (*this->psf)(row,col);
        }

        double get_psf_sigma() {
            return this->psf_sigma;
        }
        double get_psf_nu() {
            return this->psf_nu;
        }
        void print_psf_shapelets() {
            std::cout<<this->psf_shapelets<<"\n";
        }

        double get_sigma0() {
            return sigma0;
        }
        double get_sigma() {
            return gal_shapelets.getSigma();
        }
        double get_nu() {
            return shear_nu;
        }
        double get_cov11() {
            return shear_cov(0,0);
        }
        double get_cov12() {
            return shear_cov(0,1);
        }
        double get_cov22() {
            return shear_cov(1,1);
        }
        double get_shear1() {
            return real(this->shear);
        }
        double get_shear2() {
            return imag(this->shear);
        }
        double get_e1() {
            return this->gal_shapelets(3)/gal_shapelets(0)*sqrt(2);
        }
        double get_e2() {
            return this->gal_shapelets(4)/gal_shapelets(0)*sqrt(2);
        }






        // python stuff
        Image<double>* copy_numpy_image(PyObject* npyim) {
            npy_intp* dims = PyArray_DIMS(npyim);
            Image<double> *im = new Image<double>(dims[0], dims[1]);

            double *pdata;
            for (npy_intp row=0; row<dims[0]; row++) {
                for (npy_intp col=0; col<dims[1]; col++) {
                    pdata = (double*) PyArray_GETPTR2(npyim,row,col);
                    (*im)(row,col) = *pdata;
                }
            }
            return im;
        }
        void check_numpy_image(PyObject* image_obj) throw (const char*) {
            if (!PyArray_Check(image_obj) 
                    || NPY_DOUBLE != PyArray_TYPE(image_obj)
                    || 2 != PyArray_NDIM(image_obj)) {
                throw "image input must be a 2D double PyArrayObject";
            }
        }
        void hello() {
            std::cout<<"hello world\n";
        }

    private:

        double sigma0;

        int psf_order;
        Image<double> *psf;
        Position psf_cen;
        double psf_max_aperture; // in pixels?
        double psf_sigma;
        double psf_sky;
        BVec psf_shapelets; // must be initialized on construction
        double psf_nu;

        Image<double> *image;
        Position image_cen;

        std::complex<double> shear;
        DSmallMatrix22 shear_cov;
        double shear_nu;

        double image_sky;
        double image_skyvar;

        int shear_gal_order;
        int shear_gal_order2;
        int shear_min_gal_order;
        int meas_gal_order;  // will hold what we actually used
        bool shear_fix_centroid;
        bool shear_base_order_on_nu;
        double shear_aperture; // in sigma
        double shear_max_aperture; // in pixels?
        double shear_min_gal_size;  // multiples of psf size

        bool native_only;
        bool shear_do_force_sigma;
        double shear_force_sigma_val;
        double minfpsf;
        double maxfpsf;

        BVec gal_shapelets; // must be initialized on construction

        double gain;
        double xoffset;
        double yoffset;
};


class CWL {
  public:
    CWL() {};

    CWL(string config_file) throw (const char*);

    // file i/o methods
    void load_config(string file) throw (const char*);
    void load_fitsparams();

    // note only strings supported for now
    void set_param(string key, string value) throw (const char*);
    // this sets all logs to the same name
    void set_log(string logfile);

    void load_config_images_catalog(
        string config_file,
        string image_file,
        string cat_file) throw (const char*);

    void load_images(string file) throw (const char*);
    void load_catalog(string file) throw (const char*);
    void load_trans(string file)  throw (const char*);

    void write_starcat(string file, bool flush_log=true) throw (const char*);
    void load_starcat(string file) throw (const char*);

    void write_psfcat(string file, bool flush_log=true) throw (const char*);
    void load_psfcat(string file) throw (const char*);
    void write_fitpsf(string file) throw (const char*);
    void load_fitpsf(string file) throw (const char*);

    void write_shearcat(string file, bool flush_log=true) throw (const char*);
    void load_shearcat(string file) throw (const char*);

    void split_starcat(string file1, string file2) throw (const char*);

    //
    // processing methods
    //
    
    // because of the log system, which we rely on to write 
    // status entries to the header, these must be all-in-one
    // calculating and writing the data to file.  We can't
    // factor out find_stars() and then run write(file) because
    // the log must know the file name
    void find_stars(string outfile) throw (const char*);

    void measure_psf(
        string psf_file, 
        string fitpsf_file) throw (const char*);
    void measure_shear(string shear_file) throw (const char*);

    void print_config();
    void set_verbose(int verbosity);

    /*
    template<typename T, size_t n> size_t array_size(const T (&)[n]) {
        return n;
    }
    */

  private:
    ConfigFile params;

    std::auto_ptr<Image<double> > image;
    std::auto_ptr<Image<double> > weight_image;
    std::auto_ptr<Transformation> trans;

    std::auto_ptr<InputCatalog>   cat;

    std::auto_ptr<StarCatalog>    starcat;
    std::auto_ptr<PsfCatalog>     psfcat;
    std::auto_ptr<FittedPsf>      fitpsf;
    std::auto_ptr<ShearCatalog>   shearcat;


    // send logs to stdout unless logfile gets set
    string stars_logfile;
    std::auto_ptr<FindStarsLog> stars_log;

    string psf_logfile;
    std::auto_ptr<PsfLog> psf_log;

    string shear_logfile;
    std::auto_ptr<ShearLog> shear_log;

};
