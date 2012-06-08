#include <cstdlib>
#include "../src/ConfigFile.h"
#include "../src/Image.h"
#include "../src/InputCatalog.h"
#include "../src/Transformation.h"
#include "../src/StarCatalog.h"
#include "../src/PsfCatalog.h"
#include "../src/FittedPsf.h"
#include "../src/ShearCatalog.h"
#include "../src/Pixel.h"
#include "../src/Ellipse.h"
#include "../src/MeasureShearAlgo.h"
#include "errors.h"

#include <Python.h>
#include "numpy/arrayobject.h" 

// Python will abort if an exception is thrown and not caught.  for now I'm
// converting all exceptions to const char* since these are easy to work with
// in SWIG, but I plan to configure the SWIG wrapper to catch other types
// later.

using std::string;


class WLShear;

// we hide most data members because we can't access most of them
// from python anyway and don't want the swig wrapper to blow up
// Use a friend class instead
class WLObject {
    public:

        WLObject(PyObject* image_obj, 
                 double rcen, 
                 double ccen, 
                 double sky, 
                 double skyvar,
                 double max_aperture, // this is max aperture pixels
                 double sigma_guess)  // guess at sigma0
                               throw (const char*) {

            import_array();
            this->flags=0;
            this->image = this->copy_numpy_image(image_obj);

            std::complex<double> tcen(rcen,ccen);
            this->cen = tcen;

            this->sky=sky;
            this->skyvar=skyvar;
            this->weight_image=NULL;
            // gain=0 means ignore gain and don't add im/gain to variance
            this->gain=0;
            this->sigma0_guess=sigma_guess;

            this->max_aperture = max_aperture;

            // also sets flags
            this->get_pixel_list();

            this->sigma0=-9999;
            if (this->flags==0) {
                this->calculate_sigma0();
            }

        }
        ~WLObject() {
            if (this->image) delete this->image;
            this->image=NULL;
        }

        double get_sigma0() {
            return this->sigma0;
        }
        long get_flags() {
            return this->flags;
        }
        /*
        double get_max_aperture() {
            return this->max_aperture;
        }
        */

        friend class WLShear;

    private:
        Position cen;
        double sky;
        double skyvar;
        double gain;
        double max_aperture;
        double sigma0_guess;
        double sigma0;
        long flags;
        Image<double> *image;
        Image<double>* weight_image;
        Transformation trans; // defaults to identity

        // we have to use a vector here because of 
        // measureSingleShear
        std::vector<PixelList> vpix;



        void get_pixel_list() {
            this->vpix.resize(1);
            try {
                ConfigFile params;
                // in Pixel.cpp
                GetPixList(*this->image, 
                           this->vpix[0], 
                           this->cen, 
                           this->sky, 
                           this->skyvar, 
                           this->weight_image, 
                           this->trans, 
                           this->max_aperture, 
                           params,
                           this->flags);

                /*
                getPixList(*this->image, 
                           this->vpix[0], this->cen, 
                           this->sky, this->skyvar, 
                           this->gain, this->weight_image, this->trans, 
                           this->max_aperture, 0.0, 0.0, this->flags);
                */
            } catch (RangeException& e) {
                throw "distortion range error: \n";
                xdbg<<"center = "<<this->cen<<", b = "<<e.getBounds()<<std::endl;
                dbg<<"FLAG TRANSFORM_EXCEPTION\n";
                std::cerr<<"distortion range error\n";
                this->flags |= TRANSFORM_EXCEPTION;
            }
        }
        void calculate_sigma0() {
            Ellipse ell;
            ell.fixGam();
            ell.peakCentroid(this->vpix[0],this->max_aperture/3.);
            // implementation in CrudeMeasure.cpp
            ell.crudeMeasure(this->vpix[0],this->sigma0_guess);

            double mu = ell.getMu();
            this->sigma0 = this->sigma0_guess*exp(mu);
        }



        // python stuff
        Image<double>* copy_numpy_image(PyObject* npyim) throw (const char*) {

            this->check_numpy_image(npyim);
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

};

class WLShear {
    public:

        /*
        WLShear(const WLObject &obj, const WLObject &psfobj,
                int order_psf, int order_gal) 
                    : psf_order(order_psf), psf_maxm(order_psf), 
                      gal_order(order_gal), gal_maxm(order_gal),
                      gal_order2(20), 
                      min_gal_order(4),
                      min_gal_size(0.5),  // in multiple of psf size
                      shear_aperture(5.0), // in sigma
                      shear_fix_centroid(false),
                      minfpsf(2), maxfpsf(2), 
                      base_order_on_nu(false),
                      native_only(false),
                      shear_force_sigma(false),
                      shear_force_sigma_val(0),
                      shear(-9999,-9999),
                      psf_shapelets(order_psf,1.0),
                      gal_shapelets(order_gal,1.0),
                      shear_nu(0),
                      flags(0), fix_psf_cen(false) {
                      */
        WLShear(std::string config_fname,
                const WLObject &obj, 
                const WLObject &psfobj,
                int order_psf, 
                int order_gal, 
                double shear_ap_nsig) throw (const char*) // in number of sigma
                    : psf_order(order_psf), 
                      gal_order(order_gal), 
                      shear_aperture_nsig(shear_ap_nsig),
                      psf_shapelets(order_psf,1.0),
                      gal_shapelets(order_gal,1.0),
                      shear(-9999,-9999),
                      shear_nu(0),
                      flags(0) {

            import_array();
            this->load_config(config_fname);
            this->set_extra_config_pars(obj);

            shear_cov(0,0)=9999;
            shear_cov(1,1)=9999;
            this->measure_psf(psfobj);
            if (this->flags == 0) {
                this->measure_shear(obj);
            }
        }


        long get_flags() {
            return this->flags;
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
        double get_prepsf_sigma() {
            return gal_shapelets.getSigma();
        }

    private:
        ConfigFile params;

        int psf_order;
        //int psf_maxm;  // will equal psf_order
        int gal_order;
        //int gal_order2;
        //int gal_maxm;  // will equal gal_order
        //int min_gal_order;
        //double min_gal_size;// in multiple of psf size
        //bool shear_fix_centroid;
        double shear_aperture_nsig; // in sigma
        //bool shear_force_sigma;
        //double shear_force_sigma_val;
        //double minfpsf;
        //double maxfpsf;
        //bool base_order_on_nu;
        //bool native_only;
        //bool fix_psf_cen;
        BVec psf_shapelets; // must be initialized on construction
        Position psf_cen;   // we will allow cen to change from input
        BVec gal_shapelets; // must be initialized on construction

        // these will be made visible
        std::complex<double> shear;
        DSmallMatrix22 shear_cov;
        double shear_nu;
        long flags;


        void load_config(std::string cf) throw (const char*) {
            try {
                this->params.load(cf);
            } catch (FileNotFoundException& e) {
                std::string err;
                err = "file not found: " + cf + "\n";
                throw err.c_str();
            }
        }
        void set_extra_config_pars(const WLObject &obj) {
            params["ignore_edges"] = true;
            params["shear_aperture"] = this->shear_aperture_nsig; 
            params["shear_max_aperture"] = obj.max_aperture;
            params["shear_gal_order"] = this->gal_order;
            params["shear_gal_order2"] = 20;
            params["shear_maxm"] = this->gal_order;
            params["psf_maxm"] = this->psf_order;
            params["shear_min_gal_order"] = 4;
            params["shear_f_psf"] = 2; // min f psf
            params["shear_max_f_psf"] = 2;
            params["shear_min_gal_size"] = 0.5;


            /*
            this->params["psf_order"] = 
                    : psf_order(order_psf), psf_maxm(order_psf), 
                      gal_order(order_gal), gal_maxm(order_gal),
                      gal_order2(20), 
                      min_gal_order(4),
                      min_gal_size(0.5),  // in multiple of psf size
                      shear_aperture(5.0), // in sigma
                      shear_fix_centroid(false),
                      minfpsf(2), maxfpsf(2), 
                      base_order_on_nu(false),
                      native_only(false),
                      shear_force_sigma(false),
                      shear_force_sigma_val(0),
                      shear(-9999,-9999),
                      psf_shapelets(order_psf,1.0),
                      gal_shapelets(order_gal,1.0),
                      shear_nu(0),
                      flags(0), fix_psf_cen(false) {

                      */
        }


        void measure_shear(const WLObject &obj) {
            // we are forced to make a copy of our psf shapelets
            // because the function takes a vector
            std::vector<BVec> vpsf(1, this->psf_shapelets);
            ShearLog log;
            // in MeasureShearAlgo.cpp
            /*
            ConfigFile params;
            params["shear_aperture"] = this->shear_aperture;
            params["shear_max_aperture"] = obj.max_aperture;
            params["shear_gal_order"] = this->gal_order;
            params["shear_gal_order2"] = this->gal_order2;
            params["shear_maxm"] = this->gal_maxm;
            params["shear_min_gal_order"] = this->min_gal_order;
            params["shear_f_psf"] = this->minfpsf;
            params["shear_max_f_psf"] = this->maxfpsf;
            params["shear_min_gal_size"] = this->min_gal_size;
            params["shear_fix_centroid"] = this->shear_fix_centroid;

            params["shear_use_fake_pixels"] = true;
            if (this->shear_force_sigma)
                params["shear_force_sigma"] = this->shear_force_sigma_val;
            */
            int galorder_final=0; // returned
            MeasureSingleShear(
                    obj.vpix,       // const ref
                    vpsf,           // const ref
                    galorder_final, // input ref, value not used, modified
                    this->params,         // const ref
                    log,            // input ref, value not used, modified
                    this->gal_shapelets, // input ref, value not used, modified
                    this->shear,         // input ref, value not used, modified
                    this->shear_cov, // input ref, value not used, modified
                    this->shear_nu, // input ref, value not used, modified
                    this->flags); // input ref, value not used, modified

            /*
            measureSingleShear(
                    // Input data:
                    obj.vpix, vpsf,
                    // Parameters:
                    this->shear_aperture, 
                    obj.max_aperture,
                    this->gal_order, 
                    this->gal_order2, 
                    this->gal_maxm,
                    this->min_gal_order, 
                    this->base_order_on_nu,
                    this->minfpsf, 
                    this->maxfpsf, 
                    this->min_gal_size, 
                    this->shear_fix_centroid,
                    this->shear_force_sigma,
                    this->shear_force_sigma_val, 
                    this->native_only,
                    // Log information
                    log,
                    // Ouput values:
                    this->gal_shapelets, 
                    this->shear, 
                    this->shear_cov, 
                    this->shear_nu, 
                    this->flags);
            */
            return;
        }
        void measure_psf(const WLObject &psfobj) {
            // we have to repeat this because this one has fixMu
            this->psf_cen = psfobj.cen;

            Ellipse ell;
            ell.fixGam();
            ell.fixMu();
            bool fix_psf_cen = this->params.read("psf_fix_centroid",false);
            if (fix_psf_cen) ell.fixCen();
            else {
                ell.peakCentroid(psfobj.vpix[0],psfobj.max_aperture/3.);
                ell.crudeMeasure(psfobj.vpix[0],psfobj.sigma0);
            }

            // First make sure it is centered.
            int psf_maxm = this->params["psf_maxm"];
            if (!(ell.measure(psfobj.vpix,this->psf_order,this->psf_order+4,
                              psf_maxm,psfobj.sigma0,
                              this->flags,1.e-4))) {
                dbg<<"initial measure fail FLAG MEASURE_PSF_FAILED\n";
                this->flags |= MEASURE_PSF_FAILED;
                return;
            } 

            DMatrix cov(int(this->psf_shapelets.size()),
                        int(this->psf_shapelets.size()));
            // Then measure it with the accurate center.
            this->psf_shapelets.setSigma(psfobj.sigma0);
            if (!ell.measureShapelet(psfobj.vpix,
                                     this->psf_shapelets,
                                     this->psf_order,this->psf_order+4,
                                     psf_maxm,&cov)) {
                dbg<<"psf shapelet meas fail: FLAG MEASURE_PSF_FAILED\n";
                this->flags |= MEASURE_PSF_FAILED;
                return;
            }

            // also sets flags
            if (!this->check_psf(psfobj, ell))
                return;

            this->psf_shapelets.normalize();  // Divide by (0,0) element
            if (!fix_psf_cen) this->psf_cen += ell.getCen();

        }
        int check_psf(const WLObject &psfobj, Ellipse &ell) {
            // Calculate the 0th order shapelet to make sure the flux is consistent 
            // fail if a factor of 3 diff
            BVec flux(0,psfobj.sigma0);
            DMatrix fluxCov(1,1);
            if (!ell.measureShapelet(psfobj.vpix,flux,0,0,0,&fluxCov)) {
                dbg<<"psf flux meas fail: FLAG PSF_BAD_FLUX\n";
                this->flags |= PSF_BAD_FLUX;
                return 0;
            }
            double nu = flux(0) / std::sqrt(fluxCov(0,0));
            if (!(flux(0) > 0.0 &&
                        this->psf_shapelets(0) >= flux(0)/3. &&
                        this->psf_shapelets(0) <= flux(0)*3.)) {
                dbg<<"Bad flux value: \n";
                dbg<<"flux = "<<flux(0)<<std::endl;
                dbg<<"psf = "<<this->psf_shapelets.vec()<<std::endl;
                dbg<<"FLAG PSF_BAD_FLUX\n";
                this->flags |= PSF_BAD_FLUX;
                return 0;
            }
            return 1;
        }
};


/*
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

    //template<typename T, size_t n> size_t array_size(const T (&)[n]) {
    //    return n;
    //}

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
*/
