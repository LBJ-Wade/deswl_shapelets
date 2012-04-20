#include "../src/ConfigFile.h"
#include "../src/Image.h"
#include "../src/InputCatalog.h"
#include "../src/Transformation.h"
#include "../src/StarCatalog.h"
#include "../src/PsfCatalog.h"
#include "../src/FittedPsf.h"
#include "../src/ShearCatalog.h"
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

        WLQuick() : psf(NULL), image(NULL) {
            import_array();
            psf=NULL;
            image=NULL;
        };
        ~WLQuick() {
            if (psf) delete psf;
            if (image) delete image;

            psf=NULL;
            image=NULL;
        }

        void set_image(PyObject* image_obj, double rcen, double ccen) 
                       throw (const char*) {
            check_numpy_image(image_obj);
            if (image) { delete image; image=NULL; }
            image = copy_numpy_image(image_obj);

            image_cen = 
        }
        void set_psf(PyObject* psf_obj, double rcen, double ccen) 
                     throw (const char*) {
            check_numpy_image(psf_obj);
            if (psf) { delete psf; psf=NULL; }

            psf = copy_numpy_image(psf_obj);

            std::complex<double> tcen(rcen,ccen);
            psf_cen = tcen;
        }

        double get_image_val(int row, int col) throw (const char*) {
            if (!image) { throw "image is NULL"; }
            return (*image)(row,col);
        }
        double get_psf_val(int row, int col) throw (const char*) {
            if (!psf) { throw "psf is NULL"; }
            return (*psf)(row,col);
        }

        void calculate_psf_sigma() {
            double sigma=-9999;
            /*
            calculateSigma1(
                    sigma,
                    *psf, pos, sky, noise, gain, weightIm,
                    trans, psfAp, xOffset, yOffset, flag, shouldUseShapeletSigma);
                    */
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
        Image<double> *psf;
        Position psf_cen;
        Image<double> *image;
        Position image_cen;
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
