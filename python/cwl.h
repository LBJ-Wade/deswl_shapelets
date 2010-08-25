#include "../src/ConfigFile.h"
#include "../src/Image.h"
#include "../src/InputCatalog.h"
#include "../src/Transformation.h"
#include "../src/StarCatalog.h"
#include "../src/PsfCatalog.h"
#include "../src/FittedPsf.h"
#include "../src/ShearCatalog.h"
#include "errors.h"

// a wrapper class for the wl classes and functions.  This is designed to
// be simple to call from python and supports less general functionality.
// Another principle is encapsulation to simplify calls

// Also, python will abort if an exception is thrown and not caught.  for now
// I'm converting all exceptions to const char* since these are easy to work
// with in SWIG, but I plan to configure the SWIG wrapper to catch other types
// later.

class CWL {
  public:
    CWL() {};

    CWL(std::string config_file) throw (const char*);

    // file i/o methods
    void load_config(std::string file) throw (const char*);
    void load_fitsparams();

    void load_config_images_catalog(
        std::string config_file,
        std::string image_file,
        std::string cat_file) throw (const char*);

    void load_images(std::string file) throw (const char*);
    void load_catalog(std::string file) throw (const char*);
    void load_trans(std::string file)  throw (const char*);

    void write_starcat(std::string file, bool flush_log=true) throw (const char*);
    void load_starcat(std::string file) throw (const char*);

    void write_psfcat(std::string file, bool flush_log=true) throw (const char*);
    void load_psfcat(std::string file) throw (const char*);
    void write_fitpsf(std::string file) throw (const char*);
    void load_fitpsf(std::string file) throw (const char*);

    void write_shearcat(std::string file, bool flush_log=true) throw (const char*);
    void load_shearcat(std::string file) throw (const char*);

    //
    // processing methods
    //
    
    // because of the log system, which we rely on to write 
    // status entries to the header, these must be all-in-one
    // calculating and writing the data to file.  We can't
    // factor out find_stars() and then run write(file) because
    // the log must know the file name
    void find_stars(std::string outfile) throw (const char*);

    void measure_psf(
        std::string psf_file, 
        std::string fitpsf_file) throw (const char*);
    void measure_shear(std::string shear_file) throw (const char*);

    void print_config();
    void set_verbosity(bool verbosity);

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
    std::string stars_logfile;
    std::auto_ptr<FindStarsLog> stars_log;

    std::string psf_logfile;
    std::auto_ptr<PsfLog> psf_log;

    std::string shear_logfile;
    std::auto_ptr<ShearLog> shear_log;

};
