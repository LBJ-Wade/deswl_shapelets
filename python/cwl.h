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

using std::string;

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
    string stars_logfile;
    std::auto_ptr<FindStarsLog> stars_log;

    string psf_logfile;
    std::auto_ptr<PsfLog> psf_log;

    string shear_logfile;
    std::auto_ptr<ShearLog> shear_log;

};
