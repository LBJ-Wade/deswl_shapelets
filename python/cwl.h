#include "../src/ConfigFile.h"
#include "../src/Image.h"
#include "../src/InputCatalog.h"
#include "../src/Transformation.h"
#include "../src/StarCatalog.h"
#include "errors.h"

// a wrapper class for the wl classes and functions.  This is designed to
// be simple to call from python and supports less general functionality.
class CWL {
  public:
    CWL() {};

    CWL(std::string config_file) throw (const char*);

    void load_config(std::string file) throw (const char*);
    void load_fitsparams();

    void load_config_images_catalog(
        std::string config_file,
        std::string image_file,
        std::string cat_file) throw (const char*);

    void load_images(std::string file) throw (const char*);
    void load_catalog(std::string file) throw (const char*);
    void load_trans(std::string file)  throw (const char*);

    void find_stars(std::string outfile) throw (const char*);

    void write_starcat(std::string file) throw (const char*);

    void print_config();
    void set_verbosity(bool verbosity);

  private:
    ConfigFile params;

    std::auto_ptr<Image<double> > image;
    std::auto_ptr<Image<double> > weight_image;
    std::auto_ptr<InputCatalog> cat;
    Transformation trans;

    std::auto_ptr<StarCatalog> starcat;

    // send logs to stdout
    std::string logfile;
    std::auto_ptr<Log> log;

    std::stringstream err;
};
