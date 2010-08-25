#include "cwl.h"
#include "../src/StarFinder.h"
#include "../src/dbg.h"

#include "../src/fp.h"

#if defined(__GNUC__) && defined(OPENMP_LINK)
//__thread std::ostream* dbgout=&std::cout;
//__thread bool XDEBUG=true;
__thread std::ostream* dbgout=0;
__thread bool XDEBUG=false;
#else
//std::ostream* dbgout=&std::cout;
//bool XDEBUG=true;
std::ostream* dbgout=0;
bool XDEBUG=false
#endif


CWL::CWL(std::string config_file) throw (const char*) {
  this->load_config(config_file);
  this->load_fitsparams();
}

void CWL::load_config_images_catalog(
    std::string config_file,
    std::string image_file,
    std::string cat_file) throw (const char*) {

  this->load_config(config_file);
  this->load_fitsparams();

  this->load_images(image_file);
  this->load_catalog(cat_file);
}

void CWL::load_config(std::string file) throw (const char*) {
  try {
    this->params.load(file);
  }
  PY_CATCHALL 
}

void CWL::load_fitsparams() {
    std::string fp((const char*)fitsparams_config,fitsparams_config_len);
    std::istringstream is(fp);
    this->params.read(is);
}

void CWL::print_config() {
  std::cout<<this->params;
}
void CWL::set_verbosity(bool verbosity) {
  this->params["verbose"] = verbosity;
  if (params.read<int>("verbose") > 1) XDEBUG = true;
  dbgout = &std::cout;
  dbgout->setf(std::ios_base::unitbuf);
}

void CWL::load_images(std::string file) throw (const char*) {

  if (this->params.size() == 0) {
    throw "load a config file";
  }

  this->params["image_file"] = file;
  this->params["image_hdu"] = 2;

  this->params["badpix_file"] = file;
  this->params["badpix_hdu"] = 3;

  this->params["weight_file"] = file;
  this->params["weight_hdu"] = 4;


  try {
    this->image.reset( new Image<double>(this->params, this->weight_image) );
    this->load_trans(file);
  } catch(CCfits::FITS::CantOpen& e) {
    this->err<<"CCfits error: file not found: "<<file;
    throw this->err.str().c_str();
  }
  PY_CATCHALL 
}

void CWL::load_catalog(std::string file) throw (const char*) {
  if (this->params.size() == 0) {
    throw "load a config file";
  }

  this->params["cat_file"] = file;
  this->params["cat_hdu"] = 3;

  try {
    this->cat.reset( new InputCatalog(this->params) );
    this->cat->read();
  } catch(CCfits::FITS::CantOpen& e) {
    this->err<<"CCfits error: file not found: "<<file;
    throw this->err.str().c_str();
  }
  PY_CATCHALL 

}

void CWL::load_trans(std::string file) throw (const char*) {

  this->params["dist_file"] = file;
  this->params["dist_hdu"] = 2;
  try {
    trans.initFromParams(this->params);
  } catch(CCfits::FITS::CantOpen& e) {
    this->err<<"CCfits error: file not found: "<<file;
    throw this->err.str().c_str();
  }
  PY_CATCHALL 
}


void CWL::find_stars(std::string outfile) throw (const char*) {

  std::stringstream err;

  if (!this->image->loaded()) {
    throw "image not loaded";
  }
  if (!this->weight_image->loaded()) {
    throw "weight image not loaded";
  }
  if (this->cat->size() == 0) {
    throw "catalog not loaded";
  }

  this->log.reset(
      new FindStarsLog(this->params,
                       this->logfile,
                       outfile));
  // initialize the star cat from cat
  this->starcat.reset( new StarCatalog(*this->cat.get(),this->params) );

  // Update the sizes to more robust values
  this->starcat->calculateSizes(
      *this->image.get(),
      this->weight_image.get(),
      this->trans);

  try {
    
    this->starcat->findStars(  
        static_cast<FindStarsLog&>(*this->log)
    );
    
  } catch(StarFinderException& e) {
    // Need to catch this here, so we can write the output file
    // with the sizes, even though we haven't figured out which 
    // objects are stars.
    this->write_starcat(outfile);
    err<<"Caught StarFinderException: "<<e.what();
    throw err.str().c_str();
  } catch (...) {
    this->write_starcat(outfile);
    err<<"Caught unknown exception\n";
    throw err.str().c_str();
  }

  this->write_starcat(outfile);

}

void CWL::write_starcat(std::string file) throw (const char*) {
  try {
    this->starcat->writeFits(file);
  }
  PY_CATCHALL
}
