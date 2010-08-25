#include <CCfits/CCfits>
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

  std::stringstream err;
  if (this->params.size() == 0) {
    throw "load a config file";
  }

  this->params["image_file"] = file;
  this->params["image_hdu"] = 2;

  this->params["badpix_file"] = file;
  this->params["badpix_hdu"] = 3;

  this->params["weight_file"] = file;
  this->params["weight_hdu"] = 4;

  // not using CCfits for images yet, so must catch ccfits error 
  try {
    this->image.reset( new Image<double>(this->params, this->weight_image) );
    this->load_trans(file);
  } catch(CCfits::FITS::CantOpen& e) {
    err<<"CCfits error: file not found: "<<file;
    throw err.str().c_str();
  }
  PY_CATCHALL 
}

void CWL::load_catalog(std::string file) throw (const char*) {

  std::stringstream err;
  if (this->params.size() == 0) {
    throw "load a config file";
  }

  this->params["cat_file"] = file;
  this->params["cat_hdu"] = 3;

  try {
    this->cat.reset( new InputCatalog(this->params) );
    this->cat->read();
  }
  PY_CATCHALL 

}

void CWL::load_trans(std::string file) throw (const char*) {

  std::stringstream err;
  this->params["dist_file"] = file;
  this->params["dist_hdu"] = 2;
  // not using CCfits for trans yet, so must catch ccfits error 
  try {
    this->trans.reset(
        new Transformation(this->params)
    );
  } catch(CCfits::FITS::CantOpen& e) {
    err<<"CCfits error: file not found: "<<file;
    throw err.str().c_str();
  }
  PY_CATCHALL 
}

void CWL::write_starcat(std::string file, bool flush_log) throw (const char*) {
  try {
    this->starcat->writeFits(file);
    if (flush_log) {
      // to flush the log we must delete it
      this->stars_log.reset();
    }
  }
  PY_CATCHALL
}
void CWL::load_starcat(std::string file) throw (const char*) {
  this->starcat.reset( new StarCatalog( this->params ) );
  try {
    this->starcat->readFits(file);
  }
  PY_CATCHALL 
}

void CWL::write_psfcat(std::string file, bool flush_log) throw (const char*) {
  try {
    this->psfcat->writeFits(file);
    if (flush_log) {
      // to flush the log we must delete it
      this->psf_log.reset();
    }
  }
  PY_CATCHALL
}
void CWL::load_psfcat(std::string file) throw (const char*) {
  this->psfcat.reset( new PsfCatalog( this->params ) );
  try {
    this->psfcat->readFits(file);
  }
  PY_CATCHALL 
}


void CWL::write_fitpsf(std::string file) throw (const char*) {
  try {
    this->fitpsf->writeFits(file);
  }
  PY_CATCHALL
}
void CWL::load_fitpsf(std::string file) throw (const char*) {
  this->fitpsf.reset( new FittedPsf( this->params ) );
  try {
    this->fitpsf->readFits(file);
  }
  PY_CATCHALL 
}



void CWL::write_shearcat(std::string file, bool flush_log) throw (const char*) {
  try {
    this->shearcat->writeFits(file);
    if (flush_log) {
      // to flush the log we must delete it
      this->shear_log.reset();
    }
  }
  PY_CATCHALL
}
void CWL::load_shearcat(std::string file) throw (const char*) {
  this->shearcat.reset( new ShearCatalog( this->params ) );
  try {
    this->shearcat->readFits(file);
  }
  PY_CATCHALL 
}





void CWL::find_stars(std::string outfile) throw (const char*) {

  std::stringstream err;

  if (!this->image->loaded()) {
    throw "image not loaded";
  } else if (!this->weight_image->loaded()) {
    throw "weight image not loaded";
  } else  if (this->cat->size() == 0) {
    throw "catalog not loaded";
  }

  // create a log, defaults to stdout.  This log also writes results to the
  // FITS header of outfile when it is destroyed
  this->stars_log.reset(
      new FindStarsLog(this->params, this->stars_logfile, outfile)
  );

  // initialize the star cat from cat
  this->starcat.reset( new StarCatalog(*this->cat.get(),this->params) );

  // Update the sizes to more robust values
  this->starcat->calculateSizes(
      *this->image.get(), 
      this->weight_image.get(), 
      *this->trans.get());

  try {

    this->starcat->findStars(*this->stars_log.get());

  } catch(StarFinderException& e) {
    this->write_starcat(outfile);
    err<<"Caught StarFinderException: "<<e.what();
    throw err.str().c_str();
  } catch (...) {
    // we need this catchall because we *must* write the catalog
    // and I don't want to write out everything that might be caught
    this->write_starcat(outfile);
    throw "Caught unknown error";
  } 

  this->write_starcat(outfile);

}


void CWL::measure_psf(
    std::string psf_file, 
    std::string fitpsf_file) throw (const char*) {

  if (!this->image->loaded()) {
    throw "image not loaded";
  } else if (!this->weight_image->loaded()) {
    throw "weight image not loaded";
  } else  if (this->starcat->size() == 0) {
    throw "star catalog not loaded";
  }

 

  this->psfcat.reset(
      new PsfCatalog(*this->starcat.get(), this->params)
  );

  this->psf_log.reset(
      new PsfLog(this->params,
                 this->psf_logfile,
                 psf_file)
  );

  // first estimate of the sigma to use for shapelets decomp
  try {
    double sigma_p = this->psfcat->estimateSigma(
        *this->image.get(),
        this->weight_image.get(),
        *this->trans.get());
        
    int n_good_psf = this->psfcat->measurePsf(
        *this->image.get(),
        this->weight_image.get(),
        *this->trans.get(),
        sigma_p,
        *this->psf_log.get());

    // write but don't flush the log yet
    this->write_psfcat(psf_file, false);

    if (n_good_psf == 0) {
        throw "No successful PSF measurements";
    }

    this->fitpsf.reset(
        new FittedPsf(*this->psfcat.get(),
                      this->params,
                      *this->psf_log.get())
    );

    this->write_fitpsf(fitpsf_file);

    // fitpsf can change the log, so now over-write and flush the log
    this->write_psfcat(psf_file);

  } catch (...) {
    // we need this catchall because we *must* write the catalog
    // and I don't want to write out everything that might be caught
    //this->write_starcat(outfile);
    throw "Caught unknown error";
  } 


}

void CWL::measure_shear(std::string shear_file) throw (const char*) {
  if (!this->image->loaded()) {
    throw "image not loaded";
  } else if (!this->weight_image->loaded()) {
    throw "weight image not loaded";
  } else if (this->cat->size() == 0) {
    throw "catalog not loaded";
  } else if (this->starcat->size() == 0) {
    throw "star catalog not loaded";
  } else if (this->psfcat->size() == 0) {
    throw "psf catalog not loaded";
  } else if (this->fitpsf.get() == NULL) {
    throw "fitpsf not loaded";
  }


  this->shearcat.reset(
      new ShearCatalog(
        *this->cat.get(),
        *this->trans.get(),
        *this->fitpsf.get(),
        this->params)
  );

  this->shear_log.reset(
      new ShearLog(this->params,
                 this->shear_logfile,
                 shear_file)
  );

  this->shearcat->flagStars(*this->starcat.get());
 
  int n_shear = this->shearcat->measureShears(
      *this->image.get(),
      this->weight_image.get(),
      *this->shear_log);

  this->write_shearcat(shear_file);
}

/*
void CWL::check_load(const char types, size_t size) {

  for (size_t i=0; i<size; i++) {

  }
}
*/
