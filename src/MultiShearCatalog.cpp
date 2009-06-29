#include "CoaddCatalog.h"
#include "Ellipse.h"
#include "ConfigFile.h"
#include "dbg.h"
#include "Params.h"
#include "BVec.h"
#include "Ellipse.h"
#include "TMV.h"
#include "Pixel.h"
#include "Image.h"
#include "FittedPSF.h"
#include "Log.h"
#include "TimeVars.h"
#include "MultiShearCatalog.h"
#include "Form.h"

#define UseInverseTransform

//#define OnlyNImages 10

MultiShearCatalog::MultiShearCatalog(
    const CoaddCatalog& coaddcat, const ConfigFile& _params) :
  id(coaddcat.id), skypos(coaddcat.skypos), sky(coaddcat.sky),
  noise(coaddcat.noise), flags(coaddcat.flags), params(_params)
{
  Resize(coaddcat.size());

  // Fix flags to only have INPUT_FLAG set.
  // (I don't think this is necessary, but just in case.)
  for (size_t i=0; i<size(); i++) {
    flags[i] &= INPUT_FLAG;
  }

  // TODO: A lot of the parameter names that are coadd_* probably 
  // should be changed to multishear_*.  e.g. srclist

  // Read the names of the component image and catalog files from
  // the srclist file (given as params.coadd_srclist)
  ReadFileLists();

  // Loop over the files and read pixel lists for each object.
  // The Transformation and FittedPSF constructors get the name
  // information from the parameter file, so we use that to set the 
  // names of each component image here.
  for (size_t fnum=0; fnum<image_file_list.size(); fnum++) {
    std::string image_file = image_file_list[fnum];
    std::string psf_file = fitpsf_file_list[fnum];

    dbg<<"Reading image file: "<<image_file<<"\n";
    params["image_file"] = image_file;
    params["weight_file"] = image_file;
    params["dist_file"] = image_file;
    params["dist_hdu"] = 2;
    params["dist_method"] = "WCS";

    params["fitpsf_file"] = psf_file;

    GetImagePixelLists();
    dbg<<"\n";

#ifdef OnlyNImages
    if (fnum + 1 >= OnlyNImages) break;
#endif
  }
}

MultiShearCatalog::MultiShearCatalog(const ConfigFile& _params) :
  params(_params)
{
  Read();
}

MultiShearCatalog::~MultiShearCatalog()
{
  for(size_t k=0;k<trans.size();++k) {
    if (trans[k]) delete(trans[k]);
  }
  for(size_t k=0;k<fitpsf.size();++k) {
    if (fitpsf[k]) delete(fitpsf[k]);
  }
}

// ReadFileLists reads the srclist file specified in params
// and reads the names of the images and fitpsf 
void MultiShearCatalog::ReadFileLists()
{
  std::string file = params.get("coadd_srclist");
  if (!FileExists(file))
  {
    throw FileNotFound(file);
  }

  try 
  {
    dbg<<"Opening coadd srclist\n";
    std::ifstream flist(file.c_str(), std::ios::in);
    if (!flist)
    {
      throw ReadError("Unable to open source list file " + file);
    }

    image_file_list.clear();
    fitpsf_file_list.clear();

    std::string image_filename;
    std::string psf_filename;
    while (flist >> image_filename >> psf_filename) 
    {
      image_file_list.push_back(image_filename);
      fitpsf_file_list.push_back(psf_filename);
    }
  }
  catch (std::runtime_error& e)
  {
    throw ReadError("Error reading from "+file+" -- caught error\n" +
	e.what());
  }
  catch (...)
  {
    throw ReadError("Error reading from "+file+" -- caught unknown error");
  }

  Assert(image_file_list.size() == fitpsf_file_list.size());
}

void MultiShearCatalog::Resize(int n)
{
  pixlist.clear();
  pixlist.resize(n);

  image_indexlist.clear();
  image_indexlist.resize(n);

  image_cenlist.clear();
  image_cenlist.resize(n);

  input_flags.resize(n,0);
  nimages_found.resize(n, 0);
  nimages_gotpix.resize(n, 0);

  shear.resize(n);
  nu.resize(n);
  cov.resize(n);

  int gal_order = params.read<int>("shear_gal_order");
  BVec shape_default(gal_order,1.);
  shape_default.SetAllTo(DEFVALNEG);
  shape.resize(n,shape_default);
}

// Get pixel lists from the file specified in params
void MultiShearCatalog::GetImagePixelLists()
{
  std::auto_ptr<Image<double> > weight_im;
  Image<double> im(params,weight_im);

  int maxi=im.GetMaxI();
  int maxj=im.GetMaxJ();
  xdbg<<"MaxI: "<<maxi<<" MaxJ: "<<maxj<<"\n";

  // read transformation between ra/dec and x/y
  trans.push_back(new Transformation(params));
  const Transformation& trans1 = *trans.back();

  // read the psf
  fitpsf.push_back(new FittedPSF(params));
  const FittedPSF& fitpsf1 = *fitpsf.back();

  // which image index is this one?
  int image_index = fitpsf.size()-1;

#ifdef UseInverseTransform
  Transformation invtrans;
  Bounds invb = invtrans.MakeInverseOf(trans1,im.GetBounds(),4);
#endif

  // we are using the weight image so the noise and gain are 
  // dummy variables
  double noise = 0.0;
  double gain=0.0;
  // We always use the maximum aperture size here, since we don't know
  // how big the galaxy is yet, so we don't know what galap will be.
  double max_aperture = params.read<double>("shear_max_aperture");

  dbg<<"Extracting pixel lists\n";
  // loop over the the objects, if the object falls on the image get
  // the pixel list
  for (size_t i=0; i<skypos.size(); ++i) {

    // convert ra/dec to x,y in this image

#ifdef UseInverseTransform
    // Figure out a good starting point for the nonlinear solver:
    Position pxy;
    xdbg<<"skypos = "<<skypos[i]<<std::endl;
    if (!invb.Includes(skypos[i])) {
      xdbg<<"skypos "<<skypos[i]<<" not in "<<invb<<std::endl;
      continue;
    }
    invtrans.Transform(skypos[i],pxy);
    xdbg<<"invtrans(skypos) = "<<pxy<<std::endl;
#else
    Position pxy((double)maxi/2.,(double)maxj/2.);
#endif

    if (!trans1.InverseTransform(skypos[i], pxy) ) {
      std::stringstream err;
      err << "InverseTransform failed for position "<<skypos[i]<<".";
      dbg << "InverseTransform failed for position "<<skypos[i]<<".\n";
      dbg << "Initial guess was "<<pxy<<".\n";
      throw TransformationError(err.str());
    }

    double x=pxy.GetX();
    double y=pxy.GetY();
    if ( (x >= 0) && (x <= maxi) && (y >= 0) && (y <= maxj) ) {
      xdbg<<"("<<skypos[i]<<")  x="<<x<<" y="<<y<<"\n";

      nimages_found[i]++;

      // We don't actually need the psf here.  But we want to check
      // to make sure fitpsf(pxy) doesn't throw an exception:
      BVec psf1(fitpsf1.GetPSFOrder(), fitpsf1.GetSigma());
      try {
	psf1 = fitpsf1(pxy);
      } catch (Range_error& e) {
	xdbg<<"fittedpsf range error: \n";
	xdbg<<"p = "<<pxy<<", b = "<<e.b<<std::endl;
	input_flags[i] |= FITTEDPSF_EXCEPTION;
	continue;
      }

      // Make sure the use of trans in GetPixList won't throw:
      try {
	// We don't need to save skypos.  We just want to catch the range
	// error here, so we don't need to worry about it for dudx, etc.
	Position skypos1;
	trans1.Transform(pxy,skypos1);
      } catch (Range_error& e) {
	dbg<<"distortion range error: \n";
	xdbg<<"p = "<<pxy<<", b = "<<e.b<<std::endl;
	input_flags[i] |= TRANSFORM_EXCEPTION;
	continue;
      }

      long flag = 0;
      std::vector<Pixel> tmp_pixlist;
      GetPixList(
	  im,tmp_pixlist,pxy,
	  sky[i],noise,gain,weight_im.get(),trans1,max_aperture,flag);
      xdbg<<"Got pixellist, flag = "<<flag<<std::endl;

      // make sure not (edge or < 10 pixels) although edge is already
      // checked above
      if (flag == 0) {
	dbg<<"i = "<<i<<", pixlist.size = "<<pixlist.size()<<std::endl;
	Assert(i < pixlist.size());
	Assert(i < image_indexlist.size());
	Assert(i < image_cenlist.size());
	Assert(i < nimages_gotpix.size());
	pixlist[i].push_back(tmp_pixlist);
	image_indexlist[i].push_back(image_index);
	image_cenlist[i].push_back(pxy);
	nimages_gotpix[i]++;
      } else {
	input_flags[i] |= flag;
      }
      xdbg<<"added pixlist to main list\n";
    } else {
      xdbg<<"x,y not in valid bounds\n";
    }
  } // loop over objects
  dbg<<"Done extracting pixel lists\n";
}

void MeasureMultiShear1(
    const Position& cen, 
    const std::vector<std::vector<Pixel> > allpix,
    const std::vector<BVec>& psf,
    double gal_aperture, double max_aperture,
    int gal_order, int gal_order2,
    double f_psf, double min_gal_size, bool desqa,
    OverallFitTimes* times, ShearLog& log,
    std::complex<double>& shear, 
    tmv::SmallMatrix<double,2,2>& shearcov, BVec& shapelet,
    long& flag)
{
  //TODO: This function is too long.  It should really be broken 
  //      up in to smaller units...
  //      In which case most of the overlap between this and 
  //      MeasureSingleShear1 would not have to be duplicated.

  dbg<<"Start MeasureMultiShear1:\n";
  dbg<<"cen = "<<cen<<std::endl;
  dbg<<"allpix.size = "<<allpix.size()<<std::endl;
  for(size_t i=0;i<allpix.size();++i)
    dbg<<"allpix["<<i<<"].size = "<<allpix[i].size()<<std::endl;
  dbg<<"psf.size = "<<psf.size()<<std::endl;

  // Find harmonic mean of psf sizes:
  // MJ: Is it correct to use the harmonic mean of sigma^2?
  double sigma_p = 0.;
  Assert(psf.size() > 0);
  for(size_t i=0;i<psf.size();++i) sigma_p += 1./psf[i].GetSigma();
  sigma_p = double(psf.size()) / sigma_p;

  // We start with a native fit using sigma = 2*sigma_p:
  double sigma_obs = 2.*sigma_p;
  double galap = gal_aperture * sigma_obs;
  xdbg<<"galap = "<<gal_aperture<<" * "<<sigma_obs<<" = "<<galap<<std::endl;
  if (max_aperture > 0. && galap > max_aperture) {
    galap = max_aperture;
    xdbg<<"      => "<<galap<<std::endl;
  }
  xdbg<<"sigma_obs = "<<sigma_obs<<", sigma_p = "<<sigma_p<<std::endl;

  std::vector<std::vector<Pixel> > pix(allpix.size());
  int npix = 0;
  for(size_t i=0;i<pix.size();++i) {
    GetSubPixList(pix[i],allpix[i],galap,flag);
    npix += pix[i].size();
  }
  xdbg<<"npix = "<<npix<<std::endl;

  Ellipse ell;
  ell.FixGam();
  ell.CrudeMeasure(pix,sigma_obs);
  xdbg<<"Crude Measure: centroid = "<<ell.GetCen()<<", mu = "<<ell.GetMu()<<std::endl;
  //double mu_1 = ell.GetMu();

  int go = 2;
  int gal_size = (go+1)*(go+2)/2;

  if (times) ell.DoTimings();
  long flag1=0;
  if (ell.Measure(pix,go,sigma_obs,true,flag1,desqa)) {
    if (times) {
      times->ns_native++;
      times->ts_native_integ += ell.t_integ;
      times->ts_native_centroid += ell.t_centroid;
      times->ts_native_fixflux += ell.t_fixflux;
      times->ts_native_final += ell.t_final;
    }
    log.ns_native++;
    dbg<<"Successful native fit:\n";
    dbg<<"Z = "<<ell.GetCen()<<std::endl;
    dbg<<"Mu = "<<ell.GetMu()<<std::endl;
  }
  else {
    if (times) {
      times->nf_native++;
      times->tf_native_integ += ell.t_integ;
      times->tf_native_centroid += ell.t_centroid;
      times->tf_native_fixflux += ell.t_fixflux;
      times->tf_native_final += ell.t_final;
    }
    log.nf_native++;
    dbg<<"Native measurement failed\n";
    flag |= NATIVE_FAILED;
    return;
  }
  //double mu_2 = ell.GetMu();

  // Now we can do a deconvolving fit, but one that does not 
  // shear the coordinates.
  sigma_obs *= exp(ell.GetMu());
  if (sigma_obs < min_gal_size*sigma_p) {
    dbg<<"skip: galaxy is too small -- "<<sigma_obs<<" psf size = "<<sigma_p<<std::endl;
    if (times) times->nf_small++;
    log.nf_small++;
    flag |= TOO_SMALL;
    return;
  }
  galap *= exp(ell.GetMu());
  if (max_aperture > 0. && galap > max_aperture) galap = max_aperture;
  ell.SetMu(0);

  // Check - mu should now be zero:
  // This one works
  //ell.Measure(pix,go,sigma_obs,true,flag1,desqa);
  //xdbg<<"After native fit #2:\n";
  //xdbg<<"Mu = "<<ell.GetMu()<<std::endl;

  npix = 0;
  for(size_t i=0;i<pix.size();++i) {
    GetSubPixList(pix[i],allpix[i],galap,flag);
    npix += pix[i].size();
  }

  xdbg<<"npix = "<<npix<<std::endl;
  double sigpsq = pow(sigma_p,2);
  double sigma = pow(sigma_obs,2) - sigpsq;
  if (sigma < 0.1*sigpsq) sigma = 0.1*sigpsq;
  sigma = sqrt(sigma);
  xdbg<<"sigma = "<<sigma<<std::endl;
  ell.SetFP(f_psf);
  if (times) ell.ResetTimes();
  flag1 = 0;
  if (ell.Measure(pix,psf,go,sigma,false,flag1,desqa)) {
    if (times) {
      times->ns_mu++;
      times->ts_mu_integ += ell.t_integ;
      times->ts_mu_centroid += ell.t_centroid;
      times->ts_mu_fixflux += ell.t_fixflux;
      times->ts_mu_final += ell.t_final;
    }
    log.ns_mu++;
    if (flag1) {
      Assert((flag1 & ~(SHEAR_LOCAL_MIN | SHEAR_POOR_FIT)) == false);
      // This is just bumps the Ellipse flags up two spots
      // from the SHEAR_* to SHAPE_* flags.
      flag1 <<= 2;
      Assert((flag1 & ~(SHAPE_LOCAL_MIN | SHAPE_POOR_FIT)) == false);
      flag |= flag1;
    }
    dbg<<"Successful deconvolving fit:\n";
    xdbg<<"Mu = "<<ell.GetMu()<<std::endl;
  }
  else {
    if (times) {
      times->nf_mu++;
      times->tf_mu_integ += ell.t_integ;
      times->tf_mu_centroid += ell.t_centroid;
      times->tf_mu_fixflux += ell.t_fixflux;
      times->tf_mu_final += ell.t_final;
    }
    log.nf_mu++;
    dbg<<"Deconvolving measurement failed\n";
    flag |= DECONV_FAILED;
    ell.MeasureShapelet(pix,psf,shapelet);
    return;
  }
  //double mu_3 = ell.GetMu();

#if 0
  // This doesn't work.  So for now, keep mu, sigma the same
  // I'd rather have mu = 0 and the right sigma.
  // Maybe I need to change to fit to solve for sigma directly
  // rather than solve for mu.
  sigma *= exp(ell.GetMu());
  ell.SetMu(0);

  // Check - mu should now be zero:
  b_gal.SetSigma(sigma);
  dbg<<"Meausre with sigma = "<<sigma<<std::endl;
  ell.Measure(pix,psf,go,sigma,false,flag1,desqa,&b_gal);
  dbg<<"After deconvolving fit #2:\n";
  dbg<<"Mu = "<<ell.GetMu()<<std::endl;
  dbg<<"b_gal = "<<b_gal<<std::endl;
#endif

  // Measure the galaxy shape at the full order
  go = gal_order;
  gal_size = (go+1)*(go+2)/2;
  if (npix <= gal_size) {
    while (npix <= gal_size) { gal_size -= go+1; --go; }
    dbg<<"Too few pixels ("<<npix<<") for given gal_order. \n";
    dbg<<"Reduced gal_order to "<<go<<" gal_size = "<<gal_size<<std::endl;
    flag |= SHAPE_REDUCED_ORDER;
    if (go < 2) { flag |= TOO_SMALL; return; }
  }
  //shapelet = new BVec(go,sigma);
  shapelet.SetSigma(sigma);
  std::complex<double> gale = 0.;
  ell.MeasureShapelet(pix,psf,shapelet);
  xdbg<<"Measured b_gal = "<<shapelet<<std::endl;

  // Under normal circumstances, b20/b00 ~= conj(gamma)/sqrt(2)
  gale = std::complex<double>(shapelet[3],shapelet[4]);
  gale /= shapelet[0];
  if (std::abs(gale) < 0.5) {
    ell.SetGamma(conj(gale) * sqrt(2.));
    shear = ell.GetGamma();
  }

  // Finally, we calculate the shear in the deconvolved galaxy.
  //ell.FixMu();
  ell.UnFixGam();
  go = gal_order2;
  gal_size = (go+1)*(go+2)/2;
  if (npix <= gal_size) {
    while (npix <= gal_size) { gal_size -= go+1; --go; }
    dbg<<"Too few pixels ("<<npix<<") for given gal_order. \n";
    dbg<<"Reduced gal_order to "<<go<<" gal_size = "<<gal_size<<std::endl;
  }
  if (times) ell.ResetTimes();
  tmv::Matrix<double> cov(5,5);
  if (ell.Measure(pix,psf,go,sigma,false,flag,desqa,&cov)) {
    if (times) {
      times->ns_gamma++;
      times->ts_gamma_integ += ell.t_integ;
      times->ts_gamma_centroid += ell.t_centroid;
      times->ts_gamma_fixflux += ell.t_fixflux;
      times->ts_gamma_final += ell.t_final;
    }
    log.ns_gamma++;
    dbg<<"Successful Gamma fit\n";
    xdbg<<"Measured gamma = "<<ell.GetGamma()<<std::endl;
    shear = ell.GetGamma();
    tmv::SmallMatrix<double,2,2,tmv::RowMajor> cov1 = 
      cov.SubMatrix(2,4,2,4);
    if (!(cov1.Det() > 0.)) flag |= SHEAR_BAD_COVAR;
    else shearcov = cov1;
  }
  else {
    if (times) {
      times->nf_gamma++;
      times->tf_gamma_integ += ell.t_integ;
      times->tf_gamma_centroid += ell.t_centroid;
      times->tf_gamma_fixflux += ell.t_fixflux;
      times->tf_gamma_final += ell.t_final;
    }
    log.nf_gamma++;
    dbg<<"Measurement failed\n";
    flag |= SHEAR_FAILED;
    shear = ell.GetGamma();
    tmv::SmallMatrix<double,2,2,tmv::RowMajor> cov1 = 
      cov.SubMatrix(2,4,2,4);
    if (!(cov1.Det() > 0.)) flag |= SHEAR_BAD_COVAR;
    else shearcov = cov1;
    return;
  }
  //dbg<<"Stats: Mu: "<<mu_1<<"  "<<mu_2<<"  "<<mu_3<<std::endl;
  //dbg<<"       sigma = "<<sigma<<std::endl;
  //dbg<<"       gamma = "<<shear<<std::endl;

  // Finally measure the variance of the shear
  // TODO: Write function to directly measure shear variance
  // And also the BJ02 sigma_0 value.
  // (I'm not convinced that the above covariance matrix is a good estiamte.)
}

void MeasureMultiShear(
    const Position& cen, 
    const std::vector<std::vector<Pixel> >& pix,
    const std::vector<int>& image_index,
    const std::vector<Position>& image_cen,
    const std::vector<const FittedPSF*>& fitpsf,
    double gal_aperture, double max_aperture,
    int gal_order, int gal_order2,
    double f_psf, double min_gal_size, bool desqa,
    OverallFitTimes* times, ShearLog& log,
    std::complex<double>& shear, 
    tmv::SmallMatrix<double,2,2>& shearcov, BVec& shapelet,
    long& flag)
{
  Assert(fitpsf.size() > 0);
  // Calculate the psf from the fitted-psf formula:
  int psforder = fitpsf[0]->GetPSFOrder();
  double psfsigma = fitpsf[0]->GetSigma();
  for(size_t k=1;k<fitpsf.size();++k) 
    // The psforder should be the same in all the FittedPSF's,
    // but the sigma's are allowed to vary.
    Assert(fitpsf[k]->GetPSFOrder() == psforder);
  std::vector<BVec> psf(pix.size(), BVec(psforder, psfsigma));
  try {
    for(size_t i=0;i<pix.size();++i) {
      int k = image_index[i];
      psf[i] = (*fitpsf[k])(image_cen[i]);
    }
  } catch (Range_error& e) {
    // This shouldn't happen.  We already checked that fitpsf(cen) and
    // trans(cen) don't throw above in GetImagePixelLists
    dbg<<"Unexpected range error in MeasureMultiShear: \n";
    dbg<<"p = "<<e.p<<", b = "<<e.b<<std::endl;
    flag |= FITTEDPSF_EXCEPTION;
    throw;
  }

  // Do the real meat of the calculation:
  dbg<<"measure single shear cen = "<<cen<<std::endl;
  try { 
    MeasureMultiShear1(
	// Input data:
	cen, pix, psf,
	// Parameters:
	gal_aperture, max_aperture, gal_order, gal_order2, 
	f_psf, min_gal_size, desqa,
	// Time stats if desired:
	times,
	// Log information
	log,
	// Ouput values:
	shear, shearcov, shapelet, flag);
  } catch (tmv::Error& e) { 
    dbg<<"TMV Error thrown in MeasureSingleShear\n";
    dbg<<e<<std::endl;
    log.nf_tmverror++;
    flag |= TMV_EXCEPTION;
  } catch (std::runtime_error) {
    throw;
  } catch (...) {
    dbg<<"unkown exception in MeasureSingleShear\n";
    log.nf_othererror++;
    flag |= UNKNOWN_EXCEPTION;
  } 

}

int MultiShearCatalog::MeasureMultiShears(ShearLog& log)
{
  int ngals = skypos.size();
  dbg<<"ngals = "<<ngals<<std::endl;

  // Read some needed parameters
  double gal_aperture = params.read<double>("shear_aperture");
  double max_aperture = params.read("shear_max_aperture",0.);
  int gal_order = params.read<int>("shear_gal_order");
  int gal_order2 = params.read("shear_gal_order2",gal_order);
  double f_psf = params.read<double>("shear_f_psf");
  double min_gal_size = params.read<double>("shear_min_gal_size");
  bool desqa = params.read("des_qa",false);
  bool output_dots = params.read("output_dots",false);
  bool timing = params.read("timing",false);

  OverallFitTimes alltimes;

#ifdef ENDAT
  ngals = ENDAT;
#endif
  
  log.ngals = ngals;
#ifdef STARTAT
  log.ngals -= STARTAT;
#endif
#ifdef SINGLEGAL
  log.ngals = 1;
#endif
  log.ngoodin = std::count(flags.begin(),flags.end(),0);
  dbg<<log.ngoodin<<"/"<<log.ngals<<" galaxies with no input flags\n";

  // Main loop to measure shears
#ifdef _OPENMP
#pragma omp parallel 
  {
    try {
#endif
      OverallFitTimes times; // just for this thread
      ShearLog log1; // just for this thread
      log1.NoWriteLog();
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
      for(int i=0;i<ngals;++i) if (!flags[i]) {
#ifdef STARTAT
	if (i < STARTAT) continue;
#endif
#ifdef SINGLEGAL
	if (i < SINGLEGAL) continue;
	if (i > SINGLEGAL) break;
#endif

#ifdef _OPENMP
#pragma omp critical (output)
#endif
	{
	  if (output_dots) { std::cerr<<"."; std::cerr.flush(); }
	  dbg<<"galaxy "<<i<<":\n";
	  dbg<<"pos[i] = "<<skypos[i]<<std::endl;
	}

	if (pixlist[i].size() == 0) 
	{
	  dbg<<"no valid single epoch images.\n";
	  flags[i] = NO_SINGLE_EPOCH_IMAGES;
	  continue;
	}

	// Start with an error code of unknown failure, in case
	// something happens that I don't detect as an error.
	flags[i] = UNKNOWN_FAILURE;
	long flag1 = 0;

	MeasureMultiShear(
	    // Input data:
	    skypos[i], pixlist[i], image_indexlist[i], image_cenlist[i],
	    // Fitted PSF
	    fitpsf,
	    // Parameters:
	    gal_aperture, max_aperture, gal_order, gal_order2, 
	    f_psf, min_gal_size, desqa,
	    // Time stats if desired:
	    timing ? &times : 0, 
	    // Log information
	    log1,
	    // Ouput values:
	    shear[i], cov[i], shape[i], flag1);

	flags[i] = flag1;

	if (!flag1) {
	  ompdbg<<"Successful shear measurement: "<<shear[i]<<std::endl;
	}
	else {
	  ompdbg<<"Unsuccessful shear measurement\n"; 
	}

	if (timing) {
	  ompdbg<<"So far: ns = "<<times.ns_gamma<<",  nf = "<<times.nf_native;
	  ompdbg<<", "<<times.nf_mu<<", "<<times.nf_gamma<<std::endl;
	}

      }
#ifdef _OPENMP
#pragma omp critical (add_log)
#endif
      { 
	if (timing) alltimes += times;
	log += log1;
      }
#ifdef _OPENMP
    } 
    catch (...)
    {
      // This isn't supposed to happen.
      std::cerr<<"STATUS5BEG Caught some error in parallel region STATUS5END\n";
      exit(1);
    }
  }
#endif

  dbg<<log.ns_gamma<<" successful shear measurements, ";
  dbg<<ngals-log.ns_gamma<<" unsuccessful\n";
  log.ngood = std::count(flags.begin(),flags.end(),0);
  dbg<<log.ngood<<" with no flags\n";

  if (output_dots) {
    std::cerr
      <<std::endl
      <<"Success rate: "<<log.ns_gamma<<"/"<<log.ngoodin
      <<"  # with no flags: "<<log.ngood
      <<std::endl;
  }

  if (timing) {
    dbg<<"From timing structure:\n";
    dbg<<alltimes.ns_gamma<<" successful shear measurements, ";
    dbg<<alltimes.nf_native<<" + "<<alltimes.nf_mu;
    dbg<<" + "<<alltimes.nf_gamma<<" unsuccessful\n";
    std::cerr<<alltimes<<std::endl;
  }
  xdbg<<log<<std::endl;

  return log.ns_gamma;
}

void MultiShearCatalog::Write() const
{
  std::vector<std::string> files = MultiName(params, "multishear");

  for(size_t i=0; i<files.size(); ++i) {
    const std::string& file = files[i];
    dbg<<"Writing multishear catalog to file: "<<file<<std::endl;

    bool fitsio = false;
    if (params.keyExists("multishear_io")) {
      std::vector<std::string> ios = params["multishear_io"];
      Assert(ios.size() == files.size());
      fitsio = (ios[i] == "FITS");
    }
    else if (file.find("fits") != std::string::npos)
      fitsio = true;

    try
    {
      if (fitsio) {
	WriteFits(file);
      } else {
	std::string delim = "  ";
	if (params.keyExists("multishear_delim")) {
	  std::vector<std::string> delims = params["multishear_delim"];
	  Assert(delims.size() == files.size());
	  delim = delims[i];
	}
	else if (file.find("csv") != std::string::npos) delim = ",";
	WriteAscii(file,delim);
      }
    }
    catch (std::runtime_error& e)
    {
      throw WriteError("Error writing to "+file+" -- caught error\n" +
	  e.what());
    }
    catch (...)
    {
      throw WriteError("Error writing to "+file+" -- caught unknown error");
    }
  }
  dbg<<"Done Write ShearCatalog\n";
}

void MultiShearCatalog::WriteFits(std::string file) const
{
  // ! means overwrite existing file
  CCfits::FITS fits("!"+file, CCfits::Write);

  const int nfields=20;
  std::vector<string> colnames(nfields);
  std::vector<string> colfmts(nfields);
  std::vector<string> colunits(nfields);

  colnames[0] = params.get("multishear_id_col");
  colnames[1] = params.get("multishear_x_col");
  colnames[2] = params.get("multishear_y_col");
  colnames[3] = params.get("multishear_sky_col");
  colnames[4] = params.get("multishear_noise_col");
  colnames[5] = params.get("multishear_flags_col");
  colnames[6] = params.get("multishear_ra_col");
  colnames[7] = params.get("multishear_dec_col");
  colnames[8] = params.get("multishear_shear1_col");
  colnames[9] = params.get("multishear_shear2_col");
  colnames[10] = params.get("multishear_nu_col");
  colnames[11] = params.get("multishear_cov00_col");
  colnames[12] = params.get("multishear_cov01_col");
  colnames[13] = params.get("multishear_cov11_col");
  colnames[14] = params.get("multishear_order_col");
  colnames[15] = params.get("multishear_sigma_col");
  colnames[16] = params.get("multishear_coeffs_col");

  colnames[17] = params.get("multishear_nimages_found_col");
  colnames[18] = params.get("multishear_nimages_gotpix_col");
  colnames[19] = params.get("multishear_input_flags_col");

  int ncoeff = shape[0].size();
  dbg<<"ncoeff = "<<ncoeff<<std::endl;
  std::stringstream coeff_form;
  coeff_form << ncoeff << "D";
  colfmts[16] = coeff_form.str(); // shapelet coeffs
 
  colfmts[0] = "1J"; // id
  colfmts[1] = "1D"; // x
  colfmts[2] = "1D"; // y
  colfmts[3] = "1D"; // sky
  colfmts[4] = "1D"; // noise
  colfmts[5] = "1J"; // flags
  colfmts[6] = "1D"; // ra
  colfmts[7] = "1D"; // dec
  colfmts[8] = "1D"; // shear1
  colfmts[9] = "1D"; // shear2
  colfmts[10] = "1D"; // nu
  colfmts[11] = "1D"; // cov00
  colfmts[12] = "1D"; // cov01
  colfmts[13] = "1D"; // cov11
  colfmts[14] = "1J"; // order
  colfmts[15] = "1D"; // sigma

  colfmts[17] = "1J"; // nimages_found
  colfmts[18] = "1J"; // nimages_gotpix
  colfmts[19] = "1J"; // input_flags


  dbg<<"Before Create table"<<std::endl;
  CCfits::Table* table;
  table = fits.addTable("coadd_shearcat",size(),colnames,colfmts,colunits);

  // Header keywords
  std::string tmvvers = tmv::TMV_Version();
  std::string wlvers = WlVersion();

  table->addKey("tmvvers", tmvvers, "version of TMV code");
  table->addKey("wlvers", wlvers, "version of weak lensing code");

  std::string str;
  double dbl;
  int intgr;
  //CCfitsWriteParKey(params, table, "version", str);
  CCfitsWriteParKey(params, table, "noise_method", str);
  CCfitsWriteParKey(params, table, "dist_method", str);

  CCfitsWriteParKey(params, table, "shear_aperture", dbl);
  CCfitsWriteParKey(params, table, "shear_max_aperture", dbl);
  CCfitsWriteParKey(params, table, "shear_gal_order", intgr);
  CCfitsWriteParKey(params, table, "shear_gal_order2", intgr);
  CCfitsWriteParKey(params, table, "shear_min_gal_size", dbl);
  CCfitsWriteParKey(params, table, "shear_f_psf", dbl);

  // data
  // make vector copies for writing
  std::vector<double> ra(size());
  std::vector<double> dec(size());
  std::vector<double> shear1(size());
  std::vector<double> shear2(size());
  std::vector<double> cov00(size());
  std::vector<double> cov01(size());

  for(size_t i=0;i<size();i++) { 
    ra[i] = skypos[i].GetX();
    dec[i] = skypos[i].GetY();
  }
  for(size_t i=0;i<size();i++) { 
    shear1[i] = real(shear[i]);
    shear2[i] = imag(shear[i]);
  }
  std::vector<double> cov11(size());
  for(size_t i=0;i<size();i++) { 
    cov00[i] = cov[i](0,0);
    cov01[i] = cov[i](0,1);
    cov11[i] = cov[i](1,1);
  }

  int startrow=1;

  table->column(colnames[0]).write(id,startrow);
  table->column(colnames[3]).write(sky,startrow);
  table->column(colnames[4]).write(noise,startrow);
  table->column(colnames[5]).write(flags,startrow);
  table->column(colnames[6]).write(ra,startrow);
  table->column(colnames[7]).write(dec,startrow);
  table->column(colnames[8]).write(shear1,startrow);
  table->column(colnames[9]).write(shear2,startrow);
  table->column(colnames[10]).write(nu,startrow);
  table->column(colnames[11]).write(cov00,startrow);
  table->column(colnames[12]).write(cov01,startrow);
  table->column(colnames[13]).write(cov11,startrow);

  for (size_t i=0; i<size(); i++) {
    size_t row = i+1;
    long b_order = shape[i].GetOrder();
    double b_sigma = shape[i].GetSigma();

    table->column(colnames[14]).write(&b_order,1,row);
    table->column(colnames[15]).write(&b_sigma,1,row);
    double* cptr = (double *) shape[i].cptr();
    table->column(colnames[16]).write(cptr, ncoeff, 1, row);

  }

  table->column(colnames[17]).write(nimages_found,startrow);
  table->column(colnames[18]).write(nimages_gotpix,startrow);
  table->column(colnames[19]).write(input_flags,startrow);
}

void MultiShearCatalog::WriteAscii(std::string file, std::string delim) const
{
  Assert(id.size() == size());
  Assert(skypos.size() == size());
  Assert(sky.size() == size());
  Assert(noise.size() == size());
  Assert(flags.size() == size());
  Assert(shear.size() == size());
  Assert(nu.size() == size());
  Assert(cov.size() == size());
  Assert(shape.size() == size());
  Assert(nimages_found.size() == size());
  Assert(nimages_gotpix.size() == size());
  Assert(input_flags.size() == size());

  std::ofstream fout(file.c_str());
  if (!fout) {
    throw WriteError("Error opening shear file"+file);
  }

  Form hexform; hexform.hex().trail(0);

  for(size_t i=0;i<size();i++) {
    fout
      << id[i] << delim
      << skypos[i].GetX() << delim
      << skypos[i].GetY() << delim
      << sky[i] << delim
      << noise[i] << delim
      << hexform(flags[i]) << delim
      << real(shear[i]) << delim
      << imag(shear[i]) << delim
      << nu[i] << delim
      << cov[i](0,0) << delim
      << cov[i](0,1) << delim
      << cov[i](1,1) << delim
      << shape[i].GetOrder() << delim
      << shape[i].GetSigma() << delim
      << nimages_found[i] << delim
      << nimages_gotpix[i] << delim
      << input_flags[i];
    for(size_t j=0;j<shape[i].size();++j)
      fout << delim << shape[i][j];
    fout << std::endl;
  }
}

void MultiShearCatalog::Read()
{
  std::string file = Name(params,"multishear",false,true);
  // false,true = input_prefix=false, mustexist=true.
  // It is an input here, but it is in the output_prefix directory.
  dbg<< "Reading Shear cat from file: " << file << std::endl;

  bool fitsio = false;
  if (params.keyExists("multishear_io"))
    fitsio = (params["multishear_io"] == "FITS");
  else if (file.find("fits") != std::string::npos)
    fitsio = true;

  if (!FileExists(file))
  {
    throw FileNotFound(file);
  }
  try
  {
    if (fitsio) {
      ReadFits(file);
    } else {
      std::string delim = "  ";
      if (params.keyExists("multishear_delim")) delim = 
	params["multishear_delim"];
      else if (file.find("csv") != std::string::npos) delim = ",";
      ReadAscii(file,delim);
    }
  }
  catch (std::runtime_error& e)
  {
    throw ReadError("Error reading from "+file+" -- caught error\n" +
	e.what());
  }
  catch (...)
  {
    throw ReadError("Error reading from "+file+" -- caught unknown error");
  }
  dbg<<"Done Read ShearCatalog\n";
}

void MultiShearCatalog::ReadFits(std::string file)
{
  int hdu = params.read("multishear_hdu",2);

  dbg<<"Opening FITS file at hdu "<<hdu<<std::endl;
  // true means read all as part of the construction
  CCfits::FITS fits(file, CCfits::Read, hdu-1, true);

  CCfits::ExtHDU& table=fits.extension(hdu-1);

  long nrows=table.rows();

  dbg<<"  nrows = "<<nrows<<std::endl;
  if (nrows <= 0) {
    throw ReadError("ShearCatalog found to have 0 rows.  Must have > 0 rows.");
  }

  std::string id_col=params.get("multishear_id_col");
  std::string ra_col=params.get("multishear_ra_col");
  std::string dec_col=params.get("multishear_dec_col");
  std::string sky_col=params.get("multishear_sky_col");
  std::string noise_col=params.get("multishear_noise_col");
  std::string flags_col=params.get("multishear_flags_col");
  std::string shear1_col=params.get("multishear_shear1_col");
  std::string shear2_col=params.get("multishear_shear2_col");
  std::string nu_col=params.get("multishear_nu_col");
  std::string cov00_col=params.get("multishear_cov00_col");
  std::string cov01_col=params.get("multishear_cov01_col");
  std::string cov11_col=params.get("multishear_cov11_col");
  std::string order_col=params.get("multishear_order_col");
  std::string sigma_col=params.get("multishear_sigma_col");
  std::string coeffs_col=params.get("multishear_coeffs_col");
  std::string nimages_found_col=params.get("multishear_nimages_found_col");
  std::string nimages_gotpix_col=params.get("multishear_nimages_gotpix_col");
  std::string input_flags_col=params.get("multishear_input_flags_col");

  long start=1;
  long end=nrows;

  dbg<<"Reading columns"<<std::endl;
  dbg<<"  "<<id_col<<std::endl;
  table.column(id_col).read(id, start, end);

  dbg<<"  "<<ra_col<<"  "<<dec_col<<std::endl;
  skypos.resize(nrows);
  std::vector<double> ra;
  std::vector<double> dec;
  table.column(ra_col).read(ra, start, end);
  table.column(dec_col).read(dec, start, end);
  for(long i=0;i<nrows;++i) skypos[i] = Position(ra[i],dec[i]);

  dbg<<"  "<<sky_col<<std::endl;
  table.column(sky_col).read(sky, start, end);

  dbg<<"  "<<noise_col<<std::endl;
  table.column(noise_col).read(noise, start, end);

  dbg<<"  "<<flags_col<<std::endl;
  table.column(flags_col).read(flags, start, end);

  dbg<<"  "<<shear1_col<<"  "<<shear2_col<<std::endl;
  shear.resize(nrows);
  std::vector<double> shear1;
  std::vector<double> shear2;
  table.column(shear1_col).read(shear1, start, end);
  table.column(shear2_col).read(shear2, start, end);
  for(long i=0;i<nrows;++i) {
    shear[i] = std::complex<double>(shear1[i],shear2[i]);
  }

  dbg<<"  "<<nu_col<<std::endl;
  table.column(nu_col).read(nu, start, end);

  dbg<<"  "<<cov00_col<<"  "<<cov01_col<<"  "<<cov11_col<<std::endl;
  cov.resize(nrows);
  std::vector<double> cov00;
  std::vector<double> cov01;
  std::vector<double> cov11;
  table.column(cov00_col).read(cov00, start, end);
  table.column(cov01_col).read(cov01, start, end);
  table.column(cov11_col).read(cov11, start, end);
  for(long i=0;i<nrows;++i) {
    cov[i] = tmv::ListInit, cov00[i], cov01[i], cov01[i], cov11[i];
  }

  dbg<<"  "<<sigma_col<<"  "<<order_col<<std::endl;
  // temporary
  std::vector<double> sigma;
  std::vector<int> order;
  table.column(sigma_col).read(sigma, start, end);
  table.column(order_col).read(order, start, end);

  shape.reserve(nrows);
  for (size_t i=0; i<size(); i++) {
    size_t row=i+1;

    shape.push_back(BVec(order[i],sigma[i]));
    int ncoeff=(order[i]+1)*(order[i]+2)/2;

    std::valarray<double> coeffs;
    table.column(coeffs_col).read(coeffs, row);

    double* ptri = (double* ) shape[i].cptr();
    for (int j=0; j<ncoeff; ++j) {
      ptri[i] = coeffs[i];
    }
  }

  dbg<<"  "<<nimages_found_col<<std::endl;
  table.column(nimages_found_col).read(nimages_found, start, end);

  dbg<<"  "<<nimages_gotpix_col<<std::endl;
  table.column(nimages_gotpix_col).read(nimages_gotpix, start, end);

  dbg<<"  "<<input_flags_col<<std::endl;
  table.column(input_flags_col).read(input_flags, start, end);
}

void MultiShearCatalog::ReadAscii(std::string file, std::string delim)
{
  std::ifstream fin(file.c_str());
  if (!fin) {
    throw ReadError("Error opening stars file"+file);
  }

  id.clear(); skypos.clear(); sky.clear(); noise.clear(); flags.clear();
  shear.clear(); nu.clear(); cov.clear(); shape.clear();
  nimages_found.clear(); nimages_gotpix.clear(); input_flags.clear();

  if (delim == "  ") {
    ConvertibleString flag;
    long id1,b_order,nfound,ngotpix,inflag;
    double ra,dec,sky1,noise1,s1,s2,nu1,c00,c01,c11,b_sigma;
    while ( fin >> id1 >> ra >> dec >> sky1 >> noise1 >>
	flag >> s1 >> s2 >> nu1 >>
	c00 >> c01 >> c11 >> b_order >> b_sigma >>
	nfound >> ngotpix >> inflag )
    {
      id.push_back(id1);
      skypos.push_back(Position(ra,dec));
      sky.push_back(sky1);
      noise.push_back(noise1);
      flags.push_back(flag);
      shear.push_back(std::complex<double>(s1,s2));
      nu.push_back(nu1);
      cov.push_back(tmv::SmallMatrix<double,2,2>());
      cov.back() = tmv::ListInit, c00, c01, c01, c11;
      shape.push_back(BVec(b_order,b_sigma));
      for(size_t j=0;j<shape.back().size();++j)
	fin >> shape.back()[j];
      nimages_found.push_back(nfound);
      nimages_gotpix.push_back(ngotpix);
      input_flags.push_back(inflag);
    }
  } else {
    if (delim.size() > 1) {
      // getline only works with single character delimiters.
      // Since I don't really expect a multicharacter delimiter to
      // be used ever, I'm just going to throw an exception here 
      // if we do need it, and I can write the workaround then.
      throw ParameterError("ReadAscii delimiter must be a single character");
    }
    char d = delim[0];
    long b_order;
    double ra,dec,s1,s2,b_sigma,c00,c01,c11;
    ConvertibleString temp;
    while (getline(fin,temp,d))
    {
      id.push_back(temp);
      getline(fin,temp,d); ra = temp;
      getline(fin,temp,d); dec = temp;
      skypos.push_back(Position(ra,dec));
      getline(fin,temp,d); sky.push_back(temp);
      getline(fin,temp,d); noise.push_back(temp);
      getline(fin,temp,d); flags.push_back(temp);
      getline(fin,temp,d); s1 = temp;
      getline(fin,temp,d); s2 = temp;
      shear.push_back(std::complex<double>(s1,s2));
      getline(fin,temp,d); nu.push_back(temp);
      getline(fin,temp,d); c00 = temp;
      getline(fin,temp,d); c01 = temp;
      getline(fin,temp,d); c11 = temp;
      cov.push_back(tmv::SmallMatrix<double,2,2>());
      cov.back() = tmv::ListInit, c00, c01, c01, c11;
      getline(fin,temp,d); b_order = temp;
      getline(fin,temp,d); b_sigma = temp;
      shape.push_back(BVec(b_order,b_sigma));
      for(size_t j=0;j<shape.back().size()-1;++j) {
	getline(fin,temp,d); shape.back()[j] = temp;
      }
      getline(fin,temp); shape.back()[shape.back().size()-1] = temp;
      nimages_found.push_back(temp);
      nimages_gotpix.push_back(temp);
      input_flags.push_back(temp);
    }
  }
}
