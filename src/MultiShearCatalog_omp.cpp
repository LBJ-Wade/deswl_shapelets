
#include "ConfigFile.h"
#include "dbg.h"
#include "Params.h"
#include "Log.h"
#include "TimeVars.h"
#include "MultiShearCatalog.h"
#include "Ellipse.h"
#include "ShearCatalog.h"
#include "ShearCatalogTree.h"

int MultiShearCatalog::MeasureMultiShears(const Bounds& b, ShearLog& log)
{
  dbg<<"Start MeasureMultiShears for b = "<<b<<std::endl;
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

  int nsuccess = 0;

#ifdef ENDAT
  ngals = ENDAT;
#endif
  
  // Main loop to measure shears
#ifdef _OPENMP
#pragma omp parallel 
  {
    try {
#endif
      OverallFitTimes times; // just for this thread
      ShearLog log1(params); // just for this thread
      log1.NoWriteLog();
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
      for(int i=0;i<ngals;++i) 
      {
	if (!b.Includes(skypos[i])) 
	{
	  xxdbg<<"Skipping galaxy "<<i<<" because "<<skypos[i]<<"not in bounds\n";
	  continue;
	}
	if (flags[i]) 
	{
	  xxdbg<<"Skipping galaxy "<<i<<" because flag = "<<flags[i]<<std::endl;
	  continue;
	}
#ifdef STARTAT
	if (i < STARTAT) continue;
#endif
#ifdef SINGLEGAL
	if (i < SINGLEGAL) continue;
	if (i > SINGLEGAL) break;
#endif

	if (output_dots)
#ifdef _OPENMP
#pragma omp critical (output)
#endif
	{
	  std::cerr<<"."; std::cerr.flush(); 
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
	    skypos[i], pixlist[i], psflist[i],
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
	  dbg<<"Successful shear measurement: "<<shear[i]<<std::endl;
	  nsuccess++;
	}
	else {
	  dbg<<"Unsuccessful shear measurement\n"; 
	}

	if (timing) {
	  dbg<<"So far: ns = "<<times.ns_gamma<<",  nf = "<<times.nf_native;
	  dbg<<", "<<times.nf_mu<<", "<<times.nf_gamma<<std::endl;
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

  dbg<<nsuccess<<" successful shear measurements in this pass.\n";
  dbg<<log.ns_gamma<<" successful shear measurements so far.\n";

  if (timing) {
    dbg<<"From timing structure:\n";
    dbg<<alltimes.ns_gamma<<" successful shear measurements, ";
    dbg<<alltimes.nf_native<<" + "<<alltimes.nf_mu;
    dbg<<" + "<<alltimes.nf_gamma<<" unsuccessful\n";
    std::cerr<<alltimes<<std::endl;
  }
  xdbg<<log<<std::endl;

  return nsuccess;
}

static void GetImagePixList(
    std::vector<PixelList>& pixlist,
    std::vector<BVec>& psflist,
    std::vector<std::complex<double> >& se_shearlist,
    std::vector<double>& se_sizelist,
    long& input_flags, int& nimages_found, int& nimages_gotpix,
    const Position& skypos,
    const Image<double>& im,
    const Transformation& trans,
    const Transformation& invtrans,
    BVec& psf, const FittedPSF& fitpsf,
    const ShearCatalog& shearcat,
    const ShearCatalogTree& shearcat_tree,
    const Image<double>*const weight_im,
    const double noise, const double gain,
    const std::string& sky_method, const double mean_sky, 
    const Image<float>*const sky_map,
    const double gal_aperture, const double max_aperture)
{
  Assert(psflist.size() == pixlist.size());
  Assert(se_shearlist.size() == pixlist.size());
  Assert(se_sizelist.size() == pixlist.size());

  // Convert ra/dec to x,y in this image

  // First, figure out a good starting point for the nonlinear solver:
  Position pos;
  dbg<<"skypos = "<<skypos<<std::endl;
  invtrans.Transform(skypos,pos);
  xdbg<<"invtrans(skypos) = "<<pos<<std::endl;

  // Now do the full non-linear solver, which should be pretty fast
  // given the decent initial guess.
  if (!trans.InverseTransform(skypos, pos) ) {
    dbg << "InverseTransform failed for position "<<skypos<<".\n";
    dbg << "Initial guess was "<<pos<<".\n";
    input_flags |= TRANSFORM_EXCEPTION;
    //std::stringstream err;
    //err << "InverseTransform failed for position "<<skypos<<".";
    //throw TransformationError(err.str());
  }
  xdbg<<"after exact InverseTransform: pos -> "<<pos<<std::endl;
  if (!(fitpsf.GetBounds().Includes(pos))) 
  {
    xdbg<<"Reject pos "<<pos<<" not in fitpsf bounds ";
    xdbg<<fitpsf.GetBounds()<<std::endl;
    return;
  }

  nimages_found++;

  try {
    psf = fitpsf(pos);
  } catch (Range_error& e) {
    xdbg<<"fittedpsf range error: \n";
    xdbg<<"p = "<<pos<<", b = "<<e.b<<std::endl;
    input_flags |= FITTEDPSF_EXCEPTION;
    return;
  }

  // Make sure the use of trans in GetPixList won't throw:
  try {
    // We don't need to save skypos.  We just want to catch the range
    // error here, so we don't need to worry about it for dudx, etc.
    Position skypos1;
    trans.Transform(pos,skypos1);
  } catch (Range_error& e) {
    dbg<<"distortion range error: \n";
    xdbg<<"p = "<<pos<<", b = "<<e.b<<std::endl;
    input_flags |= TRANSFORM_EXCEPTION;
    return;
  }

  // Find the nearest object in the shear catalog:
  int nearest = shearcat_tree.FindNearestTo(pos);

  xxdbg<<"pixlist.size = "<<pixlist.size()<<std::endl;
  xxdbg<<"psflist.size = "<<psflist.size()<<std::endl;
  xxdbg<<"se_shearlist.size = "<<se_shearlist.size()<<std::endl;
  xxdbg<<"se_sizelist.size = "<<se_sizelist.size()<<std::endl;
  Assert(psflist.size() == pixlist.size());
  Assert(se_shearlist.size() == pixlist.size());
  Assert(se_sizelist.size() == pixlist.size());

  std::complex<double> se_shear = 0.;
  double se_size = 0.;
  double galap = max_aperture;
  if (std::abs(shearcat.pos[nearest] - pos) < 1.)
  { 
    se_shear = shearcat.shear[nearest];
    se_size = shearcat.shape[nearest].GetSigma();
    // If we have a good measurement from the single_epoch image, use
    // that sigma for the size.
    // But expand it by 30% in case we need it.
    galap = gal_aperture * se_size * 1.3;
    if (galap > max_aperture) galap = max_aperture;
  }

  // Calculate the local sky value.
  long flag = 0;
  double sky=0;
  if (sky_method == "MEAN")
    sky = mean_sky;
  else if (sky_method == "NEAREST")
    sky = shearcat.sky[nearest];
  else 
  {
    Assert(sky_method == "MAP");
    sky = GetLocalSky(*sky_map,pos,trans,galap,flag);
    // If no pixels in aperture, then
    // a) something is very wrong, but
    // b) use the NEAREST method instead.
    if (flag & BKG_NOPIX) sky = shearcat.sky[nearest];
    input_flags |= flag;
  }

  pixlist.push_back(PixelList());
#ifdef PIXELLIST_BLOCK
  pixlist.back().UseBlockMem();
#endif
  xdbg<<"pixlist.size = "<<pixlist.size()<<std::endl;

  flag = 0;
  GetPixList(
      im,pixlist.back(),pos,
      sky,noise,gain,weight_im,trans,galap,flag);
  xdbg<<"Got pixellist, flag = "<<flag<<std::endl;

  // Make sure not (edge or < 10 pixels) although edge is already
  // checked above
  if (flag == 0) {
    dbg<<"pixlist.size = "<<pixlist.size()<<std::endl;
    psflist.push_back(psf);

    // If object in the ShearCatalog is within 1 arcsec of the position
    // then assume it is the correct object.
    // Otherwise put in 0's to indicate we don't have an initial
    // guess for this object (on this image).
    if (std::abs(shearcat.pos[nearest] - pos) < 1.)
    { 
      se_shearlist.push_back(shearcat.shear[nearest]);
      se_sizelist.push_back(shearcat.shape[nearest].GetSigma());
    }
    else 
    {
      se_shearlist.push_back(0.);
      se_sizelist.push_back(0.);
    }
    xxdbg<<"se_shearlist.size = "<<se_shearlist.size()<<std::endl;
    xxdbg<<"se_sizelist.size = "<<se_sizelist.size()<<std::endl;
    nimages_gotpix++;
  } else {
    input_flags |= flag;
    pixlist.pop_back();
    xxdbg<<"pixlist.size = "<<pixlist.size()<<std::endl;
  }
}

// Get pixel lists from the file specified in params
void MultiShearCatalog::GetImagePixelLists(int se_index, const Bounds& b)
{
  dbg<<"Start GetImagePixelLists: se_index = "<<se_index<<std::endl;

  // If the skybounds for each shear catalog have been saved, then
  // we might be able to skip the ShearCatalog load.
  if (int(saved_se_skybounds.size()) > se_index)
  {
    Bounds se_skybounds = saved_se_skybounds[se_index];
    dbg<<"saved bounds for image "<<se_index<<" = "<<se_skybounds;
    if (!se_skybounds.Intersects(b)) 
    {
      dbg<<"Skipping index "<<se_index<<" because bounds don't intersect\n";
      return;
    }
  }

  // Read the shear catalog
  ShearCatalog shearcat(params);
  Bounds se_skybounds = shearcat.skybounds;
  Bounds se_bounds = shearcat.bounds;
  dbg<<"bounds for image "<<se_index<<" = "<<se_skybounds;

  // Skip this file if none of the objects in it are in this section of sky.
  if (int(saved_se_skybounds.size()) <= se_index) // Then save the se_skybounds
  {
    Assert(int(saved_se_skybounds.size()) == se_index);
    saved_se_skybounds.push_back(se_skybounds);
  }
  if (!se_skybounds.Intersects(b)) 
  {
    dbg<<"Skipping index "<<se_index<<" because bounds don't intersect\n";
    return;
  }

  // Read transformation between ra/dec and x/y
  Transformation trans(params);

  // Read the psf
  FittedPSF fitpsf(params);

  // Make a tree of the shear catalog to more easily find the nearest
  // single-epoch object to each coadd detection.
  ShearCatalogTree shearcat_tree(shearcat);

  // Figure out which method we are going to use to calculate the 
  // local sky values.
  std::string sky_method = params.get("multishear_sky_method");
  Assert(sky_method == "MEAN" || sky_method == "NEAREST"
      || sky_method == "MAP");
  double mean_sky=0.;
  std::auto_ptr<Image<float> > sky_map(0);
  if (sky_method == "MEAN")
  {
    for(size_t i=0;i<shearcat.size();++i) mean_sky += shearcat.sky[i];
    mean_sky /= shearcat.size();
  }
  if (sky_method == "MAP")
  {
    int skymap_hdu = params.read("skymap_hdu",1);
    sky_map.reset(new Image<float>(Name(params,"skymap",true),skymap_hdu));
  }

  BVec psf(fitpsf.GetPSFOrder(), fitpsf.GetSigma());

  // Make an inverse transformation that we will use as a starting 
  // point for the more accurate InverseTransform function.
  Transformation invtrans;
  Bounds invb = invtrans.MakeInverseOf(trans,se_bounds,4);
  dbg<<"skybounds = "<<skybounds<<std::endl;
  dbg<<"se_skybounds = "<<se_skybounds<<std::endl;
  dbg<<"se_bounds = "<<se_bounds<<std::endl;
  dbg<<"invb = "<<invb<<std::endl;

  // We always use the maximum aperture size here, since we don't know
  // how big the galaxy is yet, so we don't know what galap will be.
  double gal_aperture = params.get("shear_aperture");
  double max_aperture = params.get("shear_max_aperture");

  // Load the image
  // The bounds needed are 
  std::auto_ptr<Image<double> > im;
  std::auto_ptr<Image<double> > weight_im;
  if (skybounds.Includes(se_skybounds))
    im.reset(new Image<double>(params,weight_im));
  else
  {
    Bounds intersect = invb & skybounds;
    dbg<<"intersect = invb & skybounds = "<<intersect<<std::endl;

    Bounds subb;
    subb += invtrans(intersect.Get00());
    subb += invtrans(intersect.Get01());
    subb += invtrans(intersect.Get10());
    subb += invtrans(intersect.Get11());

    // Grow bounds by max_aperture
    tmv::SmallMatrix<double,2,2> D;
    trans.GetDistortion(subb.Center(),D);
    double det = std::abs(D.Det());
    double pixscale = sqrt(det); // arcsec/pixel
    subb.AddBorder(max_aperture / pixscale);

    dbg<<"subb = "<<subb<<std::endl;
    im.reset(new Image<double>(params,weight_im,subb));
  }

  // We are using the weight image so the noise and gain are dummy variables
  Assert(weight_im.get());
  double noise = 0.0;
  double gain = 0.0;

  dbg<<"Extracting pixel lists\n";
  // loop over the the objects, if the object falls on the image get
  // the pixel list
  Assert(pixlist.size() == skypos.size());
  Assert(psflist.size() == skypos.size());
  Assert(se_shearlist.size() == skypos.size());
  Assert(se_sizelist.size() == skypos.size());

  const int n_size = size();
#ifdef _OPENMP
//#pragma omp parallel for
#endif
  for (int i=0; i<n_size; ++i) 
  {
    Assert(i < int(flags.size()));
    Assert(i < int(skypos.size()));
    Assert(i < int(pixlist.size()));
    Assert(i < int(psflist.size()));
    Assert(i < int(se_shearlist.size()));
    Assert(i < int(se_sizelist.size()));
    Assert(i < int(input_flags.size()));
    Assert(i < int(nimages_found.size()));
    Assert(i < int(nimages_gotpix.size()));
    if (!flags[i] && 
	b.Includes(skypos[i]) &&
	invb.Includes(skypos[i]))
    {
      GetImagePixList(
	  pixlist[i], psflist[i], se_shearlist[i], se_sizelist[i],
	  input_flags[i], nimages_found[i], nimages_gotpix[i], skypos[i], 
	  *im, trans, invtrans, psf, fitpsf, shearcat, shearcat_tree,
	  weight_im.get(), noise, gain,
	  sky_method, mean_sky, sky_map.get(),
	  gal_aperture, max_aperture);
    }
  } 
  // Keep track of how much memory we are using.
  // TODO: introduce a parameter max_memory and check to make sure
  // we stay within the allowed memory usage.
  bool output_dots = params.read("output_dots",false);
  if (output_dots) 
  {
    std::cerr<<"Using image# "<<se_index;
    std::cerr<<"... Memory Usage in MultiShearCatalog = ";
    std::cerr<<CalcMemoryFootprint()<<" MB";
    //std::cerr<<" ["<<(MemoryFootprint(pixlist)/1024./1024.)<<"]";
    //std::cerr<<" ["<<(MemoryFootprint(shape)/1024./1024.)<<"]";
    std::cerr<<"\n";
    //std::cerr<<"Memory used in PixelList pool = "<<
      //pool_allocator<Pixel,PIXELLIST_BLOCK>::total_memory_used()/(1024.*1024.)<<
      //" MB\n";
  }

  dbg<<"Done extracting pixel lists\n";
}

