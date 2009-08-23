
#include <valarray>
#include "TMV.h"
#include <CCfits/CCfits>

#include "CoaddCatalog.h"
#include "Ellipse.h"
#include "ConfigFile.h"
#include "dbg.h"
#include "Params.h"
#include "BVec.h"
#include "Ellipse.h"
#include "Pixel.h"
#include "Image.h"
#include "FittedPSF.h"
#include "Log.h"
#include "TimeVars.h"
#include "MultiShearCatalog.h"
#include "ShearCatalog.h"
#include "Form.h"
#include "WriteParKey.h"

#define UseInverseTransform

//#define OnlyNImages 30

// TODO: Turn this into a ok_input_flag parameter and check against it
// with if ((flag & ~ok_input_flag)) { .. skip .. }
#define SKIP_INPUT_FLAGS

MultiShearCatalog::MultiShearCatalog(
    const CoaddCatalog& coaddcat, const ConfigFile& _params) :
  id(coaddcat.id), skypos(coaddcat.skypos), sky(coaddcat.sky),
  noise(coaddcat.noise), flags(coaddcat.flags),
  skybounds(coaddcat.skybounds), params(_params)
{
  Resize(coaddcat.size());

  for (size_t i=0; i<size(); i++) {
    if (flags[i]) flags[i] = INPUT_FLAG;
  }

  // Read the names of the component image and catalog files from
  // the srclist file (given as params.coadd_srclist)
  ReadFileLists();
}

std::vector<Bounds> MultiShearCatalog::SplitBounds(double side)
{
  double xrange = skybounds.GetXMax() - skybounds.GetXMin();
  double yrange = skybounds.GetYMax() - skybounds.GetYMin();

  int nx = int(floor(xrange/side))+1;
  int ny = int(floor(yrange/side))+1;
  return skybounds.Divide(nx,ny);
}

int MultiShearCatalog::GetPixels(const Bounds& b)
{
  // The pixlist object takes up a lot of memory, so at the start 
  // of this function, I clear it out, along with image_indexlist and
  // image_cenlist.  Then we loop over sections of the coaddcat and
  // only measure the shears for the objects within one section at a time.

  // MJ: pixlist used to be a vector<vector<vector<Pixel> > > object.
  // Then at the start of this function, I would just call
  // pixlist.clear();
  // This worked fine for serial code, but when I tried openmp, a memory
  // leak showed up.  Apparently, clear didn't actually release the 
  // memory in that case, so on the next pass through GetPixels,
  // the other pixel lists would be given new memory rather than reacquiring
  // the memory from the cleared vectors.
  // I don't really understand why this happened, especially since this 
  // function is called from outside the openmp parallel region, so I 
  // wouldn't have thought omp would be a factor at all.
  //
  // Update: this didn't work either.  Hmmm....  
  const int n = pixlist.size();
  //pixlist.reset(0);
  //pixlist.reset(new std::vector<std::vector<std::vector<Pixel> > >(n));

  for (int i=0;i<n;++i) 
  {
    pixlist[i].clear();
    image_indexlist[i].clear();
    image_cenlist[i].clear();
  }

  dbg<<"Start GetPixels for b = "<<b<<std::endl;
  // Loop over the files and read pixel lists for each object.
  // The Transformation and FittedPSF constructors get the name
  // information from the parameter file, so we use that to set the 
  // names of each component image here.
  for (size_t fnum=0; fnum<image_file_list.size(); fnum++) {

    dbg<<"bounds for image "<<fnum<<" = "<<image_skybounds[fnum];
    // Skip this file in none of the objects in it are in this section of sky.
    if (!image_skybounds[fnum].Intersects(b)) 
    {
      dbg<<"Skipping fnum = "<<fnum<<", because bounds don't intersect\n";
      continue;
    }

    // Get the file names
    std::string image_file = image_file_list[fnum];
    std::string shear_file = shear_file_list[fnum];
    std::string psf_file = fitpsf_file_list[fnum];

    dbg<<"Reading image file: "<<image_file<<"\n";
    // Set the appropriate parameters
    params["image_file"] = image_file;
    params["weight_file"] = image_file;
    params["dist_file"] = image_file;
    params["shear_file"] = shear_file;
    params["dist_hdu"] = 2;
    params["dist_method"] = "WCS";
    params["fitpsf_file"] = psf_file;

    // Load the pixels
    GetImagePixelLists(fnum,b);
    dbg<<"\n";

#ifdef OnlyNImages
    if (fnum + 1 >= OnlyNImages) break;
#endif
  }
  int npix=0;
  for (int i=0;i<size();++i) if (pixlist[i].size() > 0) ++npix;
  return npix;
}

MultiShearCatalog::MultiShearCatalog(const ConfigFile& _params) :
  params(_params)
{
  Read();
}

MultiShearCatalog::~MultiShearCatalog()
{
  for(size_t k=0;k<trans.size();++k) {
    if (trans[k]) delete trans[k];
  }
  for(size_t k=0;k<fitpsf.size();++k) {
    if (fitpsf[k]) delete fitpsf[k];
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
    shear_file_list.clear();
    fitpsf_file_list.clear();

    std::string image_filename;
    std::string shear_filename;
    std::string psf_filename;
    params["dist_hdu"] = 2;
    params["dist_method"] = "WCS";
    while (flist >> image_filename >> shear_filename >> psf_filename) 
    {
      image_file_list.push_back(image_filename);
      shear_file_list.push_back(shear_filename);
      fitpsf_file_list.push_back(psf_filename);
      dbg<<"Files are :\n"<<image_filename<<std::endl;
      dbg<<shear_filename<<std::endl;
      dbg<<psf_filename<<std::endl;

      // read transformation between ra/dec and x/y
      params["dist_file"] = image_filename;
      trans.push_back(new Transformation(params));

      // read the psf
      params["fitpsf_file"] = psf_filename;
      fitpsf.push_back(new FittedPSF(params));
      dbg<<"fitpsf["<<fitpsf.size()-1<<"] sigma = "<<fitpsf.back()->GetSigma()<<std::endl;

      // We need the bounds for each component image so that we can 
      // efficiently load only the images we need in each area of the sky 
      // we will be working on.  
      // The easiest way to do this is to load the shear catalog
      // and read off the skybounds from that.  But we don't need to 
      // keep the whole shear catalog in memory.
      params["shear_file"] = shear_filename;
      ShearCatalog scat(params);
      image_skybounds.push_back(scat.skybounds);
    }
  }
  catch (std::exception& e)
  {
    throw ReadError("Error reading from "+file+" -- caught error\n" +
	e.what());
  }
  catch (...)
  {
    throw ReadError("Error reading from "+file+" -- caught unknown error");
  }

  dbg<<"Done reading file lists\n";
  Assert(image_file_list.size() == fitpsf_file_list.size());
  Assert(shear_file_list.size() == fitpsf_file_list.size());
  Assert(image_skybounds.size() == fitpsf_file_list.size());

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

template <class T>
long long MemoryFootprint(const T& x)
{
  return sizeof(x); 
};

template <class T>
long long MemoryFootprint(T*const x)
{
  long long res=sizeof(x);
  res += MemoryFootprint(*x);
  return res;
}

template <class T>
long long MemoryFootprint(const T*const x)
{
  long long res=sizeof(x);
  res += MemoryFootprint(*x);
  return res;
}

template <class T>
long long MemoryFootprint(const std::auto_ptr<T>& x)
{
  long long res=sizeof(x);
  res += MemoryFootprint(*x);
  return res;
}

template <class T, class Alloc>
long long MemoryFootprint(const std::vector<T,Alloc>& x)
{
  long long res=sizeof(x);
  for(size_t i=0;i<x.size();++i) res += MemoryFootprint(x[i]);
  return res;
}

long long MemoryFootprint(const PixelList& x)
{
  long long res=sizeof(x);
  for(size_t i=0;i<x.size();++i) res += MemoryFootprint(x[i]);
  return res;
}

template <class T>
long long MemoryFootprint(const tmv::Vector<T>& x)
{
  long long res=sizeof(x);
  res += x.size() * sizeof(T);
  return res;
}

long long MemoryFootprint(const BVec& x)
{
  long long res=sizeof(x);
  res += x.size() * sizeof(double);
  return res;
}

template <class T>
long long MemoryFootprint(const tmv::Matrix<T>& x)
{
  long long res=sizeof(x);
  res += x.colsize() * x.rowsize() * sizeof(T);
  return res;
}

// Get pixel lists from the file specified in params
void MultiShearCatalog::GetImagePixelLists(int image_index, const Bounds& b)
{
  dbg<<"Start GetImagePixelLists: image_index = "<<image_index<<std::endl;
  dbg<<"trans.size = "<<trans.size()<<std::endl;
  dbg<<"pixlist.size = "<<pixlist.size()<<std::endl;
  dbg<<"Memory usage:\n";
  dbg<<"trans: "<<MemoryFootprint(trans)/1024./1024.<<" MB\n";
  dbg<<"fitpsf: "<<MemoryFootprint(fitpsf)/1024./1024.<<" MB\n";
  dbg<<"input_flags: "<<MemoryFootprint(input_flags)/1024./1024.<<" MB\n";
  dbg<<"nimages_found: "<<MemoryFootprint(nimages_found)/1024./1024.<<" MB\n";
  dbg<<"pixlist: "<<MemoryFootprint(pixlist)/1024./1024.<<" MB\n";
  dbg<<"nimages_gotpix: "<<MemoryFootprint(nimages_gotpix)/1024./1024.<<" MB\n";
  dbg<<"id: "<<MemoryFootprint(id)/1024./1024.<<" MB\n";
  dbg<<"skypos: "<<MemoryFootprint(skypos)/1024./1024.<<" MB\n";
  dbg<<"sky: "<<MemoryFootprint(sky)/1024./1024.<<" MB\n";
  dbg<<"noise: "<<MemoryFootprint(noise)/1024./1024.<<" MB\n";
  dbg<<"flags: "<<MemoryFootprint(flags)/1024./1024.<<" MB\n";
  dbg<<"shear: "<<MemoryFootprint(shear)/1024./1024.<<" MB\n";
  dbg<<"nu: "<<MemoryFootprint(nu)/1024./1024.<<" MB\n";
  dbg<<"cov: "<<MemoryFootprint(cov)/1024./1024.<<" MB\n";
  dbg<<"shape: "<<MemoryFootprint(shape)/1024./1024.<<" MB\n";
  dbg<<"image_indexlist: "<<MemoryFootprint(image_indexlist)/1024./1024.<<" MB\n";
  dbg<<"image_cenlist: "<<MemoryFootprint(image_cenlist)/1024./1024.<<" MB\n";
  dbg<<"image_file_list: "<<MemoryFootprint(image_file_list)/1024./1024.<<" MB\n";
  dbg<<"image_skybounds: "<<MemoryFootprint(image_skybounds)/1024./1024.<<" MB\n";
  dbg<<"shear_file_list: "<<MemoryFootprint(shear_file_list)/1024./1024.<<" MB\n";
  dbg<<"fitpsf_file_list: "<<MemoryFootprint(fitpsf_file_list)/1024./1024.<<" MB\n";

  bool output_dots = params.read("output_dots",false);
  if (output_dots)
  {
    long long totmem = 
      MemoryFootprint(trans) +
      MemoryFootprint(fitpsf) +
      MemoryFootprint(input_flags) +
      MemoryFootprint(nimages_found) +
      MemoryFootprint(pixlist) +
      MemoryFootprint(nimages_gotpix) +
      MemoryFootprint(id) +
      MemoryFootprint(skypos) +
      MemoryFootprint(sky) +
      MemoryFootprint(noise) +
      MemoryFootprint(flags) +
      MemoryFootprint(shear) +
      MemoryFootprint(nu) +
      MemoryFootprint(cov) +
      MemoryFootprint(shape) +
      MemoryFootprint(image_indexlist) +
      MemoryFootprint(image_cenlist) +
      MemoryFootprint(image_file_list) +
      MemoryFootprint(image_skybounds) +
      MemoryFootprint(shear_file_list) +
      MemoryFootprint(fitpsf_file_list);

    std::cerr<<"Using image# "<<image_index;
    std::cerr<<"... Memory Usage in MultiShearCatalog = ";
    std::cerr<<(totmem/1024./1024.)<<" MB";
    //std::cerr<<" ["<<(MemoryFootprint(pixlist)/1024./1024.)<<"]";
    //std::cerr<<" ["<<(MemoryFootprint(shape)/1024./1024.)<<"]";
    std::cerr<<"\n";
  }

  std::auto_ptr<Image<double> > weight_im;
  Image<double> im(params,weight_im);

  int maxi=im.GetMaxI();
  int maxj=im.GetMaxJ();
  xdbg<<"MaxI: "<<maxi<<" MaxJ: "<<maxj<<"\n";

  const Transformation& trans1 = *trans[image_index];
  const FittedPSF& fitpsf1 = *fitpsf[image_index];
  BVec psf1(fitpsf1.GetPSFOrder(), fitpsf1.GetSigma());

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
  Assert(pixlist.size() == skypos.size());
  Assert(image_indexlist.size() == skypos.size());
  Assert(image_cenlist.size() == skypos.size());
  for (size_t i=0; i<size(); ++i) {
#ifdef SKIP_INPUT_FLAGS
    if (flags[i]) continue;
#endif
    Assert(image_indexlist[i].size() == pixlist[i].size());
    Assert(image_cenlist[i].size() == pixlist[i].size());
    if (!b.Includes(skypos[i])) continue;

    // convert ra/dec to x,y in this image

#ifdef UseInverseTransform
    // Figure out a good starting point for the nonlinear solver:
    Position pxy;
    //xdbg<<"skypos = "<<skypos[i]<<std::endl;
    if (!invb.Includes(skypos[i])) {
      //xdbg<<"skypos "<<skypos[i]<<" not in "<<invb<<std::endl;
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
      Assert(i < pixlist.size());
      xdbg<<"pixlist["<<i<<"].size = "<<pixlist[i].size()<<std::endl;
      xdbg<<"image_indexlist["<<i<<"].size = "<<image_indexlist[i].size()<<std::endl;
      xdbg<<"image_cenlist["<<i<<"].size = "<<image_cenlist[i].size()<<std::endl;
      Assert(image_indexlist[i].size() == pixlist[i].size());
      Assert(image_cenlist[i].size() == pixlist[i].size());
      pixlist[i].push_back(PixelList());
      //pixlist[i].back().UseBlockMem();
      xdbg<<"pixlist["<<i<<"].size = "<<pixlist[i].size()<<std::endl;
      GetPixList(
	  im,pixlist[i].back(),pxy,
	  sky[i],noise,gain,weight_im.get(),trans1,max_aperture,flag);
      xdbg<<"Got pixellist, flag = "<<flag<<std::endl;

      // make sure not (edge or < 10 pixels) although edge is already
      // checked above
      if (flag == 0) {
	dbg<<"i = "<<i<<", pixlist.size = "<<pixlist.size()<<std::endl;
	Assert(i < image_indexlist.size());
	Assert(i < image_cenlist.size());
	Assert(i < nimages_gotpix.size());
	image_indexlist[i].push_back(image_index);
	image_cenlist[i].push_back(pxy);
	xdbg<<"image_indexlist["<<i<<"].size = "<<image_indexlist[i].size()<<std::endl;
	xdbg<<"image_cenlist["<<i<<"].size = "<<image_cenlist[i].size()<<std::endl;
	nimages_gotpix[i]++;
      } else {
	input_flags[i] |= flag;
	pixlist[i].pop_back();
	xdbg<<"pixlist["<<i<<"].size = "<<pixlist[i].size()<<std::endl;
      }
      xdbg<<"pixlist["<<i<<"].size = "<<pixlist[i].size()<<std::endl;
      xdbg<<"image_indexlist["<<i<<"].size = "<<image_indexlist[i].size()<<std::endl;
      xdbg<<"image_cenlist["<<i<<"].size = "<<image_cenlist[i].size()<<std::endl;
      Assert(image_indexlist[i].size() == pixlist[i].size());
      Assert(image_cenlist[i].size() == pixlist[i].size());
    } else {
      xdbg<<"x,y not in valid bounds\n";
    }
    xdbg<<"end loop over objects\n";
    xdbg<<"pixlist["<<i<<"].size = "<<pixlist[i].size()<<std::endl;
    xdbg<<"image_indexlist["<<i<<"].size = "<<image_indexlist[i].size()<<std::endl;
    xdbg<<"image_cenlist["<<i<<"].size = "<<image_cenlist[i].size()<<std::endl;
    Assert(image_indexlist[i].size() == pixlist[i].size());
    Assert(image_cenlist[i].size() == pixlist[i].size());
  } // loop over objects
  dbg<<"Done extracting pixel lists\n";
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
    catch (CCfits::FitsException& e)
    {
      throw WriteError("Error writing to "+file+" -- caught error\n" +
	  e.message());
    }
    catch (std::exception& e)
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

  for(size_t i=0;i<size();i++) 
  {
    ra[i] = skypos[i].GetX();
    dec[i] = skypos[i].GetY();
  }
  for(size_t i=0;i<size();i++) 
  {
    shear1[i] = real(shear[i]);
    shear2[i] = imag(shear[i]);
  }
  std::vector<double> cov11(size());
  for(size_t i=0;i<size();i++) 
  {
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
  catch (CCfits::FitsException& e)
  {
    throw ReadError("Error reading from "+file+" -- caught error\n" +
	e.message());
  }
  catch (std::exception& e)
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

