
#include "Pixel.h"
#include "Params.h"
#include "dbg.h"

void GetPixList(const Image<double>& im, PixelList& pix,
    const Position cen, double sky, double noise, double gain,
    const Image<double>* wt_im, const Transformation& trans,
    double aperture, long& flag)
{
  xdbg<<"Start GetPixList\n";
  if (wt_im) {
    xdbg<<"Using weight image for pixel noise.\n";
  } else {
    xdbg<<"noise = "<<noise<<std::endl;
    xdbg<<"gain = "<<gain<<std::endl;
  }

  tmv::SmallMatrix<double,2,2> D;
  trans.GetDistortion(cen,D);

  double det = std::abs(D.Det());
  double pixscale = sqrt(det); // arcsec/pixel
  xdbg<<"pixscale = "<<pixscale<<std::endl;

  // xap,yap are the maximum deviation from the center in x,y
  // such that u^2+v^2 = aperture^2
  double xap = aperture / det * sqrt(D(0,0)*D(0,0) + D(0,1)*D(0,1));
  double yap = aperture / det * sqrt(D(1,0)*D(1,0) + D(1,1)*D(1,1));
  xdbg<<"aperture = "<<aperture<<std::endl;
  xdbg<<"xap = "<<xap<<", yap = "<<yap<<std::endl;

  int xmin = im.GetXMin();
  int ymin = im.GetYMin();

  double xcen = cen.GetX();
  double ycen = cen.GetY();
  xdbg<<"cen = "<<xcen<<"  "<<ycen<<std::endl;
  xdbg<<"xmin, ymin = "<<xmin<<"  "<<ymin<<std::endl;
  // xcen,ycen are given on a 1-based grid.
  // ie. where the lower left corner pixel is (1,1), rather than (0,0).
  // The easiest way to do this is to just decrease xcen,ycen by 1 each:
  //--xcen; --ycen;
  // === This is now handled by x_offset, y_offset in ReadCatalog

  int i1 = int(floor(xcen-xap-xmin));
  int i2 = int(ceil(xcen+xap-xmin));
  int j1 = int(floor(ycen-yap-ymin));
  int j2 = int(ceil(ycen+yap-ymin));
  xdbg<<"i1,i2,j1,j2 = "<<i1<<','<<i2<<','<<j1<<','<<j2<<std::endl;

  if (i1 < 0) { i1 = 0; flag |= EDGE; }
  if (i2 > int(im.GetMaxI())) { i2 = im.GetMaxI(); flag |= EDGE; }
  if (j1 < 0) { j1 = 0; flag |= EDGE; }
  if (j2 > int(im.GetMaxJ())) { j2 = im.GetMaxJ(); flag |= EDGE; }
  xdbg<<"i1,i2,j1,j2 => "<<i1<<','<<i2<<','<<j1<<','<<j2<<std::endl;

  double apsq = aperture*aperture;

  // Do this next loop in two passes.  First figure out which 
  // pixels we want to use.  Then we can resize pix to the full size
  // we will need, and go back through and enter the pixels.
  // This saves us a lot of resizing calls in vector, which are
  // both slow and can fragment the memory.
  xdbg<<"nx = "<<i2-i1+1<<std::endl;
  xdbg<<"ny = "<<j2-j1+1<<std::endl;
  Assert(i2-i1+1 >= 0);
  Assert(j2-j1+1 >= 0);
  std::vector<std::vector<bool> > usepix(i2-i1+1,
      std::vector<bool>(j2-j1+1,false));
  int npix = 0;

  double chipx = xmin+i1-xcen;
  for(int i=i1;i<=i2;i++,chipx+=1.) {
    double chipy = ymin+j1-ycen;
    double u = D(0,0)*chipx+D(0,1)*chipy;
    double v = D(1,0)*chipx+D(1,1)*chipy;
    for(int j=j1;j<=j2;j++,u+=D(0,1),v+=D(1,1)) {
      // u,v are in arcsec
      double rsq = u*u + v*v;
      if (rsq <= apsq) 
      {
	usepix[i-i1][j-j1] = true;
	npix++;
      }
    }
  }

  xdbg<<"npix = "<<npix<<std::endl;
  pix.resize(npix);

  xdbg<<"pixlist size = "<<npix<<" = "<<npix*sizeof(Pixel)<<" bytes = "<<npix*sizeof(Pixel)/1024.<<" KB\n";
  int k=0;

  chipx = xmin+i1-xcen;
  for(int i=i1;i<=i2;i++,chipx+=1.) {
    double chipy = ymin+j1-ycen;
    double u = D(0,0)*chipx+D(0,1)*chipy;
    double v = D(1,0)*chipx+D(1,1)*chipy;
    for(int j=j1;j<=j2;j++,u+=D(0,1),v+=D(1,1)) {
      if (usepix[i-i1][j-j1]) {
	double I = im(i,j)-sky;
	double wt;
	if (wt_im) {
	  wt = (*wt_im)(i,j);
	} else {
	  double var = noise;
	  if (gain != 0.) var += im(i,j)/gain;
	  wt = 1./var;
	}
	if (wt > 0.0) 
	{
	  Assert(k < int(pix.size()));
	  pix[k++] = Pixel(u,v,I,wt);
	  //if (xmin > 0 || ymin > 0)
	    //xdbg<<k-1<<"  "<<i<<"  "<<j<<"  "<<u<<"  "<<v<<"  "<<I<<"  "<<wt<<std::endl;
	}
      }
    }
  }
  Assert(k <= int(pix.size()));
  // Not necessarily == because we skip pixels with 0.0 variance
  pix.resize(k);
  //pix.erase(pix.begin()+k,pix.end());
  Assert(k == int(pix.size()));
  npix = pix.size(); // may have changed.
  xdbg<<"npix => "<<npix<<std::endl;
  if (npix < 10) flag |= LT10PIX;
}

double GetLocalSky(const Image<float>& bkg, 
    const Position cen, const Transformation& trans,
    double aperture, long& flag)
{
  // This function is very similar in structure to the above GetPixList
  // function.  It does the same thing with the distortion and the 
  // aperture and such.  
  // The return value is the mean sky value within the aperture.

  xdbg<<"Start GetLocalSky\n";

  tmv::SmallMatrix<double,2,2> D;
  trans.GetDistortion(cen,D);

  double det = std::abs(D.Det());
  double pixscale = sqrt(det); // arcsec/pixel
  xdbg<<"pixscale = "<<pixscale<<std::endl;

  // xap,yap are the maximum deviation from the center in x,y
  // such that u^2+v^2 = aperture^2
  double xap = aperture / det * sqrt(D(0,0)*D(0,0) + D(0,1)*D(0,1));
  double yap = aperture / det * sqrt(D(1,0)*D(1,0) + D(1,1)*D(1,1));
  xdbg<<"aperture = "<<aperture<<std::endl;
  xdbg<<"xap = "<<xap<<", yap = "<<yap<<std::endl;

  int xmin = bkg.GetXMin();
  int ymin = bkg.GetYMin();

  double xcen = cen.GetX();
  double ycen = cen.GetY();
  // xcen,ycen are given on a 1-based grid.
  // ie. where the lower left corner pixel is (1,1), rather than (0,0).
  // The easiest way to do this is to just decrease xcen,ycen by 1 each:
  //--xcen; --ycen;
  // === This is now handled by x_offset, y_offset in ReadCatalog

  int i1 = int(floor(xcen-xap-xmin));
  int i2 = int(ceil(xcen+xap-xmin));
  int j1 = int(floor(ycen-yap-ymin));
  int j2 = int(ceil(ycen+yap-ymin));
  xdbg<<"i1,i2,j1,j2 = "<<i1<<','<<i2<<','<<j1<<','<<j2<<std::endl;
  if (i1 < 0) { i1 = 0; }
  if (i2 > int(bkg.GetMaxI())) { i2 = bkg.GetMaxI(); }
  if (j1 < 0) { j1 = 0; }
  if (j2 > int(bkg.GetMaxJ())) { j2 = bkg.GetMaxJ(); }
  xdbg<<"i1,i2,j1,j2 => "<<i1<<','<<i2<<','<<j1<<','<<j2<<std::endl;

  double apsq = aperture*aperture;

  xdbg<<"nx = "<<i2-i1+1<<std::endl;
  xdbg<<"ny = "<<j2-j1+1<<std::endl;
  Assert(i2-i1+1 >= 0);
  Assert(j2-j1+1 >= 0);

  double meansky = 0.;
  int npix = 0;

  double chipx = xmin+i1-xcen;
  for(int i=i1;i<=i2;i++,chipx+=1.) {
    double chipy = ymin+j1-ycen;
    double u = D(0,0)*chipx+D(0,1)*chipy;
    double v = D(1,0)*chipx+D(1,1)*chipy;
    for(int j=j1;j<=j2;j++,u+=D(0,1),v+=D(1,1)) {
      // u,v are in arcsec
      double rsq = u*u + v*v;
      if (rsq <= apsq) 
      {
	meansky += bkg(i,j);
	npix++;
      }
    }
  }

  xdbg<<"npix = "<<npix<<std::endl;
  if (npix == 0) { flag |= BKG_NOPIX; return 0.; }
  
  meansky /= npix;
  xdbg<<"meansky = "<<meansky<<std::endl;
  return meansky;
}

void GetSubPixList(PixelList& pix, const PixelList& allpix,
    double aperture, long& flag)
{
  // Select a subset of allpix that are within the given aperture
  xdbg<<"Start GetSubPixList\n";
  xdbg<<"allpix has "<<allpix.size()<<" objects\n";
  xdbg<<"new apertur size is "<<aperture<<std::endl;

  double apsq = aperture*aperture;

  pix.clear();
  for(size_t i=0;i<allpix.size();++i) 
  {
    xxdbg<<"allpix["<<i<<"] = "<<allpix[i].z<<std::endl;
    if (std::norm(allpix[i].z) < apsq)
    {
      pix.push_back(allpix[i]);
      xxdbg<<"added to pix\n";
    } 
    else 
    {
      xxdbg<<"not in aperture\n";
    }

  }
  xdbg<<"done: npix = "<<pix.size()<<std::endl;
  if (pix.size() < 10) flag |= LT10PIX;
}
