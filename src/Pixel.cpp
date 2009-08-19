
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
  // xcen,ycen are given on a 1-based grid.
  // ie. where the lower left corner pixel is (1,1), rather than (0,0).
  // The easiest way to do this is to just decrease xcen,ycen by 1 each:
  //--xcen; --ycen;
  // === This is now handled by x_offset, y_offset in ReadCatalog

  int i1 = int(floor(xcen-xap-xmin));
  int i2 = int(ceil(xcen+xap-xmin));
  int j1 = int(floor(ycen-yap-ymin));
  int j2 = int(ceil(ycen+yap-ymin));

  if (i1 < 0) { i1 = 0; flag |= EDGE; }
  if (i2 > int(im.GetMaxI())) { i2 = im.GetMaxI(); flag |= EDGE; }
  if (j1 < 0) { j1 = 0; flag |= EDGE; }
  if (j2 > int(im.GetMaxJ())) { j2 = im.GetMaxJ(); flag |= EDGE; }

  double apsq = aperture*aperture;

  // Do this next loop in two passes.  First figure out which 
  // pixels we want to use.  Then we can resize pix to the full size
  // we will need, and go back through and enter the pixels.
  // This saves us a lot of resizing calls in vector, which are
  // both slow and can fragment the memory.
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
	Assert(k < int(pix.size()));
	pix[k++] = Pixel(u,v,I,wt);
      }
    }
  }
  Assert(k == int(pix.size()));
  if (npix < 10) flag |= LT10PIX;
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
