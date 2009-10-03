
#include "Pixel.h"
#include "Params.h"
#include "dbg.h"

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
