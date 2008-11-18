
#include "Pixel.h"
#include "dbg.h"

void GetPixList(const Image<double>& im, std::vector<Pixel>& pix,
    const Position cen, double sky, double noise, double gain,
    const Image<double>* wt_im, const Transformation& trans,
    double aperture, int& flag)
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
  --xcen; --ycen;

  int i1 = int(floor(xcen-xap-xmin));
  int i2 = int(ceil(xcen+xap-xmin));
  int j1 = int(floor(ycen-yap-ymin));
  int j2 = int(ceil(ycen+yap-ymin));

  if (i1 < 0) { i1 = 0; flag |= EDGE; }
  if (i2 > int(im.GetMaxI())) { i2 = im.GetMaxI(); flag |= EDGE; }
  if (j1 < 0) { j1 = 0; flag |= EDGE; }
  if (j2 > int(im.GetMaxJ())) { j2 = im.GetMaxJ(); flag |= EDGE; }

  double apsq = aperture*aperture;

  double chipx = xmin+i1-xcen;
  for(int i=i1;i<=i2;i++,chipx+=1.) {
    double chipy = ymin+j1-ycen;
    double u = D(0,0)*chipx+D(0,1)*chipy;
    double v = D(1,0)*chipx+D(1,1)*chipy;
    for(int j=j1;j<=j2;j++,u+=D(0,1),v+=D(1,1)) {
      // u,v are in arcsec
      double rsq = u*u + v*v;
      xxdbg<<"i,j,u,v,rsq = "<<i<<','<<j<<"  "<<u<<','<<v<<"  "<<rsq<<std::endl;
      if (rsq > apsq) continue;
      xxdbg<<rsq<<" < "<<apsq<<std::endl;
      double I = im(i,j)-sky;
      double wt;
      if (wt_im) {
	wt = (*wt_im)(i,j);
      } else {
	double var = noise;
	if (gain != 0.) var += im(i,j)/gain;
	wt = 1./var;
      }
      pix.push_back(Pixel(u,v,I,wt));
    }
  }
  xdbg<<"done: npix = "<<pix.size()<<std::endl;
  if (pix.size() < 10) flag |= LT10PIX;
}

