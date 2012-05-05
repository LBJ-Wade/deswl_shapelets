
#include "dbg.h"
#include "Pixel.h"
#include "Params.h"

void GetPixList(
    const Image<double>& im, PixelList& pix,
    const Position cen, double sky, double noise,
    const Image<double>* weight_image, const Transformation& trans,
    double aperture, const ConfigFile& params, long& flag)
{
    double gain = params.read("image_gain",0.);
    double xOffset = params.read("cat_x_offset",0.);
    double yOffset = params.read("cat_y_offset",0.);
    bool ignore_edges = params.read("ignore_edges",false);

    xdbg<<"Start GetPixList\n";
    if (weight_image) {
        xdbg<<"Using weight image for pixel noise.\n";
    } else {
        xdbg<<"noise = "<<noise<<std::endl;
        xdbg<<"gain = "<<gain<<std::endl;
    }

    DSmallMatrix22 D;
    trans.getDistortion(cen,D);
    // This gets us from chip coordinates to sky coordinates:
    // ( u ) = D ( x )
    // ( v )     ( y )
    xdbg<<"D = "<<D<<std::endl;
    double det = std::abs(D.TMV_det());
    double pixScale = sqrt(det); // arsec/pixel
    xdbg<<"pixscale = "<<pixScale<<std::endl;

    // xAp,yAp are the maximum deviation from the center in x,y
    // such that u'^2+v'^2 = aperture^2
    double xAp = aperture * sqrt(D(1,0)*D(1,0) + D(1,1)*D(1,1))/det;
    double yAp = aperture * sqrt(D(0,0)*D(0,0) + D(0,1)*D(0,1))/det;
    xdbg<<"aperture = "<<aperture<<std::endl;
    xdbg<<"xap = "<<xAp<<", yap = "<<yAp<<std::endl;

    int xMin = im.getXMin();
    int yMin = im.getYMin();
    int xMax = im.getXMax();
    int yMax = im.getYMax();

    double xCen = cen.getX();
    double yCen = cen.getY();
    xdbg<<"cen = "<<xCen<<"  "<<yCen<<std::endl;
    xdbg<<"xmin, ymin = "<<xMin<<"  "<<yMin<<std::endl;
    // xCen,yCen are given on a 1-based grid.
    // ie. where the lower left corner pixel is (1,1), rather than (0,0).
    // The easiest way to do this is to just decrease xCen,yCen by 1 each:
    //--xCen; --yCen;
    // This is now handled by parameters xOffset and yOffset, which are
    // set from the parameters: cat_x_offset and cat_y_offset
    xCen -= xOffset;
    yCen -= yOffset;

    if (xCen < xMin || xCen > xMax || yCen < yMin || yCen > yMax) {
        dbg<<"Position "<<xCen<<" , "<<yCen<<
            " does not fall within bounds of image:\n";
        dbg<<xMin<<"  "<<xMax<<"  "<<yMin<<"  "<<yMax<<std::endl;
        dbg<<"returning with no pixels\n";
        pix.resize(0);
        if (!ignore_edges) flag |= EDGE;
        flag |= LT10PIX;
        return;
    }


    int i1 = int(floor(xCen-xAp-xMin));
    int i2 = int(ceil(xCen+xAp-xMin));
    int j1 = int(floor(yCen-yAp-yMin));
    int j2 = int(ceil(yCen+yAp-yMin));
    xdbg<<"i1,i2,j1,j2 = "<<i1<<','<<i2<<','<<j1<<','<<j2<<std::endl;

    if (i1 < 0) { i1 = 0; if (!ignore_edges) flag |= EDGE; }
    if (i2 > int(im.getMaxI())) { i2 = im.getMaxI(); if (!ignore_edges) flag |= EDGE; }
    if (j1 < 0) { j1 = 0; if (!ignore_edges) flag |= EDGE; }
    if (j2 > int(im.getMaxJ())) { j2 = im.getMaxJ(); if (!ignore_edges) flag |= EDGE; }
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
    std::vector<std::vector<bool> > use(
        i2-i1+1,std::vector<bool>(j2-j1+1,false));
    int nPix = 0;

    double chipX = xMin+i1-xCen;
    double peak = 0.;
    for(int i=i1;i<=i2;++i,chipX+=1.) {
        double chipY = yMin+j1-yCen;
        double u = D(0,0)*chipX+D(0,1)*chipY;
        double v = D(1,0)*chipX+D(1,1)*chipY;
        for(int j=j1;j<=j2;++j,u+=D(0,1),v+=D(1,1)) {
            double rsq = u*u+v*v;
            if (rsq <= apsq) {
                use[i-i1][j-j1] = true;
                ++nPix;
                if (im(i,j) > peak) peak = im(i,j);
            }
        }
    }

    xdbg<<"npix = "<<nPix<<std::endl;
    pix.resize(nPix);

    xdbg<<"pixlist size = "<<nPix<<" = "<<nPix*sizeof(Pixel)<<
        " bytes = "<<nPix*sizeof(Pixel)/1024.<<" KB\n";

    int k=0;
    chipX = xMin+i1-xCen;
    xdbg<<"Bright pixels are:\n";
    peak -= sky;
    for(int i=i1;i<=i2;++i,chipX+=1.) {
        double chipY = yMin+j1-yCen;
        double u = D(0,0)*chipX+D(0,1)*chipY;
        double v = D(1,0)*chipX+D(1,1)*chipY;
        for(int j=j1;j<=j2;++j,u+=D(0,1),v+=D(1,1)) {
            if (use[i-i1][j-j1]) {
                double flux = im(i,j)-sky;
                double inverseVariance;
                if (weight_image) {
                    inverseVariance = (*weight_image)(i,j);
                } else {
                    double var = noise;
                    if (gain != 0.) var += im(i,j)/gain;
                    inverseVariance = 1./var;
                }
                if (inverseVariance > 0.0) {
                    double inverseSigma = sqrt(inverseVariance);
                    Assert(k < int(pix.size()));
                    Pixel p(u,v,flux,inverseSigma);
                    pix[k++] = p;
                    if (flux > peak / 10.) {
                        xdbg<<p.getPos()<<"  "<<p.getFlux()<<std::endl;
                    }
                }
            }
        }
    }
    Assert(k <= int(pix.size()));
    // Not necessarily == because we skip pixels with 0.0 variance
    pix.resize(k);
    Assert(k == int(pix.size()));
    nPix = pix.size(); // may have changed.
    xdbg<<"npix => "<<nPix<<std::endl;
    if (nPix < 10) flag |= LT10PIX;
}

double GetLocalSky(
    const Image<double>& bkg, 
    const Position cen, const Transformation& trans, double aperture,
    const ConfigFile& params, long& flag)
{
    double xOffset = params.read("cat_x_offset",0.);
    double yOffset = params.read("cat_y_offset",0.);

    // This function is very similar in structure to the above GetPixList
    // function.  It does the same thing with the distortion and the 
    // aperture and such.  
    // The return value is the mean sky value within the aperture.

    xdbg<<"Start GetLocalSky\n";

    DSmallMatrix22 D;
    trans.getDistortion(cen,D);

    double det = std::abs(D.TMV_det());
    double pixScale = sqrt(det); // arcsec/pixel
    xdbg<<"pixscale = "<<pixScale<<std::endl;

    // xAp,yAp are the maximum deviation from the center in x,y
    // such that u^2+v^2 = aperture^2
    double xAp = aperture / det * 
        sqrt(D(0,0)*D(0,0) + D(0,1)*D(0,1));
    double yAp = aperture / det * 
        sqrt(D(1,0)*D(1,0) + D(1,1)*D(1,1));
    xdbg<<"aperture = "<<aperture<<std::endl;
    xdbg<<"xap = "<<xAp<<", yap = "<<yAp<<std::endl;

    int xMin = bkg.getXMin();
    int yMin = bkg.getYMin();

    double xCen = cen.getX();
    double yCen = cen.getY();
    xdbg<<"cen = "<<xCen<<"  "<<yCen<<std::endl;
    xdbg<<"xmin, ymin = "<<xMin<<"  "<<yMin<<std::endl;
    xCen -= xOffset;
    yCen -= yOffset;

    int i1 = int(floor(xCen-xAp-xMin));
    int i2 = int(ceil(xCen+xAp-xMin));
    int j1 = int(floor(yCen-yAp-yMin));
    int j2 = int(ceil(yCen+yAp-yMin));
    xdbg<<"i1,i2,j1,j2 = "<<i1<<','<<i2<<','<<j1<<','<<j2<<std::endl;
    if (i1 < 0) { i1 = 0; }
    if (i2 > int(bkg.getMaxI())) { i2 = bkg.getMaxI(); }
    if (j1 < 0) { j1 = 0; }
    if (j2 > int(bkg.getMaxJ())) { j2 = bkg.getMaxJ(); }
    xdbg<<"i1,i2,j1,j2 => "<<i1<<','<<i2<<','<<j1<<','<<j2<<std::endl;

    double apsq = aperture*aperture;

    xdbg<<"nx = "<<i2-i1+1<<std::endl;
    xdbg<<"ny = "<<j2-j1+1<<std::endl;
    Assert(i2-i1+1 >= 0);
    Assert(j2-j1+1 >= 0);

    double meanSky = 0.;
    int nPix = 0;

    double chipX = xMin+i1-xCen;
    for(int i=i1;i<=i2;++i,chipX+=1.) {
        double chipY = yMin+j1-yCen;
        double u = D(0,0)*chipX+D(0,1)*chipY;
        double v = D(1,0)*chipX+D(1,1)*chipY;
        for(int j=j1;j<=j2;++j,u+=D(0,1),v+=D(1,1)) {
            // u,v are in arcsec
            double rsq = u*u + v*v;
            if (rsq <= apsq) {
                meanSky += bkg(i,j);
                ++nPix;
            }
        }
    }

    xdbg<<"nPix = "<<nPix<<std::endl;
    if (nPix == 0) { flag |= BKG_NOPIX; return 0.; }

    meanSky /= nPix;
    xdbg<<"meansky = "<<meanSky<<std::endl;
    return meanSky;
}

void GetSubPixList(
    PixelList& pix, const PixelList& allPix,
    std::complex<double> cen_offset, std::complex<double> shear,
    double aperture, long& flag)
{
    // Select a subset of allPix that are within the given aperture
    const int ntot = allPix.size();
    xdbg<<"Start GetSubPixList\n";
    xdbg<<"allPix has "<<ntot<<" objects\n";
    xdbg<<"new aperture = "<<aperture<<std::endl;
    xdbg<<"cen_offset = "<<cen_offset<<std::endl;
    xdbg<<"shear = "<<shear<<std::endl;

    double normg = norm(shear);
    double g1 = real(shear);
    double g2 = imag(shear);
    double apsq = aperture*aperture;

    // Do this next loop in two passes.  First figure out which 
    // pixels we want to use.  Then we can resize pix to the full size
    // we will need, and go back through and enter the pixels.
    // This saves us a lot of resizing calls in vector, which are
    // both slow and can fragment the memory.
    std::vector<bool> use(ntot,false);
    int nPix = 0;

    double peak = 0.;
    for(int i=0;i<ntot;++i) {
        std::complex<double> z = allPix[i].getPos() - cen_offset;
        double u = real(z);
        double v = imag(z);
        // (1 + |g|^2) (u^2+v^2) - 2g1 (u^2-v^2) - 2g2 (2uv)
        // u,v are in arcsec
        double usq = u*u;
        double vsq = v*v;
        double rsq = (1.+normg)*(usq+vsq) - 2.*g1*(usq-vsq) - 4.*g2*u*v;
        rsq /= (1.-normg);
        if (rsq <= apsq) {
            //xdbg<<"u,v = "<<u<<','<<v<<"  rsq = "<<rsq<<std::endl;
            use[i] = true;
            ++nPix;
            if (allPix[i].getFlux() > peak) peak = allPix[i].getFlux();
        }
    }

    xdbg<<"npix = "<<nPix<<std::endl;
    pix.resize(nPix);

    xdbg<<"pixlist size = "<<nPix<<" = "<<nPix*sizeof(Pixel)<<
        " bytes = "<<nPix*sizeof(Pixel)/1024.<<" KB\n";

    int k=0;
    xdbg<<"Bright pixels are:\n";
    for(int i=0;i<ntot;++i) if(use[i]) {
        Pixel p = allPix[i];
        p.setPos(p.getPos() - cen_offset);
        pix[k++] = p;
        if (p.getFlux() > peak / 10.) {
            xdbg<<p.getPos()<<"  "<<p.getFlux()<<std::endl;
        }
    }
    Assert(k == int(pix.size()));

    if (nPix < 10) flag |= LT10PIX;
}

