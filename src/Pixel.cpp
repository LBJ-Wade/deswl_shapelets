
#include "Pixel.h"
#include "Params.h"

void GetPixList(
    const Image<double>& im, PixelList& pix,
    const Position cen, double sky, double noise,
    const Image<double>* weight_image, const Transformation& trans,
    double aperture, const ConfigFile& params, long& flag)
{
    double gain = params.read("image_gain",0.);
    double x_offset = params.read("cat_x_offset",0.);
    double y_offset = params.read("cat_y_offset",0.);
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
    double pixscale = sqrt(det); // arsec/pixel
    xdbg<<"pixscale = "<<pixscale<<std::endl;

    // xap,yap are the maximum deviation from the center in x,y
    // such that u'^2+v'^2 = aperture^2
    double xap = aperture * sqrt(D(1,0)*D(1,0) + D(1,1)*D(1,1))/det;
    double yap = aperture * sqrt(D(0,0)*D(0,0) + D(0,1)*D(0,1))/det;
    xdbg<<"aperture = "<<aperture<<std::endl;
    xdbg<<"xap = "<<xap<<", yap = "<<yap<<std::endl;

    int xmin = im.getXMin();
    int ymin = im.getYMin();
    int xmax = im.getXMax();
    int ymax = im.getYMax();

    double xcen = cen.getX();
    double ycen = cen.getY();
    xdbg<<"cen = "<<xcen<<"  "<<ycen<<std::endl;
    xdbg<<"xmin, ymin = "<<xmin<<"  "<<ymin<<std::endl;
    // xcen,ycen are given on a 1-based grid.
    // ie. where the lower left corner pixel is (1,1), rather than (0,0).
    // The easiest way to do this is to just decrease xcen,ycen by 1 each:
    //--xcen; --ycen;
    // This is now handled by parameters x_offset and y_offset, which are
    // set from the parameters: cat_x_offset and cat_y_offset
    xcen -= x_offset;
    ycen -= y_offset;

    if (xcen < xmin || xcen > xmax || ycen < ymin || ycen > ymax) {
        dbg<<"Position "<<xcen<<" , "<<ycen<<
            " does not fall within bounds of image:\n";
        dbg<<xmin<<"  "<<xmax<<"  "<<ymin<<"  "<<ymax<<std::endl;
        dbg<<"returning with no pixels\n";
        pix.resize(0);
        if (!ignore_edges) flag |= EDGE;
        flag |= LT10PIX;
        return;
    }


    int i1 = int(floor(xcen-xap-xmin));
    int i2 = int(ceil(xcen+xap-xmin));
    int j1 = int(floor(ycen-yap-ymin));
    int j2 = int(ceil(ycen+yap-ymin));
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
    int npix = 0;

    double chipx = xmin+i1-xcen;
    double peak = 0.;
    for(int i=i1;i<=i2;++i,chipx+=1.) {
        double chipy = ymin+j1-ycen;
        double u = D(0,0)*chipx+D(0,1)*chipy;
        double v = D(1,0)*chipx+D(1,1)*chipy;
        for(int j=j1;j<=j2;++j,u+=D(0,1),v+=D(1,1)) {
            double rsq = u*u+v*v;
            if (rsq <= apsq) {
                use[i-i1][j-j1] = true;
                ++npix;
                if (im(i,j) > peak) peak = im(i,j);
            }
        }
    }

    xdbg<<"npix = "<<npix<<std::endl;
    pix.resize(npix);

    xdbg<<"pixlist size = "<<npix<<" = "<<npix*sizeof(Pixel)<<
        " bytes = "<<npix*sizeof(Pixel)/1024.<<" KB\n";

    int k=0;
    chipx = xmin+i1-xcen;
    xdbg<<"Bright pixels are:\n";
    peak -= sky;
    for(int i=i1;i<=i2;++i,chipx+=1.) {
        double chipy = ymin+j1-ycen;
        double u = D(0,0)*chipx+D(0,1)*chipy;
        double v = D(1,0)*chipx+D(1,1)*chipy;
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
    npix = pix.size(); // may have changed.
    xdbg<<"npix => "<<npix<<std::endl;
    if (npix < 10) flag |= LT10PIX;
}

double GetLocalSky(
    const Image<double>& bkg, 
    const Position cen, const Transformation& trans, double aperture,
    const ConfigFile& params, long& flag)
{
    double x_offset = params.read("cat_x_offset",0.);
    double y_offset = params.read("cat_y_offset",0.);

    // This function is very similar in structure to the above GetPixList
    // function.  It does the same thing with the distortion and the 
    // aperture and such.  
    // The return value is the mean sky value within the aperture.

    xdbg<<"Start GetLocalSky\n";

    DSmallMatrix22 D;
    trans.getDistortion(cen,D);

    double det = std::abs(D.TMV_det());
    double pixscale = sqrt(det); // arcsec/pixel
    xdbg<<"pixscale = "<<pixscale<<std::endl;

    // xap,yap are the maximum deviation from the center in x,y
    // such that u^2+v^2 = aperture^2
    double xap = aperture / det * 
        sqrt(D(0,0)*D(0,0) + D(0,1)*D(0,1));
    double yap = aperture / det * 
        sqrt(D(1,0)*D(1,0) + D(1,1)*D(1,1));
    xdbg<<"aperture = "<<aperture<<std::endl;
    xdbg<<"xap = "<<xap<<", yap = "<<yap<<std::endl;

    int xmin = bkg.getXMin();
    int ymin = bkg.getYMin();

    double xcen = cen.getX();
    double ycen = cen.getY();
    xdbg<<"cen = "<<xcen<<"  "<<ycen<<std::endl;
    xdbg<<"xmin, ymin = "<<xmin<<"  "<<ymin<<std::endl;
    xcen -= x_offset;
    ycen -= y_offset;

    int i1 = int(floor(xcen-xap-xmin));
    int i2 = int(ceil(xcen+xap-xmin));
    int j1 = int(floor(ycen-yap-ymin));
    int j2 = int(ceil(ycen+yap-ymin));
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

    double mean_sky = 0.;
    int npix = 0;

    double chipx = xmin+i1-xcen;
    for(int i=i1;i<=i2;++i,chipx+=1.) {
        double chipy = ymin+j1-ycen;
        double u = D(0,0)*chipx+D(0,1)*chipy;
        double v = D(1,0)*chipx+D(1,1)*chipy;
        for(int j=j1;j<=j2;++j,u+=D(0,1),v+=D(1,1)) {
            // u,v are in arcsec
            double rsq = u*u + v*v;
            if (rsq <= apsq) {
                mean_sky += bkg(i,j);
                ++npix;
            }
        }
    }

    xdbg<<"npix = "<<npix<<std::endl;
    if (npix == 0) { flag |= BKG_NOPIX; return 0.; }

    mean_sky /= npix;
    xdbg<<"meansky = "<<mean_sky<<std::endl;
    return mean_sky;
}

void GetSubPixList(
    PixelList& pix, const PixelList& allpix,
    std::complex<double> cen_offset, std::complex<double> shear,
    double aperture, double inner_fake_ap, double outer_fake_ap,
    const ConfigFile& params, long& flag)
{
    bool use_fake = outer_fake_ap > aperture;

    // Select a subset of allpix that are within the given aperture
    const int ntot = allpix.size();
    xdbg<<"Start GetSubPixList\n";
    xdbg<<"allpix has "<<ntot<<" objects\n";
    xdbg<<"new aperture = "<<aperture<<std::endl;
    xdbg<<"cen_offset = "<<cen_offset<<std::endl;
    xdbg<<"shear = "<<shear<<std::endl;
    xdbg<<"use_fake = "<<use_fake<<std::endl;
    xdbg<<"inner_fake_ap = "<<inner_fake_ap<<std::endl;
    xdbg<<"outer_fake_ap = "<<outer_fake_ap<<std::endl;

    double normg = norm(shear);
    double g1 = real(shear);
    double g2 = imag(shear);
    double apsq = aperture*aperture;
    double inner_fake_apsq = inner_fake_ap*inner_fake_ap;
    double outer_fake_apsq = outer_fake_ap*outer_fake_ap;

    // Do this next loop in two passes.  First figure out which 
    // pixels we want to use.  Then we can resize pix to the full size
    // we will need, and go back through and enter the pixels.
    // This saves us a lot of resizing calls in vector, which are
    // both slow and can fragment the memory.
    std::vector<bool> use(ntot,false);
    std::vector<bool> fake(ntot,false);
    int npix = 0;

    double peak = 0.;
    for(int i=0;i<ntot;++i) {
        std::complex<double> z = allpix[i].getPos() - cen_offset;
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
            ++npix;
            if (allpix[i].getFlux() > peak) peak = allpix[i].getFlux();
        } else if (use_fake && 
                   rsq >= inner_fake_apsq &&
                   rsq <= outer_fake_apsq) {
            fake[i] = true;
            ++npix;
        }
    }

    xdbg<<"npix = "<<npix<<std::endl;
    xdbg<<"regular pixels = "<<std::count(use.begin(),use.end(),true)<<std::endl;
    xdbg<<"fake pixels = "<<std::count(fake.begin(),fake.end(),true)<<std::endl;
    pix.resize(npix);

    xdbg<<"pixlist size = "<<npix<<" = "<<npix*sizeof(Pixel)<<
        " bytes = "<<npix*sizeof(Pixel)/1024.<<" KB\n";

    int k=0;
    xdbg<<"Bright pixels are:\n";
    for(int i=0;i<ntot;++i) {
        if(use[i]) {
            Pixel p = allpix[i];
            p.setPos(p.getPos() - cen_offset);
            pix[k++] = p;
            if (p.getFlux() > peak / 10.) {
                xdbg<<p.getPos()<<"  "<<p.getFlux()<<std::endl;
            }
        } else if (fake[i]) {
            Pixel p = allpix[i];
            p.setPos(p.getPos() - cen_offset);
            p.setFlux(0.);
            // Keep same noise for fake pixels
            pix[k++] = p;
        }
    }
    Assert(k == int(pix.size()));

    if (npix < 10) flag |= LT10PIX;
}

