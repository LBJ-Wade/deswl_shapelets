
#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <algorithm>
#include <cstdlib>

// Weight() calculates an appropriate weight for each galaxy.
// The important part of the weight function is really the dr parameter.
// The responsivity will be Sum(dr) / Sum(w).
// There are several different weight schemes that are available,
// currently specified by #if blocks (TODO: Make this a runtime parameter).
// See BJ02 Eqn 5-33.
static double Weight(
    std::complex<double> e, double vare, double shapenoise, double& dr)
{
    // Impose a minimum on the variance to keep them from dominating
    // the weights.
    if (vare < 1.e-3) vare = 1.e-3;

    double esq = std::norm(e);
    double f = shapenoise / (vare+shapenoise);
    double k0 = f*vare;
    double k1 = f*f;

#if 1
    // No Weight
    double w = 1.;
    double edwde = 0.;
#endif

#if 0
    // Simple inverse variance weight
    double w = 1./vare;
    double edwde = 0.;
#endif

#if 0
    // Better weight if lensing shear is small 
    double w = 1./sqrt(esq + 2.25 * vare);
    //if (vare > 1.) w = 0.;
    double edwde = -w*w*w*esq;
#endif

    dr = w * (1-k0-0.5*k1*esq) + 0.5*edwde*(1.-k0-k1*esq);
    return w;
}

// ReadData() reads in the columns of a shear file as output by measureshears.
// This is currently just the Ascii output, which is written as _shear.dat.
static void ReadData(
    std::istream& fin, 
    std::vector<long>& idVec,
    std::vector<double>& xVec,
    std::vector<double>& yVec,
    std::vector<double>& raVec,
    std::vector<double>& decVec,
    std::vector<std::complex<double> >& gVec,
    std::vector<double>& vargVec)
{
    std::string line;
    while (getline(fin,line)) {
        std::istringstream strin(line);
        std::string junk;
        long id,flag;
        double x,y,ra,dec,g1,g2,nu;
        strin >> 
            // id  x  y  sky  noise  flag  ra  dec  g1  g2  nu ...
            id >> x >> y >> junk >> junk >> 
            flag >> ra >> dec >> g1 >> g2 >> nu;
        if (flag & ~0xf8000) continue;
        idVec.push_back(id);
        xVec.push_back(x);
        yVec.push_back(y);
        raVec.push_back(ra);
        decVec.push_back(dec);
        gVec.push_back(std::complex<double>(g1,g2));
        // I think this is more accurate for now than the covariance estimates,
        // but eventually we might want to switch to trace(cov)
        double varg = 4./(nu*nu);
        vargVec.push_back(varg);
    }
}

int main(int argc, char **argv)
{
    if (argc != 3) {
        std::cerr<<"Usage: shearave infile outfile\n";
        exit(1);
    }

    //
    // Open files
    //

    std::string infile = argv[1];
    std::string outfile = argv[2];
    std::ifstream fin(infile.c_str());
    if (!fin) {
        std::cerr<<"Unable to open input file "<<infile<<std::endl;
        exit(1);
    }
    std::ofstream fout(outfile.c_str());
    if (!fout) {
        std::cerr<<"Unable to open output file "<<outfile<<std::endl;
        exit(1);
    }

    //
    // Read data
    //

    std::vector<long> idVec;
    std::vector<double> xVec;
    std::vector<double> yVec;
    std::vector<double> raVec;
    std::vector<double> decVec;
    std::vector<std::complex<double> > gVec;
    std::vector<double> vargVec;

    ReadData(fin,idVec,xVec,yVec,raVec,decVec,gVec,vargVec);

    // 
    // Calculate shape noise
    //

    const int nGal = idVec.size();
    std::vector<double> vareVec(nGal);
    std::vector<std::complex<double> > eVec(nGal);
    double shapenoise=0., sumw=0.;

    for(int i=0;i<nGal;++i) {
        std::complex<double> g = gVec[i];
        double varg = vargVec[i];
        double gsq = std::norm(g);
        double absg = sqrt(gsq);
        double eta = atanh(absg)*2.;
        double e = tanh(eta);
        double esq = e*e;
        // e = tanh(2*atanh(g))
        // de/dg = 2(1-e^2)/(1-g^2)
        double dedg = 2.*(1.-esq)*(1.-gsq);
        double vare = dedg*dedg*varg;
        eVec[i] = e/absg * g;
        vareVec[i] = vare;
        // TODO: The 0.1 in the next line should be a parameter.
        if (vare > 0.1) continue; // Throw out sigma_e > 0.3
        double dr;
        double galw = Weight(e,vare,0.,dr);
        shapenoise += (esq - vare)*galw;
        sumw += galw;
    }
    shapenoise /= 2.*sumw; // 2 so shapenoise per component
    std::cout<<"sigma_SN = "<<shapenoise<<std::endl;

    //
    // Calculate responsivity
    //

    double resp=0.,sumw2e2=0.;
    sumw = 0.;
    std::complex<double> sumwe = 0.;
    for (int i=0; i<nGal; ++i) {
        std::complex<double> galpos(raVec[i],decVec[i]);
        std::complex<double> e = eVec[i];
        double esq = std::norm(e);
        double vare = vareVec[i];
        if (vare > 0.1) continue;
        double dr;
        double galw = Weight(e,vare,shapenoise,dr);
        resp += dr;
        sumw += galw;
        sumwe += galw*e;
        sumw2e2 += galw*galw*esq;
    }
    resp /= sumw;
    std::cout<<"resp = "<<resp<<", ngal = "<<nGal<<std::endl;

    // 
    // Calculate estimate of applied shear
    //

    std::complex<double> e = sumwe / sumw;
    e /= resp;
    double vare = sumw2e2 / (2.*sumw*sumw);
    double sige = sqrt(vare);
    // Now convert back to gamma:
    double abse = std::abs(e);
    double eta = atanh(abse);
    double absgamma = tanh(eta/2);
    std::complex<double> gamma = absgamma/abse * e;
    std::cout<<"e = "<<e<<" +- "<<sige<<std::endl;
    std::cout<<"gamma = "<<gamma<<std::endl;

    //
    // Write out responsivity-corrected catalog
    //

    for (int i=0; i<nGal; i++) {
        fout.precision(12);
        fout<<
            xVec[i]<<"  "<<yVec[i]<<"  "<<
            raVec[i]<<"  "<<decVec[i]<<"  ";
        fout.precision(4);
        fout<<
            real(gVec[i])/resp<<"  "<<imag(gVec[i])/resp<<"  "<<
            sqrt(vargVec[i])/resp<<std::endl;
    }
    fout.close();

    return 0;
}

