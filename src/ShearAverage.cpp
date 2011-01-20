
#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <algorithm>
#include <cstdlib>

#include "Params.h"

#if 1
static const long ok_flags = (
    EDGE |
    SHEAR_REDUCED_ORDER | 
    SHAPE_REDUCED_ORDER | 
    SHAPE_POOR_FIT |
    SHAPE_LOCAL_MIN |
    SHAPE_BAD_FLUX );
#else
static const long ok_flags = 0;
#endif


// Weight() calculates an appropriate weight for each galaxy.
// The important part of the weight function is really the dr parameter.
// The responsivity will be Sum(dr) / Sum(w).
// There are several different weight schemes that are available,
// currently specified by #if blocks (TODO: Make this a runtime parameter).
// See BJ02 Eqn 5-33.
static double Weight(
    std::complex<double> g, double varg, double shapenoise, double& dr)
{
    // Impose a minimum on the variance to keep them from dominating
    // the weights.
    if (varg < 1.e-3) varg = 1.e-3;

    double gsq = std::norm(g);
    double f = shapenoise / (varg+shapenoise);
    double k0 = f*varg;
    double k1 = f*f;

#if 0
    // No Weight
    double w = 1.;
    double gdwdg = 0.;
#endif

#if 0
    // Simple inverse variance weight
    double w = 1./varg;
    double gdwdg = 0.;
#endif

#if 0
    // Weight that makes a shear average act like a distortion average
    // e = tanh(2*atanh(g))
    //   = 2*g/(1+g^2)
    // w = e/g = 2/(1+g^2)
    // dw/dg = -4g/(1+g^2)^2
    // g dw/dg = -g^2*w^2
    double w = 2./(1.+gsq);
    double gdwdg = -gsq*w*w;
#endif

#if 0
    // Better weight if lensing shear is small 
    double w = 1./sqrt(gsq + 2.25 * varg);
    //if (varg > 1.) w = 0.;
    double gdwdg = -w*w*w*gsq;
#endif

#if 1
    // De-weight large g
    int k = 3;
    double w = 1./std::pow(1.+gsq,k);
    double gdwdg = -2.*k*gsq*w/(1.+gsq);
#endif

    dr = w + 0.5*gdwdg*(1.-k0-k1*gsq);
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
        long id,flag;
        double x,y,sky,noise,ra,dec,g1,g2,nu,c00,c01,c11;
        strin >> 
            // id  x  y  sky  noise  flag  ra  dec  g1  g2  nu c00..
            id >> x >> y >> sky >> noise >> 
            flag >> ra >> dec >> g1 >> g2 >> nu >> c00 >> c01 >> c11;

        if (flag & ~ok_flags) continue;
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
    if (!(argc == 3 || argc == 5)) {
        std::cerr<<"Usage: shearave inFile outFile\n";
        std::cerr<<"or     shearave inFile outFile fileId answerslist\n";
        std::cerr<<"The first version merely converts the standard ouput\n";
        std::cout<<"from measureshears into numbers that can be used for\n";
        std::cout<<"computing averages.\n";
        std::cout<<"The second version also computes the averages and writes\n";
        std::cout<<"them to a file (appended) along with a fileId string.\n";
        exit(1);
    }

    //
    // Open files
    //

    std::string inFile = argv[1];
    std::string outFile = argv[2];
    std::ifstream fin(inFile.c_str());
    if (!fin) {
        std::cerr<<"Unable to open input file "<<inFile<<std::endl;
        exit(1);
    }
    std::ofstream fout(outFile.c_str());
    if (!fout) {
        std::cerr<<"Unable to open output file "<<outFile<<std::endl;
        exit(1);
    }
    std::string fileId;
    std::string answersFile;
    if (argc == 5) {
        fileId = argv[3];
        answersFile = argv[4];
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
    double shapenoise=0., sumw=0.;

    for(int i=0;i<nGal;++i) {
        std::complex<double> g = gVec[i];
        double varg = vargVec[i];
        double gsq = std::norm(g);
        double dr;
        // TODO: The 0.1 in the next line should be a parameter.
        if (varg > 0.1) continue; // Throw out sigma_g > 0.3
        double galw = Weight(g,varg,0.,dr);
        shapenoise += (gsq - varg) * galw;
        sumw += galw;
        //std::cout<<i<<"  "<<g<<"  "<<gsq<<"  "<<esq<<"  "<<vare<<"  "<<galw<<"  "<<shapenoise<<"  "<<sumw<<std::endl;
    }
    shapenoise /= 2.*sumw; // 2 so shapenoise per component
    std::cout<<"sigma_SN = "<<shapenoise<<std::endl;

    //
    // Calculate responsivity
    //

    double resp=0.,sumw2g2=0.;
    sumw = 0.;
    std::complex<double> sumwg = 0.;
    for (int i=0; i<nGal; ++i) {
        std::complex<double> g = gVec[i];
        double gsq = std::norm(g);
        double varg = vargVec[i];
        if (varg > 0.1) continue;
        double dr;
        double galw = Weight(g,varg,shapenoise,dr);
        resp += dr;
        sumw += galw;
        sumwg += galw*g;
        sumw2g2 += galw*galw*gsq;
    }
    resp /= sumw;
    std::cout<<"resp = "<<resp<<", ngal = "<<nGal<<std::endl;

    // 
    // Calculate estimate of applied shear
    //

    std::complex<double> gamma = sumwg / sumw;
    gamma /= resp;
    double varg = sumw2g2 / (2.*sumw*sumw);
    double sigg = sqrt(varg);
    double siggamma = tanh(atanh(sigg)/2.);
    std::cout<<"gamma = "<<gamma<<" +- "<<siggamma<<std::endl;

    //
    // Write out responsivity-corrected catalog
    //

    for (int i=0; i<nGal; i++) {
        fout.precision(12);
        fout<<
            xVec[i]<<"  "<<yVec[i]<<"  "<<
            raVec[i]<<"  "<<decVec[i]<<"  ";
        fout.precision(6);
        fout<<
            real(gVec[i])/resp<<"  "<<imag(gVec[i])/resp<<"  "<<
            sqrt(vargVec[i])/resp<<std::endl;
    }
    fout.close();

    if (argc == 5) {
        std::ofstream answersout(answersFile.c_str(),std::ios_base::app);
        fout.precision(6);
        answersout << fileId<<"  "<<real(gamma)<<"  "<<imag(gamma)<<"  ";
        answersout << siggamma<<"  "<<siggamma<<std::endl;
    }

    return 0;
}

