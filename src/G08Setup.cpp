#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>

int main(int argc, char* argv[])
{
  if (argc < 2) {
    std::cerr<<"Usage: setup LowNoise_Known_psf.txt\n";
    exit(1);
  }
  std::ifstream fin(argv[1]);
  std::string junk;
  getline(fin,junk); // ignore first line
  std::string catnum, psfnum;
  while (fin >> catnum >> psfnum) {
    // Copy centroid file:
    std::string catfile = std::string("set") + catnum + ".cat";
    std::ifstream cen_in("centroid.dat");
    std::ofstream cen_out(catfile.c_str());
    int id = 1;
    double x,y;
    while (cen_in >> x >> y) {
        cen_out <<id++<<"  "<<x<<"  "<<y<<std::endl;
    }
    cen_in.close();
    cen_out.close();

    // Copy fitpsf file:
    std::string psffile = std::string("set") + catnum + ".fitpsf";
    std::string inpsffile = std::string("../Stars/set") + psfnum + ".fitpsf";
    std::ifstream psf_in(inpsffile.c_str());
    std::ofstream psf_out(psffile.c_str());
    psf_out<<psf_in.rdbuf();
    psf_in.close();
    psf_out.close();
  }
}
