#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>

inline double SQR(double x) { return x*x; }

// Official code from Sarah:
// tot=0;
// n=0;
// for i=1:nset
//  dg1 = my_g1(i) - true_g1(i);
//  dg2 = my_g2(i) - true_g2(i);
//  tot = tot + dg1^2 + dg2^2;
//  n = n + 2;
// end
// meandgsq=tot/n;
// Q = 1e-4 / meandgsq;


int main(int argc, char* argv[])
{
    if (argc != 3) {
        std::cerr<<"Usage: calcq correctfile myanswers\n";
        exit(1);
    }
    std::ifstream correct(argv[1]);
    std::ifstream myanswers(argv[2]);

    std::string junk;
    getline(correct,junk);

    double e1_c, e2_c, e1_m, e2_m;
    int n=0;
    double sumy1 = 0., sumy2 = 0., sumq = 0.;
    double sumx1 = 0., sumx2 = 0., sumxsq1 = 0., sumxsq2 = 0.;
    double sumxy1 = 0., sumxy2 = 0.;
    std::string id;
    // y = mx + c
    // chisq = Sum (yi - m xi - c)^2
    // d/dm = 2 Sum (yi - m xi - c) (-xi) = 0
    // d/dc = 2 Sum (yi - m xi - c) = 0
    // Sum xi yi = m Sum xi^2 + c Sum xi
    // Sum yi = m Sum xi + c N
    // c = Sum yi / N - m Sum xi / N
    // Sum xi yi / N = m Sum xi^2 / N + Sum xi Sum yi / N^2 - m (Sum xi)^2/N^2
    // m = (N Sum xi yi - Sum xi Sum yi) / (N Sum xi^2 - (Sum xi)^2)
    while (correct >> id >> e1_c >> e2_c) {
        myanswers >> id >> e1_m >> e2_m >> junk >> junk;
        std::cout<<id<<": "<<e1_c<<" "<<e2_c<<" --- "<<e1_m<<" "<<e2_m;
        sumy1 += e1_m;
        sumy2 += e2_m;
        sumx1 += e1_c;
        sumx2 += e2_c;
        sumxsq1 += e1_c*e1_c;
        sumxsq2 += e2_c*e2_c;
        sumxy1 += e1_m*e1_c;
        sumxy2 += e2_m*e2_c;
        sumq += SQR(e1_m-e1_c) + SQR(e2_m-e2_c);
        n++;
        std::cout<<"   "<<SQR(e1_m-e1_c)/1.e-4<<"  "<<SQR(e2_m-e2_c)/1.e-4;
        std::cout<<std::endl;
    }
    sumq /= 2*n;
    std::cout<<"Q = "<<(1.e-4 / sumq)<<std::endl;
    double altq = 1.e-4 / (SQR((sumy1-sumx1)/n) + SQR((sumy2-sumx2)/n));
    std::cout<<"Alt Q = "<<altq<<std::endl;
    double m1 = (sumxy1*n - sumx1 * sumy1) / (sumxsq1*n - sumx1 * sumx1)-1.;
    double m2 = (sumxy2*n - sumx2 * sumy2) / (sumxsq2*n - sumx2 * sumx2)-1.;
    double c1 = sumy1 / n - m1 * sumx1/n;
    double c2 = sumy2 / n - m2 * sumx2/n;
    //std::cout<<"avex1 = "<<sumx1<<" / "<<n<<" = "<<sumx1/n<<std::endl;
    //std::cout<<"avey1 = "<<sumy1<<" / "<<n<<" = "<<sumy1/n<<std::endl;
    //std::cout<<"varx1 = "<<sumxsq1<<" / "<<n<<" - "<<sumx1/n<<"^2 = "<<sumxsq1/n-(sumx1/n)*(sumx1/n)<<std::endl;
    //std::cout<<"m1 = ("<<sumxy1<<" - "<<sumx1<<" * "<<sumy1<<") / ("<<sumxsq1<<" - "<<sumx1<<"^2) - 1\n";
    std::cout<<"e1: m = "<<m1<<", c = "<<c1<<std::endl;
    //std::cout<<"avex2 = "<<sumx2<<" / "<<n<<" = "<<sumx2/n<<std::endl;
    //std::cout<<"avey2 = "<<sumy2<<" / "<<n<<" = "<<sumy2/n<<std::endl;
    //std::cout<<"varx2 = "<<sumxsq2<<" / "<<n<<" - "<<sumx2/n<<"^2 = "<<sumxsq2/n-(sumx2/n)*(sumx2/n)<<std::endl;
    //std::cout<<"m2 = ("<<sumxy2<<" - "<<sumx2<<" * "<<sumy2<<") / ("<<sumxsq2<<" - "<<sumx2<<"^2) - 1\n";
    std::cout<<"e2: m = "<<m2<<", c = "<<c2<<std::endl;
}
