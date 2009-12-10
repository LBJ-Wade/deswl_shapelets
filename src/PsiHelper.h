#ifndef PSIHELPER_H
#define PSIHELPER_H

#include "TMV.h"

// Here is the order of p,q along the indices of psi:
//
// k 0 1,2 3,4 5 6,7 8,9 10,11 12,13 14 15,16 17,18 19,20 21,22 23,24 25,26 27
// p 0  1   2  1  3   2    4     3    2   5     4     3     6     5     4    3
// q 0  0   0  1  0   1    0     1    2   0     1     2     0     1     2    3
// n 0  1   2  2  3   3    4     4    4   5     5     5     6     6     6    6
// m 0  1   2  0  3   1    4     2    0   5     3     1     6     4     2    0

void makePsi1(
    const tmv::VectorView<double>& psi, 
    std::complex<double> z, int order);

void makePsi(
    const tmv::MatrixView<double>& psi, 
    const tmv::Vector<std::complex<double> >& z, int order,
    const tmv::DiagMatrix<double>* coeff=0);

void augmentPsi(
    tmv::Matrix<double>& psi,
    const tmv::Vector<std::complex<double> >& z, int order);

void setupGx(tmv::Matrix<double>& Gx, int order1, int order2);
void setupGy(tmv::Matrix<double>& Gy, int order1, int order2);
void setupGg1(tmv::Matrix<double>& Gg1, int order1, int order2);
void setupGg2(tmv::Matrix<double>& Gg2, int order1, int order2);
void setupGmu(tmv::Matrix<double>& Gmu, int order1, int order2);
void setupGth(tmv::Matrix<double>& Gth, int order1, int order2);

#endif
