#include "Legendre2D.h"
#include "Function2D.cpp"
#include <algorithm>

template <class T>
void Legendre2D<T>::SetFunction(int _xorder, int _yorder, 
	const tmv::Vector<T>& fvect)
{
  if (xorder != _xorder || yorder != _yorder) {
    xorder = _xorder; yorder = _yorder;
    coeffs.reset(new tmv::Matrix<T>(xorder+1,yorder+1,0.));
  }
  int k=0;
  for(int m=0; m <= std::max(xorder,yorder); m++)
    for(int i=std::min(m,xorder);m-i<=std::min(m,yorder);i--) { 
      (*coeffs)(i,m-i) = fvect(k++);
    }
  Assert(k==(int)fvect.size());
}

template <class T>
Legendre2D<T>::Legendre2D(std::istream& fin) : Function2D<T>()
// Order of parameters is same as for Polynomial2D.  Difference
// is that instead of x,x^2,x^3,etc., we use P1(x),P2(x),P3(x),etc.
{
  fin >> xorder >> yorder >> bounds;
  if (!fin) throw std::runtime_error("reading order, bounds");
  coeffs.reset(new tmv::Matrix<T>(xorder+1,yorder+1,0.));
  int maxorder = std::max(xorder,yorder);
  for(int m = 0; m <= maxorder; m++) {
    for(int i=std::min(m,xorder);m-i<=std::min(m,yorder);i--) {
      fin >> (*coeffs)(i,m-i);
    }
  }
  if (!fin) throw std::runtime_error("reading (legendre) function");
}

template <class T>
void Legendre2D<T>::Write(std::ostream& fout) const
{
  int oldprec = fout.precision(6);
  std::ios::fmtflags oldf = fout.setf(std::ios::scientific,std::ios::floatfield);
  int maxorder = std::max(xorder,yorder);
  if (maxorder == 0) fout << "C " << (*coeffs)(0,0) << std::endl;
  else {
    fout << "L " << xorder << ' ' << yorder << ' ' << bounds << ' ';
    for(int m = 0; m <= maxorder; m++) {
      for(int i=std::min(m,xorder);m-i<=std::min(m,yorder);i--) {
	fout << (*coeffs)(i,m-i) << ' ';
      }
    }
    fout << std::endl;
  }
  if (!fout) throw std::runtime_error("writing (legendre) function");
  fout.precision(oldprec);
  fout.flags(oldf);
}

template <class T>
void Legendre2D<T>::AddLinear(T a, T b, T c)
{
  double xave = (GetXMin() + GetXMax())/2.;
  double yave = (GetYMin() + GetYMax())/2.;
  double xrange = GetXMax() - GetXMin();
  double yrange = GetYMax() - GetYMin();

  (*coeffs)(0,0) += a + b*xave + c*yave;
  (*coeffs)(1,0) += b*xrange/2.;
  (*coeffs)(0,1) += c*yrange/2.;
}

template <class T>
void Legendre2D<T>::LinearPreTransform(T , T , T , T , T , T )
{
  // Not implemented yet.
  Assert(false);
}

template <class T>
void Legendre2D<T>::operator+=(const Function2D<T>& rhs)
{
  const Legendre2D<T>* lrhs = dynamic_cast<const Legendre2D<T>*>(&rhs);
  Assert(lrhs);
  Assert(GetBounds() == lrhs->GetBounds());
  if (xorder == lrhs->xorder && yorder == lrhs->yorder) {
    *coeffs += *lrhs->coeffs;
  } else {
    int newxorder = std::max(xorder,lrhs->xorder);
    int newyorder = std::max(yorder,lrhs->yorder);
    std::auto_ptr<tmv::Matrix<T> > newc(new tmv::Matrix<T>(newxorder+1,newyorder+1,0.));
    newc->SubMatrix(0,xorder+1,0,yorder+1) = *coeffs;
    newc->SubMatrix(0,lrhs->xorder+1,0,lrhs->yorder+1) += *lrhs->coeffs;
    coeffs = newc;
    xorder = newxorder;
    yorder = newyorder;
  }
}

template <class T>
std::auto_ptr<Function2D<T> > Legendre2D<T>::DFDX() const 
// dP_2n(x)/dx = Sum_k=0..n-1 (4k+3) P_2k+1(x)
// dP_2n+1(x)/dx = Sum_k=0..n (4k+1) P_2k(x)
{
  if (xorder == 0) 
    return std::auto_ptr<Function2D<T> >(new Constant2D<T>());
  if (xorder == 1 && yorder == 0) 
    return std::auto_ptr<Function2D<T> >(new Constant2D<T>(
	  (*coeffs)(1,0)*2./(GetXMax()-GetXMin())));

  int newxorder = xorder-1;
  int newyorder = xorder > yorder ? yorder : yorder-1;

  std::auto_ptr<Legendre2D<T> > temp(
    new Legendre2D<T>(newxorder,newyorder,bounds));
  // initialized to 0's

  int maxorder = std::max(xorder,yorder);
  for(int i=xorder;i>=1;i--) 
    for(int j=std::min(maxorder-i,yorder);j>=0;j--) {
      if (i%2 == 0) for(int k=0;k<=i/2-1;k++) 
        (*temp->coeffs)(2*k+1,j) += (4.*k+3.)*(*coeffs)(i,j);
      else for(int k=0;k<=(i-1)/2;k++) 
        (*temp->coeffs)(2*k,j) += (4.*k+1.)*(*coeffs)(i,j);
    }
  maxorder = std::max(newxorder,newyorder);
  for(int i=newxorder;i>=0;i--) for(int j=std::min(maxorder-i,newyorder);j>=0;j--) {
    (*temp->coeffs)(i,j) *= 2./(GetXMax()-GetXMin());
  }
  return std::auto_ptr<Function2D<T> >(temp);
}

template <class T>
std::auto_ptr<Function2D<T> > Legendre2D<T>::DFDY() const 
// dP_2n(x)/dx = Sum_k=0..n-1 (4k+3) P_2k+1(x)
// dP_2n+1(x)/dx = Sum_k=0..n-1 (4k+1) P_2k(x)
{
  if (yorder == 0) 
    return std::auto_ptr<Function2D<T> >(new Constant2D<T>());
  if (yorder == 1 && xorder == 0) 
    return std::auto_ptr<Function2D<T> >(new Constant2D<T>(
	  (*coeffs)(0,1)*2./(GetYMax()-GetYMin())));

  int newxorder = yorder > xorder ? xorder : xorder-1;
  int newyorder = yorder-1;

  std::auto_ptr<Legendre2D<T> > temp(
    new Legendre2D<T>(newxorder,newyorder,bounds));
  // initialized to 0's

  int maxorder = std::max(xorder,yorder);
  for(int j=yorder;j>=1;j--) 
    for(int i=std::min(maxorder-j,xorder);i>=0;i--) {
      if (j%2 == 0) for(int k=0;k<=(j-2)/2;k++) 
        (*temp->coeffs)(i,2*k+1) += (4.*k+3.)*(*coeffs)(i,j);
      else for(int k=0;k<=(j-1)/2;k++) 
        (*temp->coeffs)(i,2*k) += (4.*k+1.)*(*coeffs)(i,j);
    }
  maxorder = std::max(newxorder,newyorder);
  for(int j=newyorder;j>=0;j--) for(int i=std::min(maxorder-j,newxorder);i>=0;i--) {
    (*temp->coeffs)(i,j) *= 2./(GetYMax()-GetYMin());
  }
  return std::auto_ptr<Function2D<T> >(temp);
}

template class Legendre2D<std::complex<double> >;
template class Legendre2D<double>;
