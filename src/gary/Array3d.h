// 	$Id: Array3d.h,v 1.8 2007-03-14 15:37:33 reiko3 Exp $	
// Cubic 3-dimensional array of doubles.  Adapted from Jarvis
//---------------------------------------------------------------------------
#ifndef Array3dH
#define Array3dH

#include "Matrix.h"

namespace mv {

//---------------------------------------------------------------------------

template <class T> class Cube {

public:
  Cube() : m(0),va() {}
  Cube(size_t mm, T val=0) 
    : m(mm),va(val,m*m*m) {}
  Cube(const Cube& rhs)
    : m(rhs.m),va(rhs.va) {}
  Cube(size_t mm, const valarray<T>& vv)
    : m(mm),va(vv)
    { MVAssert(va.size() == m*m*m); }
  Cube<T>& operator=(const Cube<T>& rhs)
    { if (&rhs == this) return *this;
      if(va.size()==0) {m=rhs.m;va.resize(m*m*m);}
      MVAssert(rhs.m == m);
      va = rhs.va; 
      return *this; }
  Cube<T>& operator=(const T rhs) {va=rhs; return *this;}
  ~Cube() {}
  void resize(size_t mm) 
    { m = mm; va.resize(m*m*m);}

  size_t GetM() const {return m;}
  const valarray<T>& GetValArray() const {return va;}

  slice slice1(size_t j, size_t k) const {return slice(j*m+k,m,m*m); }
  slice slice2(size_t i, size_t k) const {return slice(i*m*m+k,m,m); }
  slice slice3(size_t i, size_t j) const {return slice((i*m+j)*m,m,1); }
  size_t index(size_t i,size_t j, size_t k) const { return (i*m+j)*m+k; }


  T& operator()(size_t i,size_t j, size_t k) 
    { MVAssert(i<m && j<m && k<m);
      return va[index(i,j,k)]; }
  T operator()(size_t i,size_t j, size_t k) const
    { MVAssert(i<m && j<m && k<m);
      return va[index(i,j,k)]; }
  T Get(size_t i,size_t j, size_t k) const
    { return (*this)(i,j,k); }
  Slice_ref<T> GetSlice1(size_t j, size_t k) const
    { MVAssert(j<m && k<m);
      return Slice_ref<T>(va,slice1(j,k)); }
  Slice_ref<T> GetSlice2(size_t j, size_t k) const
    { MVAssert(j<m && k<m);
      return Slice_ref<T>(va,slice2(j,k)); }
  Slice_ref<T> GetSlice3(size_t j, size_t k) const
    { MVAssert(j<m && k<m);
      return Slice_ref<T>(va,slice3(j,k)); }

  void Set(size_t i,size_t j, size_t k, T t)
    { MVAssert(i<m && j<m && k<m);
      va[index(i,j,k)] = t; }
  void Zero()
    { va = 0; }

  Cube& operator*=(T x)
    { va *= x; 
      return *this; }
  Cube& operator/=(T x)
    { va /= x; 
      return *this; }

 private:
  size_t m;
  valarray<T> va;

}; // Cube

} // namespace mv

#endif
