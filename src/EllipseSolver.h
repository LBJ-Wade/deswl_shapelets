#ifndef ELLIPSESOLVER_H
#define ELLIPSESOLVER_H

#include "Pixel.h"
#include "BVec.h"
#include <vector>
#include "NLSolver.h"

class BaseEllipseSolver : public NLSolver
{
  public :

    BaseEllipseSolver() {}
    virtual void UseNumericJ() = 0;
    virtual const BVec& GetB() const = 0;
    virtual void GetBCov(tmv::Matrix<double>& bcov) const = 0;
};

class ESImpl;

class EllipseSolver : public BaseEllipseSolver
{
  public :

    EllipseSolver(const std::vector<std::vector<Pixel> >& _pix, 
	int _order, double _sigma,
	bool _fixcen=false, bool _fixgam=false, bool _fixmu=false,
	bool _useflux=false);
    EllipseSolver(const std::vector<std::vector<Pixel> >& _pix, 
	const std::vector<BVec>& _psf, double _fp,
	int _order, double _sigma,
	bool _fixcen=false, bool _fixgam=false, bool _fixmu=false,
	bool _useflux=false);
    ~EllipseSolver();

    void F(const tmv::Vector<double>& x, tmv::Vector<double>& f) const;
    void J(const tmv::Vector<double>& x, const tmv::Vector<double>& f,
	tmv::Matrix<double>& df) const;

    void UseNumericJ();
    const BVec& GetB() const;
    void GetBCov(tmv::Matrix<double>& bcov) const;

    bool Solve(tmv::Vector<double>& x, tmv::Vector<double>& f,
	tmv::Matrix<double>* cov=0) const;
    bool TestJ(const tmv::Vector<double>& x, tmv::Vector<double>& f,
	std::ostream* os=0, double relerr=0.) const;

  private :

    ESImpl* pimpl;
};

class ESImpl2;

// Use the integration method, rather than least-squares, to find b.
class EllipseSolver2 : public BaseEllipseSolver
{
  public :

    EllipseSolver2(const std::vector<std::vector<Pixel> >& _pix,
	int _order, double _sigma, double _pixscale, 
	bool _fixcen=false, bool _fixgam=false, bool _fixmu=false,
	bool _useflux=false);
    EllipseSolver2(const std::vector<std::vector<Pixel> >& _pix,
	const std::vector<BVec>& _psf, double _fp,
	int _order, double _sigma, double _pixscale, 
	bool _fixcen=false, bool _fixgam=false, bool _fixmu=false,
	bool _useflux=false);
    ~EllipseSolver2();

    void F(const tmv::Vector<double>& x, tmv::Vector<double>& f) const;

    void J(const tmv::Vector<double>& x, const tmv::Vector<double>& f,
	tmv::Matrix<double>& df) const;

    bool Solve(tmv::Vector<double>& x, tmv::Vector<double>& f,
	tmv::Matrix<double>* invcov=0) const;
    bool TestJ(const tmv::Vector<double>& x, tmv::Vector<double>& f,
	std::ostream* os=0, double relerr=0.) const;

    void UseNumericJ();
    const BVec& GetB() const;
    void GetBCov(tmv::Matrix<double>& bcov) const;

  private :

    ESImpl2* pimpl;
};

#endif
