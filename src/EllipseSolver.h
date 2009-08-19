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
    virtual ~BaseEllipseSolver() {}
    virtual void UseNumericJ() = 0;
    virtual const BVec& GetB() const = 0;
    virtual void CallF(const tmv::Vector<double>& x,
	tmv::Vector<double>& f) const = 0;
    virtual void GetBCov(tmv::Matrix<double>& bcov) const = 0;
};

class ESImpl;

class EllipseSolver : public BaseEllipseSolver
{
  public :

    EllipseSolver(const std::vector<PixelList>& pix, 
	int order, double sigma, bool desqa,
	bool fixcen=false, bool fixgam=false, bool fixmu=false,
	bool useflux=false);
    EllipseSolver(const std::vector<PixelList>& pix, 
	const std::vector<BVec>& psf, double fp,
	int order, double sigma, bool desqa,
	bool fixcen=false, bool fixgam=false, bool fixmu=false,
	bool useflux=false);
    ~EllipseSolver();

    void F(const tmv::Vector<double>& x, tmv::Vector<double>& f) const;
    void J(const tmv::Vector<double>& x, const tmv::Vector<double>& f,
	tmv::Matrix<double>& df) const;

    void UseNumericJ();
    const BVec& GetB() const;
    void GetBCov(tmv::Matrix<double>& bcov) const;

    // CallF takes x and f of length 5, rather than whatever shorter
    // length that F takex (depending on if things are fixed).
    void CallF(const tmv::Vector<double>& x, tmv::Vector<double>& f) const;
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

    EllipseSolver2(const std::vector<PixelList>& pix,
	int order, double sigma, double pixscale,
	bool fixcen=false, bool fixgam=false, bool fixmu=false,
	bool useflux=false);
    EllipseSolver2(const std::vector<PixelList>& pix,
	const std::vector<BVec>& psf, double fp,
	int order, double sigma, double pixscale,
	bool fixcen=false, bool fixgam=false, bool fixmu=false,
	bool useflux=false);
    ~EllipseSolver2();

    void F(const tmv::Vector<double>& x, tmv::Vector<double>& f) const;
    void J(const tmv::Vector<double>& x, const tmv::Vector<double>& f,
	tmv::Matrix<double>& df) const;

    void CallF(const tmv::Vector<double>& x, tmv::Vector<double>& f) const;
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
