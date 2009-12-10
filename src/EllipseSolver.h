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
    virtual void useNumericJ() = 0;
    virtual const BVec& getB() const = 0;
    virtual void callF(
        const tmv::Vector<double>& x, tmv::Vector<double>& f) const = 0;
    virtual void getBCov(tmv::Matrix<double>& bcov) const = 0;
};

class EllipseSolver : public BaseEllipseSolver
{
public :

    EllipseSolver(
        const std::vector<PixelList>& pix, 
        int order, double sigma, 
        bool fixcen=false, bool fixgam=false, bool fixmu=false,
        bool useflux=false);
    EllipseSolver(
        const std::vector<PixelList>& pix, 
        const std::vector<BVec>& psf, double fp,
        int order, double sigma, 
        bool fixcen=false, bool fixgam=false, bool fixmu=false,
        bool useflux=false);
    ~EllipseSolver();

    void calculateF(const tmv::Vector<double>& x, tmv::Vector<double>& f) const;
    void calculateJ(
        const tmv::Vector<double>& x, const tmv::Vector<double>& f,
        tmv::Matrix<double>& df) const;

    void useNumericJ();
    const BVec& getB() const;
    void getBCov(tmv::Matrix<double>& bcov) const;
    void getCovariance(tmv::Matrix<double>& cov) const;
    void getInverseCovariance(tmv::Matrix<double>& invcov) const;

    // CallF takes x and f of length 5, rather than whatever shorter
    // length that F takex (depending on if things are fixed).
    void callF(const tmv::Vector<double>& x, tmv::Vector<double>& f) const;
    bool solve(tmv::Vector<double>& x, tmv::Vector<double>& f) const;
    bool testJ(
        const tmv::Vector<double>& x, tmv::Vector<double>& f,
        std::ostream* os=0, double relerr=0.) const;
    void calculateNumericH(
        const tmv::Vector<double>& x, const tmv::Vector<double>& f,
        tmv::SymMatrix<double>& h) const;

private :

    struct ESImpl;

    ESImpl* _pimpl;
};

// Use the integration method, rather than least-squares, to find b.
class EllipseSolver2 : public BaseEllipseSolver
{
public :

    EllipseSolver2(
        const std::vector<PixelList>& pix,
        int order, double sigma, double pixscale,
        bool fixcen=false, bool fixgam=false, bool fixmu=false,
        bool useflux=false);
    EllipseSolver2(
        const std::vector<PixelList>& pix,
        const std::vector<BVec>& psf, double fp,
        int order, double sigma, double pixscale,
        bool fixcen=false, bool fixgam=false, bool fixmu=false,
        bool useflux=false);
    ~EllipseSolver2();

    void calculateF(
        const tmv::Vector<double>& x, tmv::Vector<double>& f) const;
    void calculateJ(
        const tmv::Vector<double>& x, const tmv::Vector<double>& f,
        tmv::Matrix<double>& df) const;

    void callF(const tmv::Vector<double>& x, tmv::Vector<double>& f) const;
    bool solve(tmv::Vector<double>& x, tmv::Vector<double>& f) const;
    bool testJ(
        const tmv::Vector<double>& x, tmv::Vector<double>& f,
        std::ostream* os=0, double relerr=0.) const;

    void useNumericJ();
    const BVec& getB() const;
    void getBCov(tmv::Matrix<double>& bcov) const;

private :

    struct ESImpl2;

    ESImpl2* _pimpl;
};

#endif
