// $Id: Solve.h,v 2.1 2003/05/06 21:36:47 garyb Exp $
// Template to find the zero of an equation
// Currently uses bisection method, no solution caching.

#ifndef SOLVE_H
#define SOLVE_H

#include "Std.h"
#include <cmath>

namespace solve {

    class SolveError: public laguerre::MyException 
    {
    public:
        SolveError(const std::string m) :
            laguerre::MyException("Solve error; "+m) {}
    };

    const double defaultTolerance=1.e-7;
    const int defaultMaxSteps=40;

    template <class F, class T=double>
    class Solve 
    {
    private:
        T	uBound;
        T	lBound;
        T	xTolerance;
        int	maxSteps;
        const F&  func;
    public:
        Solve(const F& func_, T lb_=0., T ub_=1.):
            func(func_), lBound(lb_), uBound(ub_), xTolerance(defaultTolerance),
            maxSteps(defaultMaxSteps) {}
        void setMaxSteps(int m) {maxSteps=m;}
        T    getXTolerance() const {return xTolerance;}
        void setXTolerance(T tol) {xTolerance=tol;}
        void setBounds(T lb, T ub) {lBound=lb; uBound=ub;}

        T    root() const 
        {
            T dx,f,fmid,xmid,rtb;

            f=func(lBound);
            fmid=func(uBound);
            if (f*fmid >= 0.0) 
                throw SolveError("Root is not bracketed");
            rtb = f < 0.0 ?
                (dx=uBound-lBound,lBound) :
                (dx=lBound-uBound,uBound);
            for (int j=1;j<=maxSteps;j++) {
                fmid=func(xmid=rtb+(dx *= 0.5));
                if (fmid <= 0.0) rtb=xmid;
                if ( (std::abs(dx) < xTolerance) || fmid == 0.0) return rtb;
            }
            throw SolveError("Too many bisections");
            return 0.0;
        }
    };

} // namespace solve

#endif
