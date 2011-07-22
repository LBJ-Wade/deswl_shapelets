
#include "NLSolver.h"
#include <cstdlib>
#include <iostream>

#define DOLINEAR
#define DONONLINEAR
#define DOPOWELL
#define DOROSENBROCK

#define SETMETHOD(solver) \
    do { \
        solver.useDogleg(); \
        solver.setOutput(std::cout); \
        solver.useVerboseOutput(); \
    } while (false)


#define USEBEST

#ifdef USE_TMV

const int N=2;
const int M=5;

class Linear : public NLSolver 
{
public :

    Linear() : b(M,0.), A(M,N,0.)
    {
        tmv::Vector<double> x(N,1.);
        for(int i=0;i<std::min(M,N);++i) {
            A(i,i) = 0.5*i - 6.;
            if (i>0 && 2*i<M) A(2*i,i) = 9.+i;
            if (i>0 && 2*i<M) A(2*i-1,i) = 1.-2*i;
            if (i>0 && 2*i<N) A(i,2*i) = -3.-i;
            if (i>0 && 2*i<N) A(i,2*i-1) = -12.+3*i;
        }
        A.col(0).addToAll(2.);
        A.row(0).addToAll(3.);
        b = A*x;
        //b.addToAll(100.);
    }

    void calculateF(const tmv::Vector<double>& x, tmv::Vector<double>& f) const
    { f = A*x-b; }

    void calculateJ(const tmv::Vector<double>& , const tmv::Vector<double>& ,
           tmv::Matrix<double>& j) const
    { j = A; }

    void calculateH(
        const tmv::Vector<double>& , const tmv::Vector<double>& , 
        const tmv::Matrix<double>& j, tmv::SymMatrix<double>& h) const
    { h = j.adjoint() * j; }

protected:

    tmv::Vector<double> b;
    tmv::Matrix<double> A;
};

const int MM=50;
class NonLinear : public NLSolver 
{
public : 

    NonLinear(double sigma) : b(MM), t(MM), _sigma(sigma)
    {
        tmv::Vector<double> x(4);
        x << 0.5, 3.5, 1.3, 1.7;

        for(int i=0;i<MM;++i) {
            t(i) = (2.0*(i+0.5)-1.0*MM)/MM;
            b(i) = func(t(i),x) + _sigma*(double(rand())/RAND_MAX -0.5);
        }
    }

    double func(double t, const tmv::Vector<double>& x) const
    {
        return x(2)*exp(x(0)*t) + x(3)*sin(x(1)*t);
    }

    void calculateF(const tmv::Vector<double>& x, tmv::Vector<double>& f) const
    {
        for(int i=0;i<MM;++i) f(i) = (func(t(i),x)-b(i));
        f /= _sigma;
    }

    void calculateJ(const tmv::Vector<double>& x, const tmv::Vector<double>& ,
           tmv::Matrix<double>& j) const
    {
        for(int i=0;i<MM;++i) {
            j(i,0) = x(2)*t(i)*exp(x(0)*t(i));
            j(i,1) = x(3)*t(i)*cos(x(1)*t(i));
            j(i,2) = exp(x(0)*t(i));
            j(i,3) = sin(x(1)*t(i));
        }
        j /= _sigma;
    }

protected:

    tmv::Vector<double> b;
    tmv::Vector<double> t;
    double _sigma;
};

class NonLinear2 : public NLSolver 
{
public : 

    NonLinear2(double sigma) : b(MM), t(MM), A(MM,2), y(2), _sigma(sigma)
    {
        double xx[2] = {0.5,3.5};
        double yy[2] = {1.3,1.7};

        for(int i=0;i<MM;++i) {
            t(i) = (2.0*(i+0.5)-1.0*MM)/MM;
            b(i) = yy[0]*exp(xx[0]*t(i)) + yy[1]*sin(xx[1]*t(i)) +
                _sigma*(2.*double(rand())/RAND_MAX - 1.);
        }
        A.saveDiv();
    }

    double func(double t, const tmv::Vector<double>& x, 
                const tmv::Vector<double>& y) const
    {
        return y(0)*exp(x(0)*t) + y(1)*sin(x(1)*t);
    }

    void calculateF(const tmv::Vector<double>& x, tmv::Vector<double>& f) const
    {
        for(int i=0;i<MM;++i) {
            A(i,0) = exp(x(0)*t(i));
            A(i,1) = sin(x(1)*t(i));
        }
        A.resetDiv();
        y = b/A;
        f = A*y-b;
        f /= _sigma;
    }

    void calculateJ(const tmv::Vector<double>& x, const tmv::Vector<double>& f,
           tmv::Matrix<double>& j) const
    {
        // At A y = At b
        // df = A dy + dA y
        // (A+dA)t (A+dA) (y+dy) = (A+dA)t b
        // dAt A y + At dA y + At A dy = dAt b
        // At A dy = -dAt (Ay-b) - At dA y
        // At A = Rt R where A = QR is the QR decomposition of A
        // dy = R^-1 Rt^-1 (-dAt f - At dA y)
        // df = A dy + dA y
        //    = A R^-1 Rt^-1 (-dAt f - At dA y) + dA y
        tmv::ConstUpperTriMatrixView<double> R = A.qrd().getR();

        tmv::Vector<double> dA = DiagMatrixViewOf(t)*A.col(0);
        tmv::Vector<double> dy = -A.transpose()*dA*y(0);
        dy(0) -= dA*f*_sigma;
        dy /= R.transpose();
        dy /= R;
        j.col(0) = A*dy + dA*y(0);

        for(int i=0;i<MM;++i) {
            dA(i) = t(i)*cos(x(1)*t(i));
        }
        dy = -A.transpose()*dA*y(1);
        dy(1) -= dA*f*_sigma;
        dy /= R.transpose();
        dy /= R;
        j.col(1) = A*dy + dA*y(1);

        j /= _sigma;
    }

    const tmv::Vector<double>& getY() const { return y; }

protected:

    tmv::Vector<double> b;
    tmv::Vector<double> t;
    mutable tmv::Matrix<double> A;
    mutable tmv::Vector<double> y;
    double _sigma;
};

class Powell : public NLSolver 
{
public :

    Powell() {}

    void calculateF(const tmv::Vector<double>& x, tmv::Vector<double>& f) const
    {
        f(0) = x(0);
        f(1) = 10.*x(0)/(0.1+x(0)) + 2.0*std::pow(x(1),2);
    }
    void calculateJ(const tmv::Vector<double>& x, const tmv::Vector<double>& ,
           tmv::Matrix<double>& j) const
    {
        j(0,0) = 1.0;
        j(0,1) = 0.0;
        j(1,0) = std::pow(0.1+x(0),-2);
        j(1,1) = 4.0*x(1);
    }
    void calculateH(const tmv::Vector<double>& x, const tmv::Vector<double>& f,
           const tmv::Matrix<double>& j, tmv::SymMatrix<double>& h) const
    {
        // H = JT J + Sum_k f(k) d^2f/(dxi dxj)
        h = j.transpose() * j;
        h(0,0) += -2.0*std::pow(0.1+x(0),-3)*f(0);
        h(1,1) += 4.0*f(1);
    }

};

class Rosen : public NLSolver 
{
public :

    Rosen(double _lam) : lambda(_lam) {}

    void calculateF(const tmv::Vector<double>& x, tmv::Vector<double>& f) const
    {
        f(0) = 10.*(x(1) - x(0)*x(0));
        f(1) = 1.-x(0);
        f(2) = lambda;
    }
    void calculateJ(const tmv::Vector<double>& x, const tmv::Vector<double>& ,
           tmv::Matrix<double>& j) const
    {
        j(0,0) = -20.*x(0);
        j(1,0) = -1.;
        j(2,0) = 0.;
        j(0,1) = 10.;
        j(1,1) = 0.;
        j(2,1) = 0.;
    }

    void calculateH(const tmv::Vector<double>& , const tmv::Vector<double>& f,
           const tmv::Matrix<double>& j, tmv::SymMatrix<double>& h) const
    {
        // H = JT J + Sum_k f(k) d^2f/(dxi dxj)
        h = j.transpose() * j;
        h(0,0) += -20.*f(0);
    }



private :

    double lambda;
};

int main() try 
{
#ifdef DOLINEAR
    {
        std::cout<<"Test Linear:\n";
        Linear lin;
        SETMETHOD(lin);
#ifdef USEBEST
        lin.useNewton();
#endif
        tmv::Vector<double> x(N,0.);
        tmv::Vector<double> f(M);
        tmv::Matrix<double> cov1(N,N);
        tmv::SymMatrix<double> h(N);

        //std::cout<<"xinit = "<<x<<std::endl;
        //if (!(lin.testJ(x,f,&std::cout))) return 1;
        if (lin.solve(x,f))
            std::cout<<"Success: \n";
        else 
            std::cout<<"Failure: \n";
        lin.getCovariance(cov1);
        std::cout<<"xfinal = "<<x<<std::endl;
        //std::cout<<"f = "<<f<<std::endl;
        std::cout<<"Norm(f) = "<<Norm(f)<<std::endl;

        // This one usually finds an exact f == 0, since the problem is 
        // so easy.
        if (!(Norm(f) <= 1.e-10)) {
            std::cout<<"Linear Test failed.\n";
            return 1;
        }

        // Check the covariance matrix:
        std::cout<<"Covariance matrix = \n"<<cov1<<std::endl;
        tmv::Matrix<double> j(M,N);
        lin.calculateJ(x,f,j);
        lin.calculateH(x,f,j,h);
        std::cout<<"exact h = "<<h<<std::endl;
        lin.calculateNumericH(x,f,h);
        std::cout<<"numeric h = "<<h<<std::endl;
        tmv::SymMatrix<double> cov2 = h.inverse();
        if (!(Norm(cov1-cov2) <= 1.e-6*(1.+Norm(cov1)))) {
            std::cout<<"cov1 = "<<cov1<<std::endl;
            std::cout<<"cov2 = "<<cov2<<std::endl;
            std::cout<<"Covariance matrices don't match\n";
            return 1;
        }

    }
#endif

#ifdef DONONLINEAR
    {
        std::cout<<"Test NonLinear:\n";
        double sigma = 1.e-3;
        NonLinear nlin(sigma);
        SETMETHOD(nlin);
#ifdef USEBEST
        nlin.useLM();
#endif
        nlin.setMinStep(1.e-15);
        tmv::Vector<double> x(4);
        x << 1., 1., 1., 1.;
        tmv::Vector<double> f(MM);
        tmv::Matrix<double> cov1(4,4);
        tmv::SymMatrix<double> h(4);

        //if (!(nlin.testJ(x,f,&std::cout))) return 1;
        if (nlin.solve(x,f))
            std::cout<<"Success: \n";
        else 
            std::cout<<"Failure: \n";
        nlin.getCovariance(cov1);
        std::cout<<"xfinal = "<<x<<std::endl;
        //if (!(nlin.testJ(x,f,&std::cout))) return 1;
        // This one is looking for a local minimum, not an exact answer.
        // The minimum should have Norm(f) ~= 2*sigma
        if (!(Norm(f) <= 3.0)) {
            std::cout<<"f = "<<f<<std::endl;
            std::cout<<"Norm(f) = "<<Norm(f)<<std::endl;
            std::cout<<"NonLinear Test failed.\n";
            return 1;
        }

        // Check the inverse covariance matrix:
        std::cout<<"Covariance matrix = \n"<<cov1<<std::endl;
        nlin.calculateNumericH(x,f,h);
        tmv::SymMatrix<double> cov2 = h.inverse();
        if (!(Norm(cov1-cov2) <= 1.e-6*(1.+Norm(cov1)))) {
            std::cout<<"cov1 = "<<cov1<<std::endl;
            std::cout<<"cov2 = "<<cov2<<std::endl;
            std::cout<<"Covariance matrices don't match\n";
            return 1;
        }

    }
    {
        std::cout<<"Test NonLinear2:\n";
        double sigma = 1.e-3;
        NonLinear2 nlin(sigma);
        SETMETHOD(nlin);
#ifdef USEBEST
        nlin.useLM();
#endif
        nlin.setMinStep(1.e-15);
        tmv::Vector<double> x(2);
        x << 1., 1.;
        tmv::Vector<double> f(MM);
        tmv::Matrix<double> cov1(2,2);
        tmv::SymMatrix<double> h(2);

        //if (!(nlin.testJ(x,f,&std::cout))) return 1;
        if (nlin.solve(x,f))
            std::cout<<"Success: \n";
        else 
            std::cout<<"Failure: \n";
        nlin.getCovariance(cov1);
        //if (!(nlin.testJ(x,f,&std::cout))) return 1;
        tmv::Vector<double> xtot(4);
        xtot.subVector(0,2) = x;
        xtot.subVector(2,4) = nlin.getY();
        std::cout<<"xfinal = "<<xtot<<std::endl;
        // This one is looking for a local minimum, not an exact answer.
        // The minimum should have Norm(f) ~= 4.
        if (!(Norm(f) <= 5.0)) {
            std::cout<<"f = "<<f<<std::endl;
            std::cout<<"Norm(f) = "<<Norm(f)<<std::endl;
            std::cout<<"NonLinear2 Test failed.\n";
            return 1;
        }

        // Check the inverse covariance matrix:
        std::cout<<"Covariance matrix = \n"<<cov1<<std::endl;
        nlin.calculateNumericH(x,f,h);
        tmv::SymMatrix<double> cov2 = h.inverse();
        if (!(Norm(cov1-cov2) <= 1.e-6*(1.+Norm(cov1)))) {
            std::cout<<"cov1 = "<<cov1<<std::endl;
            std::cout<<"cov2 = "<<cov2<<std::endl;
            std::cout<<"Covariance matrices don't match\n";
            return 1;
        }

    }
#endif

#ifdef DOPOWELL
    {
        std::cout<<"Test Powell\n";
        Powell pow;
        pow.setTol(1.e-20,1.e-15);
        pow.setMinStep(1.e-15);
        pow.setTau(1.0);
        SETMETHOD(pow);
#ifdef USEBEST
        pow.useDogleg();
#endif
        tmv::Vector<double> x(2); x(0) = 3.; x(1) = 1.;
        tmv::Vector<double> f(2);
        tmv::Matrix<double> cov1(2,2);
        tmv::SymMatrix<double> h(2);

        //std::cout<<"xinit = "<<x<<std::endl;
        if (pow.solve(x,f))
            std::cout<<"Success: \n";
        else 
            std::cout<<"Failure: \n";
        std::cout<<"xfinal = "<<x<<std::endl;
        //std::cout<<"f = "<<f<<std::endl;
        std::cout<<"Norm(f) = "<<Norm(f)<<std::endl;
        // This test is hard for most solvers, but Dogleg should get it.
        // It will probably come back as finding the local minimum instead 
        // of the exact answer, but the two are equivalent in this case.
        if (!(Norm(f) <= 1.e-15)) {
            std::cout<<"Powell Test failed.\n";
            return 1;
        }

        // Check the inverse covariance matrix:
        // The hessian becomes singular at the solution, so the regular
        // solution covariance does not produce a reliable covariance matrix.
        // We need to switch on svd to get a better answer.
        // We start with the correct answer, so the solution is actually
        // trivial.  We just need the better covariance.
        pow.useSVD();
        pow.getCovariance(cov1);
        std::cout<<"Covariance matrix = \n"<<cov1<<std::endl;
        tmv::Matrix<double> j(2,2);
        pow.calculateJ(x,f,j);
        pow.calculateH(x,f,j,h);
        h.divideUsing(tmv::SV);
        tmv::SymMatrix<double> cov2 = h.inverse();
        if (!(Norm(cov1-cov2) <= 1.e-6*(1.+Norm(cov1)))) {
            std::cout<<"cov1 = "<<cov1<<std::endl;
            std::cout<<"cov2 = "<<cov2<<std::endl;
            std::cout<<"Covariance matrices don't match\n";
            return 1;
        }

    }
#endif

#ifdef DOROSENBROCK
    {
        std::cout<<"Test Rosenbrock\n";
        const int ntest = 4;
        double lambda[ntest] = {0., 1.e-3, 1., 1.e4};
        for(int i=0;i<ntest;++i) {
            Rosen ros(lambda[i]);
            ros.setTol(1.e-20,1.e-10);
            ros.setMinStep(1.e-14);
            ros.setTau(1.e-3);
            SETMETHOD(ros);
#ifdef USEBEST
            ros.useHybrid();
#endif
#ifdef ROSENBROCKDIRECTH
            ros.useDirectH();
#endif

            tmv::Vector<double> x(2); x(0) = -1.2; x(1) = 1.;
            tmv::Vector<double> x0(2); x0(0) = 1.; x0(1) = 1.;
            tmv::Vector<double> f(3);
            tmv::Matrix<double> cov1(2,2);
            tmv::SymMatrix<double> h(2);

            //std::cout<<"xinit = "<<x<<std::endl;
            if (ros.solve(x,f))
                std::cout<<"Success: \n";
            else 
                std::cout<<"Failure: \n";
            ros.getCovariance(cov1);
            std::cout<<"xfinal = "<<x<<std::endl;
            //std::cout<<"Norm(xf-x0) = "<<Norm(x0-x)<<std::endl;
            //std::cout<<"f = "<<f<<std::endl;
            std::cout<<"Norm(f) = "<<Norm(f)<<std::endl;
            // Lambda > 0 is a challenge for most solvers.
            // Hybrid should get it even when lambda is very large.
            if (!(Norm(x-x0) <= 1.e-8)) {
                std::cout<<"Rosenbrock Test failed for lambda = "<<lambda<<".\n";
                return 1;
            }

            // Check the inverse covariance matrix:
            std::cout<<"Covariance matrix = \n"<<cov1<<std::endl;
            //ros.calculateNumericH(x,f,h);
            // The numeric H is too noisy to get an accurate inverse.
            // So use the direct calculation.
            tmv::Matrix<double> j(3,2);
            ros.calculateJ(x,f,j);
            ros.calculateH(x,f,j,h);
            tmv::SymMatrix<double> cov2 = h.inverse();
            if (!(Norm(cov1-cov2) <= 1.e-6*(1.+Norm(cov1)))) {
                std::cout<<"cov1 = "<<cov1<<std::endl;
                std::cout<<"cov2 = "<<cov2<<std::endl;
                std::cout<<"Covariance matrices don't match\n";
                return 1;
            }
        }
    }
#endif

    std::cout<<"\n\nNLSolver passed all tests.\n\n";
    return 0;
} 
catch (tmv::Error& e)
{
    std::cout<<"Caught "<<e<<std::endl;
    return 1;
}
catch (std::exception& e)
{
    std::cout<<"Caught "<<e.what()<<std::endl;
    return 1;
}
catch (...)
{
    std::cout<<"Caught unknown error\n";
    return 1;
}

#else

const int N=2;
const int M=5;

class Linear : public NLSolver 
{
public :

    Linear() : b(M), A(M,N)
    {
        b.setZero();
        A.setZero();
        DVector x(N);
        x.setConstant(1.);
        for(int i=0;i<std::min(M,N);++i) {
            A(i,i) = 0.5*i - 6.;
            if (i>0 && 2*i<M) A(2*i,i) = 9.+i;
            if (i>0 && 2*i<M) A(2*i-1,i) = 1.-2*i;
            if (i>0 && 2*i<N) A(i,2*i) = -3.-i;
            if (i>0 && 2*i<N) A(i,2*i-1) = -12.+3*i;
        }
        A.col(0).TMV_addToAll(2.);
        A.row(0).TMV_addToAll(3.);
        b = A*x;
        //b.TMV_addToAll(100.);
    }

    void calculateF(const DVector& x, DVector& f) const
    { f = A*x-b; }

    void calculateJ(const DVector& , const DVector& ,
           DMatrix& j) const
    { j = A; }

protected:

    DVector b;
    DMatrix A;
};

const int MM=50;
class NonLinear : public NLSolver 
{
public : 

    NonLinear(double sigma) : b(MM), t(MM), _sigma(sigma)
    {
        DVector x(4);
        x << 0.5, 3.5, 1.3, 1.7;

        for(int i=0;i<MM;++i) {
            t(i) = (2.0*(i+0.5)-1.0*MM)/MM;
            b(i) = func(t(i),x) + _sigma*(double(rand())/RAND_MAX -0.5);
        }
    }

    double func(double t, const DVector& x) const
    {
        return x(2)*exp(x(0)*t) + x(3)*sin(x(1)*t);
    }

    void calculateF(const DVector& x, DVector& f) const
    {
        for(int i=0;i<MM;++i) f(i) = (func(t(i),x)-b(i));
        f /= _sigma;
    }

    void calculateJ(const DVector& x, const DVector& ,
           DMatrix& j) const
    {
        for(int i=0;i<MM;++i) {
            j(i,0) = x(2)*t(i)*exp(x(0)*t(i));
            j(i,1) = x(3)*t(i)*cos(x(1)*t(i));
            j(i,2) = exp(x(0)*t(i));
            j(i,3) = sin(x(1)*t(i));
        }
        j /= _sigma;
    }

protected:

    DVector b;
    DVector t;
    double _sigma;
};

class NonLinear2 : public NLSolver 
{
public : 

    NonLinear2(double sigma) : b(MM), t(MM), A(MM,2), y(2), _sigma(sigma)
    {
        double xx[2] = {0.5,3.5};
        double yy[2] = {1.3,1.7};

        for(int i=0;i<MM;++i) {
            t(i) = (2.0*(i+0.5)-1.0*MM)/MM;
            b(i) = yy[0]*exp(xx[0]*t(i)) + yy[1]*sin(xx[1]*t(i)) +
                _sigma*(2.*double(rand())/RAND_MAX - 1.);
        }
    }

    double func(double t, const DVector& x, 
                const DVector& y) const
    {
        return y(0)*exp(x(0)*t) + y(1)*sin(x(1)*t);
    }

    void calculateF(const DVector& x, DVector& f) const
    {
        for(int i=0;i<MM;++i) {
            A(i,0) = exp(x(0)*t(i));
            A(i,1) = sin(x(1)*t(i));
        }
        A.qr().solve(b,&y);
        f = A*y-b;
        f /= _sigma;
    }

    void calculateJ(const DVector& x, const DVector& f,
           DMatrix& j) const
    {
        // At A y = At b
        // df = A dy + dA y
        // (A+dA)t (A+dA) (y+dy) = (A+dA)t b
        // dAt A y + At dA y + At A dy = dAt b
        // At A dy = -dAt (Ay-b) - At dA y
        // At A = Rt R where A = QR is the QR decomposition of A
        // dy = R^-1 Rt^-1 (-dAt f - At dA y)
        // df = A dy + dA y
        //    = A R^-1 Rt^-1 (-dAt f - At dA y) + dA y
        
        Eigen::QR<DMatrix> QR_Solver = A.qr();

        DVector dA = TMV_DiagMatrixViewOf(t)*A.col(0);
        DVector dy = -A.transpose()*dA*y(0);
        dy(0) -= (dA.transpose()*f)(0,0)*_sigma;
        //dy /= R.transpose();
        //dy /= R;
        QR_Solver.matrixR().transpose().solveTriangularInPlace(dy);
        QR_Solver.matrixR().solveTriangularInPlace(dy);
        j.col(0) = A*dy + dA*y(0);

        for(int i=0;i<MM;++i) {
            dA(i) = t(i)*cos(x(1)*t(i));
        }
        dy = -A.transpose()*dA*y(1);
        dy(1) -= (dA.transpose()*f)(0,0)*_sigma;
        //dy /= R.transpose();
        //dy /= R;
        QR_Solver.matrixR().transpose().solveTriangularInPlace(dy);
        QR_Solver.matrixR().solveTriangularInPlace(dy);
        j.col(1) = A*dy + dA*y(1);

        j /= _sigma;
    }

    const DVector& getY() const { return y; }

protected:

    DVector b;
    DVector t;
    mutable DMatrix A;
    mutable DVector y;
    double _sigma;
};

class Powell : public NLSolver 
{
public :

    Powell() {}

    void calculateF(const DVector& x, DVector& f) const
    {
        f(0) = x(0);
        f(1) = 10.*x(0)/(0.1+x(0)) + 2.0*std::pow(x(1),2);
    }
    void calculateJ(const DVector& x, const DVector& ,
           DMatrix& j) const
    {
        j(0,0) = 1.0;
        j(0,1) = 0.0;
        j(1,0) = std::pow(0.1+x(0),-2);
        j(1,1) = 4.0*x(1);
    }

};

class Rosen : public NLSolver 
{
public :

    Rosen(double _lam) : lambda(_lam) {}

    void calculateF(const DVector& x, DVector& f) const
    {
        f(0) = 10.*(x(1) - x(0)*x(0));
        f(1) = 1.-x(0);
        f(2) = lambda;
    }
    void calculateJ(const DVector& x, const DVector& ,
           DMatrix& j) const
    {
        j(0,0) = -20.*x(0);
        j(1,0) = -1.;
        j(2,0) = 0.;
        j(0,1) = 10.;
        j(1,1) = 0.;
        j(2,1) = 0.;
    }

private :

    double lambda;
};

int main() try 
{
#ifdef DOLINEAR
    {
        std::cout<<"Test Linear:\n";
        Linear lin;
        SETMETHOD(lin);
#ifdef USEBEST
        lin.useDogleg();
#endif
        
        DVector x(N);
        x.setZero();
        DVector f(M);

        //std::cout<<"xinit = "<<x<<std::endl;
        //if (!(lin.testJ(x,f,&std::cout))) return 1;
        if (lin.solve(x,f))
            std::cout<<"Success: \n";
        else 
            std::cout<<"Failure: \n";
        std::cout<<"xfinal = "<<x<<std::endl;
        //std::cout<<"f = "<<f<<std::endl;
        std::cout<<"Norm(f) = "<<f.norm()<<std::endl;

        // This one usually finds an exact f == 0, since the problem is 
        // so easy.
        if (!(f.norm() <= 1.e-10)) {
            std::cout<<"Linear Test failed.\n";
            return 1;
        }
    }
#endif

#ifdef DONONLINEAR
    {
        std::cout<<"Test NonLinear:\n";
        double sigma = 1.e-3;
        NonLinear nlin(sigma);
        SETMETHOD(nlin);
#ifdef USEBEST
        nlin.useHybrid();
#endif

        nlin.setMinStep(1.e-15);
        DVector x(4);
        x << 1., 1., 1., 1.;
        DVector f(MM);

        //if (!(nlin.testJ(x,f,&std::cout))) return 1;
        if (nlin.solve(x,f))
            std::cout<<"Success: \n";
        else 
            std::cout<<"Failure: \n";
        std::cout<<"xfinal = "<<x<<std::endl;
        //if (!(nlin.testJ(x,f,&std::cout))) return 1;
        // This one is looking for a local minimum, not an exact answer.
        // The minimum should have Norm(f) ~= 2*sigma
        if (!(f.norm() <= 3.0)) {
            std::cout<<"f = "<<f<<std::endl;
            std::cout<<"Norm(f) = "<<f.norm()<<std::endl;
            std::cout<<"NonLinear Test failed.\n";
            return 1;
        }
    }
    {
        std::cout<<"Test NonLinear2:\n";
        double sigma = 1.e-3;
        NonLinear2 nlin(sigma);
        SETMETHOD(nlin);
#ifdef USEBEST
        nlin.useHybrid();
#endif

        nlin.setMinStep(1.e-15);
        DVector x(2);
        x << 1., 1.;
        DVector f(MM);

        //if (!(nlin.testJ(x,f,&std::cout))) return 1;
        if (nlin.solve(x,f))
            std::cout<<"Success: \n";
        else 
            std::cout<<"Failure: \n";
        //if (!(nlin.testJ(x,f,&std::cout))) return 1;
        DVector xtot(4);
        xtot.TMV_subVector(0,2) = x;
        xtot.TMV_subVector(2,4) = nlin.getY();
        std::cout<<"xfinal = "<<xtot<<std::endl;
        // This one is looking for a local minimum, not an exact answer.
        // The minimum should have Norm(f) ~= 4.
        if (!(f.norm() <= 5.0)) {
            std::cout<<"f = "<<f<<std::endl;
            std::cout<<"Norm(f) = "<<f.norm()<<std::endl;
            std::cout<<"NonLinear2 Test failed.\n";
            return 1;
        }
    }
#endif

#ifdef DOPOWELL
    {
        std::cout<<"Test Powell\n";
        Powell pow;
        pow.setTol(1.e-20,1.e-15);
        pow.setMinStep(1.e-15);
        pow.setTau(1.0);
        SETMETHOD(pow);
#ifdef USEBEST
        pow.useDogleg();
#endif

        DVector x(2); x(0) = 3.; x(1) = 1.;
        DVector f(2);

        //std::cout<<"xinit = "<<x<<std::endl;
        if (pow.solve(x,f))
            std::cout<<"Success: \n";
        else 
            std::cout<<"Failure: \n";
        std::cout<<"xfinal = "<<x<<std::endl;
        //std::cout<<"f = "<<f<<std::endl;
        std::cout<<"Norm(f) = "<<f.norm()<<std::endl;
        // This test is hard for most solvers, but Dogleg should get it.
        // It will probably come back as finding the local minimum instead 
        // of the exact answer, but the two are equivalent in this case.
        if (!(f.norm() <= 1.e-15)) {
            std::cout<<"Powell Test failed.\n";
            return 1;
        }
    }
#endif

#ifdef DOROSENBROCK
    {
        std::cout<<"Test Rosenbrock\n";
        const int ntest = 4;
        double lambda[ntest] = {0., 1.e-3, 1., 1.e4};
        for(int i=0;i<ntest;++i) {
            Rosen ros(lambda[i]);
            ros.setTol(1.e-20,1.e-10);
            ros.setMinStep(1.e-14);
            ros.setTau(1.e-3);
            SETMETHOD(ros);
#ifdef USEBEST
            ros.useHybrid();
#endif

            DVector x(2); x(0) = -1.2; x(1) = 1.;
            DVector x0(2); x0(0) = 1.; x0(1) = 1.;
            DVector f(3);

            //std::cout<<"xinit = "<<x<<std::endl;
            if (ros.solve(x,f))
                std::cout<<"Success: \n";
            else 
                std::cout<<"Failure: \n";
            std::cout<<"xfinal = "<<x<<std::endl;
            //std::cout<<"Norm(xf-x0) = "<<Norm(x0-x)<<std::endl;
            //std::cout<<"f = "<<f<<std::endl;
            std::cout<<"Norm(f) = "<<f.norm()<<std::endl;
            // Lambda > 0 is a challenge for most solvers.
            // Hybrid should get it even when lambda is very large.
            if (!((x-x0).norm() <= 1.e-8)) {
                std::cout<<"Rosenbrock Test failed for lambda = "<<lambda<<".\n";
                return 1;
            }
        }
    }
#endif

    std::cout<<"\n\nNLSolver passed all tests.\n\n";
    return 0;
} 
catch (std::exception& e)
{
    std::cout<<"Caught "<<e.what()<<std::endl;
    return 1;
}
catch (...)
{
    std::cout<<"Caught unknown error\n";
    return 1;
}

#endif

