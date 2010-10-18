// 	$Id: Laguerre.h,v 1.40 2007-05-11 20:34:16 garyb Exp $	
//---------------------------------------------------------------------------
#ifndef LAGUERRE_H
#define LAGUERRE_H
//---------------------------------------------------------------------------
// Manipulation of Laguerre-decomposition-vector representations of images.

#include <string>
#include <iostream>
#include <sstream>

#include "Std.h"
#include "Matrix.h"
#include "Shear.h"
#include <sstream>

namespace laguerre {

    using mv::DMatrix;
    using mv::CMatrix;
    using mv::SqDMatrix;
    using mv::SqCMatrix;
    using mv::DVector;
    using mv::CVector;
    using mv::Vector;
    using mv::Matrix;
    using mv::SqMatrix;

    class LaguerreError: public MyException 
    {
    public: 
        LaguerreError(const string &m=""): MyException("Laguerre Error: " + m) {}
    };
    class LaguerreInsufficientOrder: public LaguerreError
    {
    public: 
        LaguerreInsufficientOrder(const string &m=""): 
            LaguerreError("Requested order is beyond available order" + m) {}
    };
    class LaguerreNonConvergent: public LaguerreError
    {
    public: 
        LaguerreNonConvergent(const string &m=""): 
            LaguerreError(
                "Failure to converge to centroid/size/roundness soln " + m) 
        {}
    };

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    // First define an index taking p & q of the eigenfunction, and
    // derived vector and matrix classes which can be indexed this way.

    class PQIndex 
    {
        // Index vectors/matrices of Laguerre coefficients using the
        // p & q quantum numbers (or others).
        // When we want to store the coefficients as a complex array,
        //  the storage order will be increasing in N=p+q, with decreasing
        //  m=p-q within a given N:
        //  pq:      00;  10, 01; 20, 11, 02; 30, 21, 12, 03; ...
        //  cIndex:   0    1   2   3   4  5   6   7    8   9  ...
        //
        // Also will need to get the index into a DVector that stores
        // the real degrees of freedom of the (Hermitian) Laguerre vector
        // for a real-valued image.  The storage order here is:
        //  pq:        Re(00); Re(10), Im(10); Re(20), Im(20), Re(11); ...
        //  rIndex:     0        1       2       3       4       5     ...
        // 
        // Methods include increments for p, q, N, and m.  This object can 
        // be used as an index into any of our Laguerre arrays.
    private:
        int p;
        int q;
        // Formulae for mapping into linear arrays:
        static int makeCIndex(const int p_, const int q_) {
            return (p_+q_)*(p_+q_+1)/2+q_;
        }
        static int makeRIndex(const int p_, const int q_) {
            return (p_+q_)*(p_+q_+1)/2 + 2*MIN(p_,q_);
        }

    public:
        PQIndex(): p(0), q(0) {} 
        PQIndex(const int p_, const int q_) {setPQ(p_,q_);}

        bool pqValid() const {return p>=0 && q>=0;}

        int getP() const {return p;}
        int getQ() const {return q;}
        PQIndex& setPQ(const int p_=0, const int q_=0) 
        {
            p=p_; q=q_; 
            return *this;
        }
        int N() const {return p+q;}
        int m() const {return p-q;}
        PQIndex& setNm(const int N, const int m) 
        {
            //Assert(abs(m)<=N && (N-m)%2==0);
            p=(N+m)/2; q=(N-m)/2; 
            return *this;
        }

        int cIndex() const {return makeCIndex(p,q);}
        int rIndex() const {return makeRIndex(p,q);}
        bool needsConjugation() const {return p<q;}
        bool isReal() const {return p==q;}
        int iSign() const {if (p<q) return -1; else return (p==q? 0 : 1);}

        // Operations that update the indices:
        void incP() {p++;}
        void incQ() {q++;}
        void incN() {p++; q++;}  //raise N by 2, same m
        void decN() {p--; q--;}  //lower N by 2, same m (could be invalid)
        void incm() {p++; q--;} // raise m by 2, same N (could be invalid)
        void decm() {p--; q++;} // lower m by 2, same N (could be invalid)
        // get next one in complex sequence
        PQIndex& operator++()  
        { 
            if (p==0) {p=q+1; q=0;}
            else {--p; ++q;}
            return *this;
        }

        // get next pq index that has m>=0
        PQIndex& nextDistinct()  
        { 
            if (p-q<2) {p=p+q+1; q=0;}
            else {--p; ++q;}
            return *this;
        }


        // Functions to report incremented/decremented indices without
        // updating this index:
        PQIndex swapPQ() const {return PQIndex(q,p);} 
        PQIndex pp1() const {return PQIndex(p+1,q);}
        PQIndex qp1() const {return PQIndex(p,q+1);}
        PQIndex pm1() const {return PQIndex(p-1,q);}
        PQIndex qm1() const {return PQIndex(p,q-1);}

        bool operator==(const PQIndex rhs) const {return p==rhs.p && q==rhs.q;}

        // Other useful things:
        static int size(int order) 
        {
            // Size of a CVector to this order N, same as number of real DOF:
            //Assert(order>=0);
            return (order+1)*(order+2)/2;
        }

        // write and ??? read
        void write(ostream& os) const;

        // Returns true if index has advanced past order:
        bool pastOrder(const int order) const {return p+q>order;}
    };

    inline ostream& operator<<(ostream& os, const laguerre::PQIndex& pq) 
    { pq.write(os); return os; }

    //--------------------------------------------------------------
    // Next object is a vector of Laguerre coefficients.  Note this is
    // a HANDLE to the coefficient vector, so it can be passed into
    // subroutines without referencing.  Copy/assignment create a new link; 
    // for fresh copy, use duplicate() method.
    //
    // LVectors are assumed to be Hermitian complex (b_qp=b_pq*), and the
    // internal storage currently enforces this and actually stores the
    // data as a DVector of the real degrees of freedom.
    // So when you change b_qp, you are also changing b_pq.

    class LVector 
    {
    private:
        int *order;
        DVector *v;
        mutable int	*pcount;

        static DComplex getElement(DVector* v, PQIndex(pq)) 
        {
            int i=pq.rIndex();
            if (pq.isReal()) 
                return DComplex( (*v)[i]);
            else if (pq.needsConjugation())
                return DComplex( (*v)[i], -(*v)[i+1]);
            else 
                return DComplex( (*v)[i], (*v)[i+1]);
        }

        static void setElement(DVector* v, PQIndex(pq), DComplex z) 
        {
            int i=pq.rIndex();
            (*v)[i] = z.real();
            if (pq.isReal()) {
                if (z.imag()!=0.) {
                    std::ostringstream oss;
                    oss << "Imaginary part " << z.imag()
                        << " in assignment to real LVector element";
                    throw LaguerreError(oss.str());
                }
                return;
            }
            i++;
            if (pq.needsConjugation())
                (*v)[i] = -z.imag();
            else 
                (*v)[i] = z.imag();
        }
        static DVector dotVector;
        static int dotVectorOrder;
        static void checkDotVector(int order);
    public:
        LVector(int ord_=0) :  
            order(new int(ord_)), pcount(new int(1)),
            v(new DVector(PQIndex::size(*order),0.)) {}
        LVector(const LVector& rhs) : 
            v(rhs.v), pcount(rhs.pcount), 
            order(rhs.order) {(*pcount)++;}
        LVector(const DVector& rhs, int order_);
        LVector(const CVector& cv);
        LVector& operator=(const LVector& rhs) 
        {
            if (v==rhs.v) return *this;
            if (--(*pcount)==0) {delete v; delete pcount; delete order;}
            v = rhs.v; pcount = rhs.pcount; order=rhs.order;
            (*pcount)++;
            return *this;
        }
        ~LVector() 
        {
            if (--(*pcount)==0) {
                delete v;
                delete pcount;
                delete order;
            }
        }
        LVector duplicate() const {
            LVector fresh(*order); 
            *(fresh.v) = *v;
            return fresh;
        }

        int getOrder() const {return *order;}
        int rSize() const {return v->size();}
        void clear() {*v = 0.;}
        void resize(int neworder) {
            *order = neworder;
            v->resize(PQIndex::size(*order));
            *v = 0.;
        }

        // Access the real-representation vector directly.
        DVector& rVector() {return *v;}
        const DVector& rVector() const {return *v;}
        // Return the real vector that can be RHS of a dot product:
        DVector realRHS() const;
        // Build a CVector representation and return ptr to it:
        CVector* cVector() const;
        CVector* cVectorConjugate() const;
        CVector* kSpaceCVector() const;
        CVector* kSpaceCVectorConjugate() const;

        // Element access.  Index checking if Assert() is active.
        // Access is either by a PQIndex or by (p,q) integer pairs.

        // Set elements (no op= available!)
        DComplex set(const PQIndex pq, DComplex z) 
        {
            //Assert(pq.pqValid() && !pq.pastOrder(*order));
            setElement(v,pq,z);
            return z;
        }
        DComplex set(int p_, int q_, DComplex z) 
        {
            return set(PQIndex(p_,q_),z);
        }

        // Now read-only element access:
        DComplex operator[](const PQIndex& pq) const 
        {
            //Assert(pq.pqValid() && !pq.pastOrder(*order));
            return getElement(v,pq);
        }
        DComplex operator()(int p_, int q_) const 
        {
            return operator[](PQIndex(p_,q_));
        }

        LVector& operator*=(double s) {*v *= s; return *this;}
        LVector  operator*(double s) const 
        {
            LVector fresh(*order); 
            *(fresh.v) = *v * s;
            return fresh;
        }
        LVector& operator/=(double s) {*v /= s; return *this;}
        LVector  operator/(double s) const 
        {
            LVector fresh(*order); 
            *(fresh.v) = *v / s;
            return fresh;
        }
        LVector& operator+=(const LVector& rhs) 
        {
            //Assert(*order== *(rhs.order));
            *v += *(rhs.v); 
            return *this;
        }
        LVector& operator-=(const LVector& rhs) 
        {
            //Assert(*order== *(rhs.order));
            *v -= *(rhs.v); 
            return *this;
        }

        // write and ??? read
        void write(ostream& os, int maxorder=-1) const;
        friend ostream& operator<<(ostream& os, const LVector& lv);
        void read(istream& is);
        friend istream& operator>>(istream& is, LVector& lv);

        // Rotate represented object by theta:
        // (note this is faster than using RotationLTransform)
        void rotate(double theta); 

        // Get the total flux or flux within an aperture of size R*sigma
        // Use all monopole terms unless maximum is specified by maxP.
        double flux(int maxP=-1) const;
        double apertureFlux(double R, int maxP=-1) const;

        // Things that use Laguerre functions:
        void fillLaguerre(const double x) ;	// Fill with L_q^m(x)
        // fill with psi_pq(z)
        // First version assumes x/y already mapped to unit circle,
        // but psi needs to be normalized by 1/sigma^2
        void fillPsi(double xunit, double yunit, double sigma=1.);
        // Second version assumes x/y are in pixel coords, and the
        // psi is to have a basis given by Ellipse, so coord map
        // and psi normalization are taken from Ellipse:
        void fillPsi(double xraw, double yraw, const Ellipse& e);
        double dot(const LVector& rhs) const;

        // fill with FTransformed Psi (xunit = x*sigma, etc)
        void fillKPsi(double xunit, double yunit);


        // This one returns the wave functions as a real vector to
        // be dotted with the real-valued b vector to produce I(x)
        static DVector realPsi(
            int N, double xx, double yy, const Ellipse& e);

        // Same thing, but fill a matrix with the psi values at a
        // vector of points:
        static void realPsi(
            int N, const DVector& xx, const DVector& yy, 
            const Ellipse& e, DMatrix* psi);
        // Return reference to a matrix that generates realPsi transformations
        // under infinitesimal point transforms (translate, dilate, shear).
        // Returned matrix is at least as large as needed to go order x (order+2)
        static const DMatrix& Generator(int iparam, int order);
        // The choices for above:
        static const  int iX = 0;
        static const  int iY = 1;
        static const  int iMu = 2;
        static const  int iE1 = 3;
        static const  int iE2 = 4;
    };

    ostream& operator<<(ostream& os, const LVector& lv);
    istream& operator>>(istream& is, LVector& lv);

    // This function finds the innermost radius at which the integrated flux
    // of the LVector's shape crosses the specified threshold, using the first
    // maxP monopole terms (or all, if maxP omitted)
    extern double fluxRadius(const LVector& lv, double threshold, int maxP=-1);

    //--------------------------------------------------------------
    //
    // Next class is a covariance matrix for Laguerre vector.  Internal 
    // storage is as a real covariance matrix for the real/imag parts of 
    // the vector.  Interface gives you the (complex) covariances of chosen
    // pqIndex pairs.

    // Again this is a HANDLE to the coefficient vector, so it can be passed into
    // subroutines without referencing.  Copy/assignment create a new link; 
    // for fresh copy, use duplicate() method.
    class LCovar 
    {
    private:
        int *order;
        mutable int	*pcount;
        SqDMatrix *m;
    public:
        LCovar(int ord_=0): order(new int(ord_)), pcount(new int(1)),
        m(new SqDMatrix(PQIndex::size(*order),0.)) {}
        LCovar(const LCovar& rhs): m(rhs.m), pcount(rhs.pcount), 
        order(rhs.order) {(*pcount)++;}
        // Build an LCovar from a SqDMatrix for the real degrees of freedom.
        // SqDMatrix must have same dimension as needed for LCovar of chosen order
        // **Note we are not checking for symmetry of input matrix.
        LCovar(const SqDMatrix& rhs, int order_);
        LCovar& operator=(const LCovar& rhs) 
        {
            if (m==rhs.m) return *this;
            if (--(*pcount)==0) {delete m; delete pcount; delete order;}
            m = rhs.m; pcount = rhs.pcount; order=rhs.order;
            (*pcount)++;
            return *this;
        }
        ~LCovar() 
        {
            if (--(*pcount)==0) {
                delete m;
                delete pcount;
                delete order;
            }
        }
        LCovar duplicate() const 
        {
            LCovar fresh(*order); 
            *(fresh.m) = *m;
            return fresh;
        }

        int getOrder() const {return *order;}
        int rSize() const {return m->size();}
        void clear() {*m = 0.;}
        void resize(int neworder) 
        {
            *order = neworder;
            m->resize(PQIndex::size(*order));
            *m = 0.;
        }

        // Access the real-representation vector directly.
        SqDMatrix& rMatrix() {return *m;}
        const SqDMatrix& rMatrix() const {return *m;}

        // Element read
        DComplex operator()(PQIndex pq1, PQIndex pq2) const;
        DComplex operator()(int p1, int q1, int p2, int q2) const 
        { return operator()(PQIndex(p1,q1),PQIndex(p2,q2)); }

        // Element write.  Note that it is necessary to give two complex
        // simultaneously to allow writing the real version of the matrix:
        void set(
            PQIndex pq1, PQIndex pq2,
            DComplex covpq1pq2, DComplex covpq1qp2);
    };

    //--------------------------------------------------------------
    //
    // Next class is a transformation matrix for Laguerre vector.  Internal 
    // storage is as a matrix over the real degrees of freedom.
    // Interface gives you the (complex) matrix elements of  pqIndex pairs.

    // Again this is a HANDLE, so it can be passed into
    // subroutines without referencing.  Copy/assignment create a new link; 
    // for fresh copy, use duplicate() method.
    class LTransform 
    {
    private:
        int *orderIn, *orderOut;
        mutable int	*pcount;
        DMatrix *m;
    public:
        LTransform(int orderOut_=0, int orderIn_=0): 
            orderIn(new int(orderIn_)), orderOut(new int(orderOut_)), 
            pcount(new int(1)),
            m(new DMatrix(PQIndex::size(*orderOut),PQIndex::size(*orderIn),0.))
        {}
        LTransform(const LTransform& rhs): 
            m(rhs.m), pcount(rhs.pcount), 
            orderIn(rhs.orderIn), orderOut(rhs.orderOut) {(*pcount)++;}
        // Build an LTransform from a DMatrix for the real degrees of freedom.
        // Matrix must have correct dimensions.
        LTransform(const DMatrix& rhs, int orderOut_, int orderIn_);
        LTransform& operator=(const LTransform& rhs) 
        {
            if (m==rhs.m) return *this;
            if (--(*pcount)==0) {
                delete m;
                delete pcount;
                delete orderIn;
                delete orderOut;
            }
            m = rhs.m; pcount = rhs.pcount; 
            orderIn=rhs.orderIn; orderOut=rhs.orderOut;
            (*pcount)++;
            return *this;
        }
        ~LTransform() 
        {
            if (--(*pcount)==0) {
                delete m;
                delete pcount;
                delete orderIn;
                delete orderOut;
            }
        }
        LTransform duplicate() const 
        {
            LTransform fresh(*orderOut, *orderIn); 
            *(fresh.m) = *m;
            return fresh;
        }

        int getOrderIn() const {return *orderIn;}
        int getOrderOut() const {return *orderOut;}
        int rSizeIn() const {return m->getN();}
        int rSizeOut() const {return m->getM();}
        void clear() {*m = 0.;}
        void identity(); // Set to identity transformation
        void resize(int newOrderOut, int newOrderIn) 
        {
            *orderIn = newOrderIn;
            *orderOut = newOrderOut;
            m->resize(PQIndex::size(*orderOut),PQIndex::size(*orderIn));
            *m = 0.;
        }

        // Access the real-representation vector directly.
        DMatrix& rMatrix() {return *m;}
        const DMatrix& rMatrix() const {return *m;}

        // Element read
        DComplex operator()(PQIndex pq1, PQIndex pq2) const;
        DComplex operator()(int p1, int q1, int p2, int q2) const 
        {
            return operator()(PQIndex(p1,q1),PQIndex(p2,q2));
        }

        // Element write.  Note that it is necessary to give two complex
        // simultaneously to allow writing the real version of the matrix:
        void set(
            PQIndex pq1, PQIndex pq2,
            DComplex Cpq1pq2, DComplex Cqp1pq2);

        // Operate on other Laguerre vectors/matrices
        LVector operator*(const LVector rhs) const;
        LCovar operator*(const LCovar rhs) const;
        LTransform operator*(const LTransform rhs) const;
        LTransform& operator*=(const LTransform rhs);

        // Precede this tranformation by a rotation (faster than
        // using RotationLTransform)
        void preRotation(double theta);
    };

    // Here are the primary types of transformations:
    // For the point transforms, set coordShift=false if we want
    // to transform the FLUX on a fixed coordinate grid.  Set true
    // if want to describe the same flux on a transformed COORD system.

    // Shear:
    LTransform MakeLTransform(
        Shear eta, int orderOut, int orderIn, bool coordShift=false);
    // Dilation:
    LTransform MakeLTransform(
        double mu, int orderOut, int orderIn, bool coordShift=false);

    // Translation:
    LTransform MakeLTransform(
        Position<double> x0, int orderOut, int orderIn, bool coordShift=false);

    // Rotation:
    LTransform RotationLTransform(
        double theta, int orderOut, int orderIn, bool coordShift=false);

    // Combination of above 3 (specify intermediate precision if 
    //     the default of MAX(out,in) is not wanted):
    LTransform MakeLTransform(
        const Ellipse& e, int orderOut, int orderIn,
        bool coordShift=false, int orderIntermediate=-1);

    // Convolution with PSF:
    LTransform MakeLTransform(
        const LVector psf, const double D,
        const int orderOut, const int orderIn, const int orderStar);

}  // namespace laguerre

#endif
