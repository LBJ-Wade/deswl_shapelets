#ifndef FormH
#define FormH

// Formatting for simplified stream output
// see Stroustrup(2000), p. 635

#include <complex>
#include <sstream>
#include <iostream>
#include <string>
#include "dbg.h"

template <typename T> class BoundForm;

class Form 
{
public:
    Form() : 
        _prc(6), _wdt(0), _fmt(), _base(std::ios_base::dec),
        _just(std::ios_base::left), _new_fill_ch(0),
        _upper(0), _plus(0), _trail(0), _bool_alpha(0),
        _n_trail(1), _trail_ch(' ')  
    {}

    template <typename T> BoundForm<T> operator()(T val) const;

    Form& prec(int p) { _prc = p; return *this; }

    Form& sci() { _fmt = std::ios_base::scientific; return *this; }
    Form& fix() { _fmt = std::ios_base::fixed; return *this; }
    Form& gen() { _fmt = ~std::ios_base::floatfield; return *this; }

    Form& width(int w) { _wdt = w; return *this; }
    Form& fill(char c) { _new_fill_ch = c; return *this; }

    Form& dec() { _base = std::ios_base::dec; return *this; }
    Form& oct() { _base = std::ios_base::oct; return *this; }
    Form& hex() { _base = std::ios_base::hex; return *this; }

    Form& left() { _just = std::ios_base::left; return *this; }
    Form& right() { _just = std::ios_base::right; return *this; }
    Form& internal() { _just = std::ios_base::internal; return *this; }

    Form& uppercase(bool b=true) { _upper = b?1:-1; return *this; }
    Form& showpos(bool b=true) { _plus = b?1:-1; return *this; }
    Form& showpoint(bool b=true) { _trail = b?1:-1; return *this; }
    Form& boolalpha(bool b=true) { _bool_alpha = b?1:-1; return *this; }

    Form& trail(int n, char ch=' ') 
    { _n_trail = n; _trail_ch = ch; return *this; }

private:
    template <typename T> 
    friend std::ostream& operator<<(std::ostream&, const BoundForm<T>&);

    friend void SetupFloat(std::ostream&, const Form&);
    friend void SetupInt(std::ostream&, const Form&);

    int _prc; // precision
    int _wdt; // width, 0 means as wide as necessary
    std::ios_base::fmtflags _fmt; // general sci, or fixed
    std::ios_base::fmtflags _base; // dec, hex, oct
    std::ios_base::fmtflags _just; // left, right, internal fill
    char _new_fill_ch; // fill character
    int _upper; // +1 to have uppercase E,X, -1 turn off, (0 leave as is)
    int _plus; // +1 to have explicit plus for positive values, -1 off, 0 same
    int _trail; // +1 to write trailing zeros, -1 off, 0 same
    int _bool_alpha; // +1 to write "true","false", -1 off, 0 same
    int _n_trail; // number of spaces after output
    char _trail_ch; // character of trailing "spaces"

};

template <typename T>
class BoundForm 
{
public:
    const Form& f;
    T val;
    BoundForm(const Form& f_, T val_) : f(f_), val(val_) {}
};

template <typename T>
inline BoundForm<T> Form::operator()(T val) const 
{ return BoundForm<T>(*this,val); }

inline void SetupFloat(std::ostream& s, const Form& f)
{
    s.precision(f._prc);
    s.setf(f._fmt,std::ios_base::floatfield);
    s.setf(f._just,std::ios_base::adjustfield);
    if (f._wdt) s.width(f._wdt);
    if (f._new_fill_ch) s.fill(f._new_fill_ch);
    if (f._upper && f._fmt == std::ios_base::scientific) {
        if (f._upper>0) s.setf(std::ios_base::uppercase);
        else s.unsetf(std::ios_base::uppercase); 
    }
    if (f._plus) {
        if (f._plus>0) s.setf(std::ios_base::showpos); 
        else s.unsetf(std::ios_base::showpos); 
    }
    if (f._trail) {
        if (f._trail>0) s.setf(std::ios_base::showpoint); 
        else s.unsetf(std::ios_base::showpoint); 
    }
}

inline void SetupInt(std::ostream& s, const Form& f)
{
    s.setf(f._just,std::ios_base::adjustfield);
    s.setf(f._base,std::ios_base::basefield);
    if (f._wdt) s.width(f._wdt);
    if (f._new_fill_ch) s.fill(f._new_fill_ch);
    if (f._upper && f._base == std::ios_base::hex) {
        if (f._upper>0) s.setf(std::ios_base::uppercase); 
        else s.unsetf(std::ios_base::uppercase); 
    }
    if (f._plus) {
        if (f._plus>0) s.setf(std::ios_base::showpos); 
        else s.unsetf(std::ios_base::showpos); 
    }
    if (f._base != std::ios_base::dec) s.setf(std::ios_base::showbase);
}

inline void Setup(std::ostream& os, const BoundForm<double>& bf)
{ SetupFloat(os,bf.f); }

inline void Setup(std::ostream& os, const BoundForm<long double>& bf)
{ SetupFloat(os,bf.f); }

inline void Setup(std::ostream& os, const BoundForm<float>& bf)
{ SetupFloat(os,bf.f); }

inline void Setup(std::ostream& os, const BoundForm<std::complex<double> >& bf)
{ SetupFloat(os,bf.f); }

inline void Setup(
    std::ostream& os, const BoundForm<std::complex<long double> >& bf)
{ SetupFloat(os,bf.f); }

inline void Setup(std::ostream& os, const BoundForm<std::complex<float> >& bf)
{ SetupFloat(os,bf.f); }

inline void Setup(std::ostream& os, const BoundForm<int>& bf)
{ SetupInt(os,bf.f); }

inline void Setup(std::ostream& os, const BoundForm<short>& bf)
{ SetupInt(os,bf.f); }

inline void Setup(std::ostream& os, const BoundForm<long>& bf)
{ SetupInt(os,bf.f); }

inline void Setup(std::ostream& os, const BoundForm<unsigned int>& bf)
{ SetupInt(os,bf.f); }

inline void Setup(std::ostream& os, const BoundForm<unsigned short>& bf)
{ SetupInt(os,bf.f); }

inline void Setup(std::ostream& os, const BoundForm<unsigned long>& bf)
{ SetupInt(os,bf.f); }

template <typename T>
inline void Setup(std::ostream& os, const BoundForm<T>& bf)
{ SetupFloat(os,bf.f); }

template <typename T>
inline std::ostream& operator<<(std::ostream& os, const BoundForm<T>& bf)
{
    std::ostringstream s;
    Setup(s,bf);
    s << bf.val;
    if (bf.f._n_trail>0) s << std::string(bf.f._n_trail,bf.f._trail_ch);
    os << s.str();
    return os;
}

#endif
