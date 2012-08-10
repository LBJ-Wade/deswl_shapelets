//---------------------------------------------------------------------------
#ifndef PotentialStarH
#define PotentialStarH
//---------------------------------------------------------------------------

#include <string>
#include "dbg.h"
#include "Bounds.h"

class PotentialStar {

public:
  PotentialStar(Position pos, double mag, double nu,
		double size, double sg, long index,
                  const std::string& line) :
    _pos(pos), _mag(mag), _nu(nu), _size(size), _index(index), _sg(sg), 
    _line(line) 
    {}

    ~PotentialStar() {}

    const Position& getPos() const { return _pos; }

    double getMag() const { return _mag; }

    double getSg() const { return _sg; }

    double getNu() const { return _nu; }

    double getSize() const { return _size; }

    long getIndex() const { return _index; }

    const std::string& getLine() const { return _line; }

    void setSize(double newsize) { _size = newsize; }

    bool isBrighterThan(const PotentialStar* rhs) const
    { return _mag < rhs->_mag; }

    bool isSmallerThan(const PotentialStar* rhs) const
    { return _size < rhs->_size; }

private:

    Position _pos;
    double _mag;
    double _size;
    long _index;
    double _sg;
    double _nu;
    std::string _line;

};


#endif

