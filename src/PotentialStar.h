//---------------------------------------------------------------------------
#ifndef PotentialStarH
#define PotentialStarH
//---------------------------------------------------------------------------

#include "Bounds.h"
#include <string>

class PotentialStar {

public:
  PotentialStar(Position _pos, double _mag, double _size, long _index,
      std::string _line) :
    itspos(_pos),itsmag(_mag),itssize(_size),itsindex(_index),
    itsline(_line) {}
  ~PotentialStar() {}

  const Position& GetPos() const {return itspos;}
  double GetMag() const { return itsmag; }
  double GetSize() const { return itssize; }
  long GetIndex() const { return itsindex; }
  const std::string& GetLine() const {return itsline;}

  void SetSize(double newsize) {itssize = newsize;}

  bool MagCompare(const PotentialStar* rhs) const
    {return itsmag < rhs->itsmag;}
  bool SizeCompare(const PotentialStar* rhs) const
    {return itssize < rhs->itssize;}

private:

  Position itspos;
  double itsmag,itssize;
  long itsindex;
  std::string itsline;

};


#endif

