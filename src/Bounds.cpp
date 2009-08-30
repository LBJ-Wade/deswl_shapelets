#include "Bounds.h"
#include <vector>

/*

int Position::ChipNum(Position *rawcoord) const
// Finds the chipnum (aka ampid) corresponding to chippos pos
{
  if (rawcoord) *rawcoord = *this;
  if(x < 2400) {
    if(y < 2400) return 4;
    else {
      if (rawcoord) rawcoord->y -= 2808;
      return 1;
    }
  }
  else {
    if (rawcoord) rawcoord->x -= 2808;
    if(y < 2400) return 3;
    else {
      if (rawcoord) rawcoord->y -= 2808;
      return 2;
    }
  }
}

*/

void Bounds::operator+=(const Position &pos)
// Expand the bounds to include the given position.
{
  if (defined) {
    if (pos.GetX() < xmin) xmin = pos.GetX();
    else if (pos.GetX() > xmax) xmax = pos.GetX();
    if (pos.GetY() < ymin) ymin = pos.GetY();
    else if (pos.GetY() > ymax) ymax = pos.GetY();
  }
  else {
    xmin = xmax = pos.GetX();
    ymin = ymax = pos.GetY();
    defined = true;
  }
}

void Bounds::operator+=(const Bounds &rec)
// Expand the bounds to include the given rectangle (ie. bounds)
{
  if (!rec.IsDefined()) return;
  if (defined) {
    if (rec.GetXMin() < xmin) xmin = rec.GetXMin();
    if (rec.GetXMax() > xmax) xmax = rec.GetXMax();
    if (rec.GetYMin() < ymin) ymin = rec.GetYMin();
    if (rec.GetYMax() > ymax) ymax = rec.GetYMax();
  }
  else {
    *this = rec;
    defined = true;
  }
}

Position Bounds::Center() const
{
  return Position((xmin + xmax)/2.,(ymin + ymax)/2.);
}

Bounds Bounds::operator&(const Bounds &rhs) const
{
  Bounds temp(xmin<rhs.xmin ? rhs.xmin : xmin,
    xmax>rhs.xmax ? rhs.xmax : xmax,
    ymin<rhs.ymin ? rhs.ymin : ymin,
    ymax>rhs.ymax ? rhs.ymax : ymax);
  if (temp.xmin>temp.xmax || temp.ymin>temp.ymax) return Bounds();
  else return temp;
}

void Bounds::Snap(double d)
{
  if (defined) {
    xmax = d*ceil(xmax/d);
    xmin = d*floor(xmin/d);
    ymax = d*ceil(ymax/d);
    ymin = d*floor(ymin/d);
  }
}

void Bounds::AddXBorder(double d)
{
  if (defined) { xmax += d; xmin -= d; }
}

void Bounds::AddYBorder(double d)
{
  if (defined) { ymax += d; ymin -= d; }
}

void Bounds::AddBorder(double d)
{
  if (defined) { xmax += d; xmin -= d; ymax += d; ymin -= d; }
}

std::vector<Bounds> Bounds::Divide(size_t nx,size_t ny) const
{
  if (!defined) return std::vector<Bounds>(nx*ny);
  std::vector<Bounds> temp;
  temp.reserve(nx*ny);
  std::vector<double> x(nx+1);
  std::vector<double> y(ny+1);
  x[0] = xmin;  x[nx] = xmax;
  y[0] = ymin;  y[ny] = ymax;
  double xstep = (xmax-xmin)/nx;
  double ystep = (ymax-ymin)/ny;
  for(size_t i=1;i<nx;i++) x[i] = x[0]+i*xstep;
  for(size_t j=1;j<ny;j++) y[j] = y[0]+j*ystep;
  for(size_t i=0;i<nx;i++) for(size_t j=0;j<ny;j++)
    temp.push_back(Bounds(x[i],x[i+1],y[j],y[j+1]));
  return temp;
}
