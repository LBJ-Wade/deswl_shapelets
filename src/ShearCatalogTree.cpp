

#include "ShearCatalogTree.h"

struct PosWithIndex
{
  Position pos;
  int index;
};

struct CompareX
{
  bool operator()(const PosWithIndex& p1, const PosWithIndex& p2)
  { return p1.pos.GetX() < p2.pos.GetX(); }
};

struct CompareY
{
  bool operator()(const PosWithIndex& p1, const PosWithIndex& p2)
  { return p1.pos.GetY() < p2.pos.GetY(); }
};

struct ShearCatalogTree::Node
{
  typedef std::vector<PosWithIndex>::iterator it;
  Node(const it& start, const it& end)
  {
    int n = end-start;
    if (n == 1)
    { 
      b += start->pos;
      left = 0;
      right = 0;
      pos = start->pos;
      index = start->index;
    }
    else
    {
      for(it i=start;i!=end;++i) b += i->pos;
      it middle = start + n/2;
      if (b.IsWide())
	std::nth_element(start,middle,end,CompareX());
      else
	std::nth_element(start,middle,end,CompareY());
      left = new Node(start,middle);
      right = new Node(middle,end);
      pos = b.Center();
      index = -1;
    }
  }

  ~Node() 
  {
    if (left) delete left;
    if (right) delete right;
  }

  void FindNearestTo(const Position& p, int& besti, double& best)
  {
    xxdbg<<"Node::FindNearestTo: "<<p<<"  "<<besti<<"  "<<best<<std::endl;
    xxdbg<<"bounds = "<<b<<std::endl;
    if (!left)
    {
      Assert(!right);
      xxdbg<<"leaf node\n";
      // If we are on a leaf node, then this becomes a simple check:
      double d = std::abs(p-pos);
      if (best < 0. || d < best)
      { 
	if (best >= 0.)
	{
	  xxdbg<<"Updating best.\n";
	}
	besti = index; 
	best = d; 
      }
      if (d < 1.)
      {
	xdbg<<"Found close match.  d = "<<d<<std::endl;
      }
      else if (d < 5.)
      {
	xdbg<<"Found moderately close match.  d = "<<d<<std::endl;
      }
    }
    else 
    {
      xxdbg<<"regular node\n";
      // Otherwise need to recurse down to sub-nodes.
      
      // First see if we can trivially stop here.
      if (best >= 0.)
      {
	if (p.GetX() < b.GetXMin() - best) return;
	if (p.GetX() > b.GetXMax() + best) return;
	if (p.GetY() < b.GetYMin() - best) return;
	if (p.GetY() > b.GetYMax() + best) return;
      }
      xxdbg<<"no early exit\n";

      // Recurse into sub-nodes smartly -- try the node that is more
      // likely to have the nearest location first.
      // Use NY distances, since they are faster, and don't need to be exact.
      double d1 = std::abs(p.GetX()-left->pos.GetX());
      d1 += std::abs(p.GetY()-left->pos.GetY());
      double d2 = std::abs(p.GetX()-right->pos.GetX());
      d2 += std::abs(p.GetY()-right->pos.GetY());
      xxdbg<<"d1,d2 = "<<d1<<"  "<<d2<<std::endl;
      if (d1 < d2)
      {
	left->FindNearestTo(p,besti,best);
	right->FindNearestTo(p,besti,best);
      }
      else
      {
	right->FindNearestTo(p,besti,best);
	left->FindNearestTo(p,besti,best);
      }
    }
  }
    
  Node* left;
  Node* right;
  Position pos;
  int index;
  Bounds b;
};

ShearCatalogTree::ShearCatalogTree(const ShearCatalog& _incat) : incat(_incat)
{
  std::vector<PosWithIndex> data(incat.size());
  for(int i=0;i<incat.size();++i)
  { 
    data[i].pos = incat.pos[i];
    data[i].index = i;
  }
  top = new Node(data.begin(),data.end());
}

ShearCatalogTree::~ShearCatalogTree()
{
  delete top;
}

int ShearCatalogTree::FindNearestTo(const Position& pos)
{
  xdbg<<"FindNearest for pos = "<<pos<<std::endl;
  int index = -1;
  double best = -1.;
  top->FindNearestTo(pos,index,best);
  xdbg<<"Found: index = "<<index<<", best = "<<best<<std::endl;
  Assert(index >= 0);
  Assert(index < incat.size());
  xdbg<<"incat.pos["<<index<<"] = "<<incat.pos[index]<<std::endl;
  xdbg<<"actual distance = "<<std::abs(pos-incat.pos[index])<<std::endl;
  Assert(std::abs(pos-incat.pos[index]) == best);
}

