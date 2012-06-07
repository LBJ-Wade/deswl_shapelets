

#include "ShearCatalogTree.h"

struct PosWithIndex
{
    Position pos;
    int index;
};

struct CompareX
{
    bool operator()(const PosWithIndex& p1, const PosWithIndex& p2)
    { return p1.pos.getX() < p2.pos.getX(); }
};

struct CompareY
{
    bool operator()(const PosWithIndex& p1, const PosWithIndex& p2)
    { return p1.pos.getY() < p2.pos.getY(); }
};

struct ShearCatalogTree::Node
{
    typedef std::vector<PosWithIndex>::iterator it;
    Node(const it& start, const it& end)
    {
        int n = end-start;
        if (n == 1) { 
            b += start->pos;
            left = 0;
            right = 0;
            pos = start->pos;
            index = start->index;
        } else {
            for(it i=start;i!=end;++i) b += i->pos;
            it middle = start + n/2;
            if (b.isWide())
                std::nth_element(start,middle,end,CompareX());
            else
                std::nth_element(start,middle,end,CompareY());
            left = new Node(start,middle);
            right = new Node(middle,end);
            pos = b.getCenter();
            index = -1;
        }
    }

    ~Node() 
    {
        if (left) delete left;
        if (right) delete right;
    }

    void findNearestTo(const Position& p, int& iBest, double& best) const
    {
        if (!left) {
            Assert(!right);
            // If we are on a leaf node, then this becomes a simple check:
            double d = std::abs(p-pos);
            if (best < 0. || d < best) { 
                iBest = index; 
                best = d; 
            }
            if (d < 1.) {
                xdbg<<"Found close match.  d = "<<d<<std::endl;
            } else if (d < 5.) {
                xdbg<<"Found moderately close match.  d = "<<d<<std::endl;
            }
        } else {
            // Otherwise need to recurse down to sub-nodes.

            // First see if we can trivially stop here.
            if (best >= 0.) {
                if (p.getX() < b.getXMin() - best) return;
                if (p.getX() > b.getXMax() + best) return;
                if (p.getY() < b.getYMin() - best) return;
                if (p.getY() > b.getYMax() + best) return;
            }

            // Recurse into sub-nodes smartly -- try the node that is more
            // likely to have the nearest location first.
            // Use NY distances, since they are faster,
            // and don't need to be exact.
            double d1 = std::abs(p.getX()-left->pos.getX());
            d1 += std::abs(p.getY()-left->pos.getY());
            double d2 = std::abs(p.getX()-right->pos.getX());
            d2 += std::abs(p.getY()-right->pos.getY());
            if (d1 < d2) {
                left->findNearestTo(p,iBest,best);
                right->findNearestTo(p,iBest,best);
            } else {
                right->findNearestTo(p,iBest,best);
                left->findNearestTo(p,iBest,best);
            }
        }
    }

    Node* left;
    Node* right;
    Position pos;
    int index;
    Bounds b;
};

ShearCatalogTree::ShearCatalogTree(const ShearCatalog& cat) : _cat(cat)
{
    const int ngal = _cat.size();
    std::vector<PosWithIndex> data(ngal);
    for(int i=0;i<ngal;++i) { 
        data[i].pos = _cat.getPos(i);
        data[i].index = i;
    }
    _top = new Node(data.begin(),data.end());
}

ShearCatalogTree::~ShearCatalogTree()
{
    delete _top;
}

int ShearCatalogTree::findNearestTo(const Position& pos) const
{
    xdbg<<"FindNearest for pos = "<<pos<<std::endl;
    int index = -1;
    double best = -1.;
    _top->findNearestTo(pos,index,best);
    xdbg<<"Found: index = "<<index<<", best = "<<best<<std::endl;
    Assert(index >= 0);
    Assert(index < int(_cat.size()));
    xdbg<<"incat.pos["<<index<<"] = "<<_cat.getPos(index)<<std::endl;
    xdbg<<"actual distance = "<<std::abs(pos-_cat.getPos(index))<<std::endl;
    xdbg<<"diff = "<<std::abs(std::abs(pos-_cat.getPos(index))-best)<<std::endl;
    Assert(std::abs(std::abs(pos-_cat.getPos(index)) - best) < 1.e-5);
    return index;
}

