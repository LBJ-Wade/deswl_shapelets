#ifndef ShearCatalogTree_H
#define ShearCatalogTree_H

#include <vector>
#include "dbg.h"
#include "ShearCatalog.h"
#include "Bounds.h"

class ShearCatalogTree
{
public :

    // Make from a ShearCatalog
    ShearCatalogTree(const ShearCatalog& cat);
    ~ShearCatalogTree();

    int findNearestTo(const Position& pos) const;

private :

    const ShearCatalog& _cat;
    struct Node;
    Node* _top;
};

#endif
