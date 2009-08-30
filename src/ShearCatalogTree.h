#ifndef ShearCatalogTree_H
#define ShearCatalogTree_H

#include <vector>
#include "ShearCatalog.h"
#include "Bounds.h"

class ShearCatalogTree
{
  public :

    // Make from a ShearCatalog
    ShearCatalogTree(const ShearCatalog& _incat);
    ~ShearCatalogTree();

    int FindNearestTo(const Position& pos);

  private :

    const ShearCatalog& incat;
    struct Node;
    Node* top;
};

#endif
