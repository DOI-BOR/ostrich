/******************************************************************************
File      : QuadTree.h
Author    : L. Shawn Matott
Copyright : 2004, L. Shawn Matott

Encapsulates a Quadtree. A set of parameters whose values at each level of the
tree are evenly spaced between an upper and lower limit.

Version History
03-12-04    lsm   created
08-17-04    lsm   RAM fragmentation fixes, Added reporting of memory allocations
******************************************************************************/
#ifndef QUAD_TREE_H
#define QUAD_TREE_H

#include "MyHeaderInc.h"

//Quad tree node definition
typedef struct QUAD_NODE_STRUCT
{
   int lvl;
   double upr;
   double lwr;
   double mid;
   struct QUAD_NODE_STRUCT * pLeft;
   struct QUAD_NODE_STRUCT * pRight;
}QuadNodeStruct;

/******************************************************************************
class QuadTree

Encapsulates a Quad tree.
******************************************************************************/
class QuadTree
{
   private:
      QuadNodeStruct * m_pTree;
      int m_NumLvls;

   public:      
      QuadTree(void);
      QuadTree(double lwr, double upr);
      ~QuadTree(void){ DBG_PRINT("QuadTree::DTOR"); Destroy(); }
      void Destroy(void);

      void Init(double lwr, double upr);
      void Expand(void);
      double * GetLevel(int lvl);      
}; /* end class QuadTree */

extern "C" {
double * GetTreeCombo(int lvl, int idx, QuadTree * pList, int size);
}

#endif /* QUAD_TREE_H */

