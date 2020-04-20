/******************************************************************************
File      : QuadTree.cpp
Author    : L. Shawn Matott
Copyright : 2004, L. Shawn Matott

Encapsulates a quad tree; a set of parameters whose values at each level of the
tree are evenly spaced between an upper and lower limit.

Version History
03-12-04    lsm   created
08-17-04    lsm   RAM fragmentation fixes, Added reporting of memory allocations
******************************************************************************/
#include <stdio.h>
#include <math.h>

#include "QuadTree.h"

#include "Exception.h"

void DestroyQuadTree(QuadNodeStruct * pTree);
void ExpandQuadTree(QuadNodeStruct * pTree);
int  GetQuadTreeLevel(int lvl, int idx, double * pVals);

/******************************************************************************
DestroyQuadTree()

Recursively frees up memory used in the tree.
******************************************************************************/
void DestroyQuadTree(QuadNodeStruct * pTree)
{
   if(pTree == NULL){ return;}
   DestroyQuadTree(pTree->pLeft);
   DestroyQuadTree(pTree->pRight);
   delete pTree->pLeft;
   delete pTree->pRight;
}/* end DestroyQuadTree() */

/******************************************************************************
ExpandQuadTree()

Add one level to the tree.
******************************************************************************/
void ExpandQuadTree(QuadNodeStruct * pTree)
{
   QuadNodeStruct * pLeaf;

   if(pTree == NULL){ return;}
   if(pTree->pLeft == NULL) //expand terminal leaf
   {
      NEW_PRINT("QuadNodeStruct", 1);
      pTree->pLeft = new QuadNodeStruct;

      NEW_PRINT("QuadNodeStruct", 1);
      pTree->pRight = new QuadNodeStruct;
      MEM_CHECK(pTree->pRight);

      pLeaf = pTree->pLeft;
      pLeaf->pLeft  = NULL;
      pLeaf->pRight = NULL;
      pLeaf->upr    = pTree->mid;
      pLeaf->lwr    = pTree->lwr;
      pLeaf->mid    = 0.5*(pTree->mid + pTree->lwr);
      pLeaf->lvl    = pTree->lvl + 1;

      pLeaf = pTree->pRight;
      pLeaf->pLeft  = NULL;
      pLeaf->pRight = NULL;
      pLeaf->upr    = pTree->upr;
      pLeaf->lwr    = pTree->mid;
      pLeaf->mid    = 0.5*(pTree->mid + pTree->upr);
      pLeaf->lvl    = pTree->lvl + 1;
   }
   else //perform recursion until terminal leaf is found
   {
      ExpandQuadTree(pTree->pLeft);
      ExpandQuadTree(pTree->pRight);
   }
}/* end ExpandQuadTree() */

/******************************************************************************
GetQuadTreeLevel()

Recursivelt copy the mid values of the tree at the specified level into pVals.
******************************************************************************/
int  GetQuadTreeLevel(QuadNodeStruct * pTree, int theLvl, double * pVals, int curIdx)
{
   if(pTree == NULL){ return 0;}
   if(pTree->lvl == theLvl)
   {
      pVals[curIdx] = pTree->mid;
      return (curIdx+1);
   }
   else
   {
      curIdx = GetQuadTreeLevel(pTree->pLeft, theLvl, pVals, curIdx);
      curIdx = GetQuadTreeLevel(pTree->pRight, theLvl, pVals, curIdx);
      return curIdx;
   }
}/* end GetQuadTreeLevel() */

/******************************************************************************
Destroy()
******************************************************************************/
void QuadTree::Destroy(void)
{
   DestroyQuadTree(m_pTree);
   delete m_pTree;
   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR
******************************************************************************/
QuadTree::QuadTree(void)
{
   m_pTree = NULL;
   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
CTOR
******************************************************************************/
QuadTree::QuadTree(double lwr, double upr)
{
   NEW_PRINT("QuadNodeStruct", 1);
   m_pTree = new QuadNodeStruct;
   MEM_CHECK(m_pTree);

   m_pTree->lwr = lwr;
   m_pTree->upr = upr;
   m_pTree->mid = 0.5 * (lwr + upr);
   m_pTree->lvl = 0;
   m_pTree->pLeft = NULL;
   m_pTree->pRight = NULL;
   m_NumLvls = 0;

   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
Init()

Initialize the Quad Tree. If already exists, tree will be destroyed.
******************************************************************************/
void QuadTree::Init(double lwr, double upr)
{
   DestroyQuadTree(m_pTree);

   NEW_PRINT("QuadNodeStruct", 1);
   m_pTree = new QuadNodeStruct;
   MEM_CHECK(m_pTree);

   m_pTree->lwr = lwr;
   m_pTree->upr = upr;
   m_pTree->mid = 0.5 * (lwr + upr);
   m_pTree->lvl = 0;
   m_pTree->pLeft = NULL;
   m_pTree->pRight = NULL;
   m_NumLvls = 0;   
} /* end default Init() */

/******************************************************************************
Expand()

Adds a level to the tree.
******************************************************************************/
void QuadTree::Expand(void)
{
   ExpandQuadTree(m_pTree);
   m_NumLvls++;
} /* end Expand() */

/******************************************************************************
GetLevel()

Returns an array stuffed with the mid points of the desired level.
******************************************************************************/
double * QuadTree::GetLevel(int lvl)
{
   double * pLvl;
   int size;

   if(lvl > m_NumLvls){ return NULL;}

   size = 1 << lvl; //2^lvl

   NEW_PRINT("double", size);
   pLvl = new double[size];
   MEM_CHECK(pLvl);

   GetQuadTreeLevel(m_pTree, lvl, pLvl, 0);

   return pLvl;
} /* end GetLevel() */

/******************************************************************************
GetTreeCombo()

Returns an array containing a combination of parameters from the 'lvl' level 
of the tree.
******************************************************************************/
double * GetTreeCombo(int lvl, int idx, QuadTree * pList, int size)
{
   int j, num_combos, a, place;
   double * tmp, * combo;
   
   num_combos = (int)(pow((double)(1 << lvl), size));   

   //check index....
   if(idx >= num_combos){ return NULL;};

   NEW_PRINT("double", size);
   combo = new double[size];
   MEM_CHECK(combo);
   
   //compute combination 'idx' at the given level 'lvl'
   //printf("COMBO %d = (", idx);
   //for each parameter, determine appropriate index
   for(j = (size-1); j >= 0; j--)
   {
      tmp = pList[j].GetLevel(lvl);
      a = 0;
      place = (int)((double)pow((double)(1 << lvl), j));
      while((a * place) <= idx){a++;}
      a--;
      combo[j] = tmp[a];
      idx -= (a * place);
      delete [] tmp;
      //printf("%+lf", combo[j]);
      //if(j > 0){printf(", ");}
   }/* end for() */      
   //printf(")\n");
   
   return combo;
}/* end GetCombination() */

