/******************************************************************************
File      : VertexList.cpp
Author    : L. Shawn Matott
Copyright : 2004, L. Shawn Matott

Encapsulates a list of vertices which make up a geometric shape. Two kinds of
vertex lists are defined, one contains constant values for all vertices, while
the other (augmented) list can contain parameters in place of actual values.

Version History
11-29-04    lsm   Created
******************************************************************************/
#include "VertexList.h"
#include "ParameterABC.h"
#include "TiedParamABC.h"

/******************************************************************************
ConvertAugVertex()

Converts an augmented vertex to a normal vertex.
******************************************************************************/
void ConvertAugVertex(AugVertexList * pVert, VertexList * pConv)
{
   if(pVert->px != NULL){ pConv->x = pVert->px->GetTransformedVal();}
   else if(pVert->tx != NULL){pConv->x = pVert->tx->GetEstVal();}
   else (pConv->x = pVert->x);

   if(pVert->py != NULL){ pConv->y = pVert->py->GetTransformedVal();}
   else if(pVert->ty != NULL){pConv->y = pVert->ty->GetEstVal();}
   else (pConv->y = pVert->y);

   if(pVert->pz != NULL){ pConv->z = pVert->pz->GetTransformedVal();}
   else if(pVert->tz != NULL){pConv->z = pVert->tz->GetEstVal();}
   else (pConv->z = pVert->z);
}/* end ConvertAugVertex() */

/******************************************************************************
ConvertAugCircle()

Converts an augmented circle to a normal circle.
******************************************************************************/
void ConvertAugCircle(AugCircle * pAug, Circle2D * pCirc, double * pZ)
{
   if     (pAug->px != NULL){ pCirc->x = pAug->px->GetTransformedVal();}
   else if(pAug->tx != NULL){ pCirc->x = pAug->tx->GetEstVal();}
   else   (pCirc->x = pAug->x);

   if     (pAug->py != NULL){ pCirc->y = pAug->py->GetTransformedVal();}
   else if(pAug->ty != NULL){ pCirc->y = pAug->ty->GetEstVal();}
   else   (pCirc->y = pAug->y);

   if     (pAug->pr != NULL){ pCirc->r = pAug->pr->GetTransformedVal();}
   else if(pAug->tr != NULL){ pCirc->r = pAug->tr->GetEstVal();}
   else   (pCirc->r = pAug->r);

   if     (pAug->pz != NULL){ *pZ = pAug->pz->GetTransformedVal();}
   else if(pAug->tz != NULL){ *pZ = pAug->tz->GetEstVal();}
   else   (*pZ = pAug->z);
}/* end ConvertAugCircle() */

/******************************************************************************
DestroyVertexList()

Recursively frees up memory in a vertex list.
******************************************************************************/
void DestroyVertexList(VertexList * pList)
{
   if(pList == NULL){ return;}
   DestroyVertexList(pList->pNxt);
   delete pList;
}/* end DestroyVertexList()*/

/******************************************************************************
DestroyAugVertexList()

Recursively frees up memory in an augmented vertex list.
******************************************************************************/
void DestroyAugVertexList(AugVertexList * pList)
{
   if(pList == NULL){ return;}
   DestroyAugVertexList(pList->pNxt);
   delete pList;
}/* end DestroyAugVertexList()*/

