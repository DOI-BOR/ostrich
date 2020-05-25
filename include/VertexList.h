/******************************************************************************
File      : VertexList.h
Author    : L. Shawn Matott
Copyright : 2004, L. Shawn Matott

Encapsulates a list of vertices which make up a geometric shape. Two kinds of
vertex lists are defined, one contains constant values for all vertices, while
the other (augmented) list can contain parameters in place of actual values.

Version History
11-29-04    lsm   Created
******************************************************************************/
#ifndef VERTEX_LIST_H
#define VERTEX_LIST_H

#include "MyHeaderInc.h"

//forward decs
class ParameterABC;
class TiedParamABC;

extern "C" {
   /* define a structure to contain geometry vertices */
   typedef struct VERTEX_LIST_STRUCT
   { 
      double x, y, z;
      struct VERTEX_LIST_STRUCT * pNxt;
   }VertexList;

   /* define a structure to contain augmented geometry vertices */
   typedef struct AUG_VERT_LIST_STRUCT
   { 
      ParameterABC * px, * py, * pz;
      TiedParamABC * tx, * ty, * tz;
      double x, y, z;
      struct AUG_VERT_LIST_STRUCT * pNxt;
   }AugVertexList;

   /* define a structure to contain augmented circle vertices */
   typedef struct AUG_CIRCLE_STRUCT
   { 
      ParameterABC * px, * py, * pz, * pr;
      TiedParamABC * tx, * ty, * tz, * tr;
      double x, y, z, r;
   }AugCircle;

   void ConvertAugCircle(AugCircle * pAug, Circle2D * pCirc, double * pZ);
   void ConvertAugVertex(AugVertexList * pVert, VertexList * pConv);
   void DestroyVertexList(VertexList * pList);
   void DestroyAugVertexList(AugVertexList * pList);
}/* end extern */

#endif /* VERTEX_LIST_H */
