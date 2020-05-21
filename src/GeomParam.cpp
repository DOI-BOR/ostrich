/******************************************************************************
File      : GeomParam.cpp
Author    : L. Shawn Matott
Copyright : 2004, L. Shawn Matott

Encapsulates a 'geometry' parameter. Geometry parameters are variables in the model 
which are composed on one or more spatial vertices. The ABC for the geometry 
parameters encapsulates the interface used by other Ostrich modules, allowing 
various specific geometry parameter relationships (line3, poly2, poly3, etc.) to be 
implemented as needed with minimal code change (just need to add the specific 
geometry parameter class and some additional input file parsing). The purpose of a
geometry parameter is to facilitate changes to geometric model properties as part of
a calibration or optimization exercise. Ostrich will ensure that all geometry parameters
are topologically correct in that:
   1) vertices will be automatically inserted if elements overlap
   2) if a given ordering of vertices is not valid, the polygon vertices will be 
   randomly reordered until a valid polygon is found

These specific geometry-parameter classes are supported:

GeomParamLine3 : a polyline containing a set of (x,y,z) values, where (x,y) are 
the spatial coordinates and z is a non-geometric value (i.e. head). When 
vertices are inserted, the z-value of the new vertex is interpolated.

GeomParamPoly3  : a polygon containing a set of (x,y,z) values, where (x,y) are 
the spatial coordinates and z is a non-geometric value (i.e. head). When 
vertices are inserted, the z-value of the new vertex is interpolated.

GeomParamPoly2  : a polygon containing a set of (x,y) values, where (x,y) are 
the spatial coordinates.

Version History
11-29-04    lsm   Created
10-19-05    lsm   Replaced calls to rand() with MyRand()
******************************************************************************/
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "GeomParamABC.h"
#include "ParameterABC.h"
#include "VertexList.h"

#include "Exception.h"
#include "Utility.h"

/******************************************************************************
GeomParamLine3::Destroy()
******************************************************************************/
void GeomParamLine3::Destroy(void)
{
   delete [] m_pName;

   //destroy vertices
   DestroyAugVertexList(m_pInit);
   DestroyVertexList(m_pFixed);

   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR (GeomParamLine3)
******************************************************************************/
GeomParamLine3::GeomParamLine3(void)
{
   m_pName = NULL;
   m_pInit = NULL;
   m_pFixed = NULL;
   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
CTOR (GeomParamLine3)
******************************************************************************/
GeomParamLine3::GeomParamLine3(IroncladString name)
{
   int len;

   m_pInit = NULL;
   m_pFixed = NULL;

   len = (int)strlen(name) + 10;
   NEW_PRINT("char", len);
   m_pName = new char[len];
   MEM_CHECK(m_pName);

   strcpy(m_pName, name);

   IncCtorCount();
} /* end CTOR */

/******************************************************************************
GeomParamLine3::Convert()

Convert augmented vertices to normal vertices. Augmented vertices may contain
parameter names instead of constant values, such as when geometry of a given
element is being calibrated.
******************************************************************************/
void GeomParamLine3::Convert(void)
{
   AugVertexList * pCur;
   VertexList * pTmp;

   //first time? then allocated space for fixed list
   if((m_pFixed == NULL) && (m_pInit != NULL))
   {
      NEW_PRINT("VertexList", 1);
      m_pFixed = new VertexList;
      m_pFixed->pNxt = NULL;
      MEM_CHECK(m_pFixed);
   }

   //convert each element in list
   pCur = m_pInit;
   pTmp = m_pFixed;
   while(pCur != NULL)
   {
      ConvertAugVertex(pCur, pTmp);
      pCur = pCur->pNxt;
      if(pCur != NULL)          
      {
         if(pTmp->pNxt == NULL)
         {
            NEW_PRINT("VertexList", 1);
            pTmp->pNxt = new VertexList;         
            pTmp = pTmp->pNxt;
            MEM_CHECK(pTmp);
            pTmp->pNxt = NULL;
         }
         else
         {
            pTmp = pTmp->pNxt;
         }
      }
      else
      {
         DestroyVertexList(pTmp->pNxt);
         pTmp->pNxt = NULL;
      }      
   }/* end while() */   
} /* end Convert() */

/******************************************************************************
GeomParamLine3::Reorder()

Reorder vertices if they overlap one another. Returns true if reorder is 
successful, false otherwise.
******************************************************************************/
bool GeomParamLine3::Reorder(void)
{
   Point2D pt;
   Segment2D seg1, seg2;
   AugVertexList * pCurAug;   
   VertexList * pCur, * pTmp;
   double xupr, xlwr, yupr, ylwr, tmpx, tmpy;
   bool mustCheck, backProp, done;
   int count, numTries, maxTries, i, r, j;

   pCurAug = m_pInit;
   pCur = m_pFixed;

   yupr = xupr = -1.00;
   ylwr = xlwr = +1.00;
   //first, determine if reordering check is even needed
   mustCheck = false;
   backProp = true;
   count = 0;
   while(pCurAug != NULL)
   {
      count++;
      //only check shapes with parameters attached
      if((pCurAug->px != NULL) || (pCurAug->py != NULL) || 
         (pCurAug->tx != NULL) || (pCurAug->ty != NULL))
      {
         mustCheck = true;
      }
      /* only back-propagate if all vertices are params w/ same bounds */
      if((pCurAug->px == NULL) || (pCurAug->py == NULL))
      {
         backProp = false;
      }
      else //check if same bounds for all
      {
         //first time?
         if(xupr < xlwr){ xupr = pCurAug->px->GetUprBnd(); xlwr = pCurAug->px->GetLwrBnd();}
         if(yupr < ylwr){ yupr = pCurAug->py->GetUprBnd(); ylwr = pCurAug->py->GetLwrBnd();}

         if((pCurAug->px->GetUprBnd() != xupr) || (pCurAug->px->GetLwrBnd() != xlwr)){ backProp = false;}
         if((pCurAug->py->GetUprBnd() != yupr) || (pCurAug->py->GetLwrBnd() != ylwr)){ backProp = false;}
      }
      pCurAug = pCurAug->pNxt;   
   }/* end while() */

   if(mustCheck == false){ return true;}   
   
   done = false;
   numTries = 0;

   maxTries = 1;
   for(i = count; i > 1; i--){ maxTries *= i;}
   maxTries *= 2;

   while((done == false) && (numTries < maxTries))
   {
      done = true;      
      //compare all segments to see if there is overlap
      for(pCur = m_pFixed; pCur->pNxt != NULL; pCur = pCur->pNxt)      
      {
         seg1.p1.x = pCur->x;
         seg1.p1.y = pCur->y;
         seg1.p2.x = pCur->pNxt->x;
         seg1.p2.y = pCur->pNxt->y;

         for(pTmp = pCur->pNxt; pTmp != NULL; pTmp = pTmp->pNxt)
         {
            seg2.p1.x = pTmp->x;
            seg2.p1.y = pTmp->y;
            if(pTmp->pNxt != NULL)
            {
               seg2.p2.x = pTmp->pNxt->x;
               seg2.p2.y = pTmp->pNxt->y;

               if(SegIntersect(&seg1, &seg2, &pt) == BOTHSEG)
               {
                  done = false;
                  break;
               }
            }/* end if() */
         }/* end for() */
         if(done == false){ break;}      
      }/* end for() */

      //randomly reorder, if overlap detected      
      if(done == false)
      {
         pCur = m_pFixed->pNxt;
         for(i = 1; i < (count-1); i++)
         {
            r = MyRand()%(count-i);
            pTmp = m_pFixed->pNxt;
            for(j = 0; j < r; j++){ pTmp = pTmp->pNxt;}
            tmpx = pTmp->x;
            tmpy = pTmp->y;
            pTmp->x = pCur->x;
            pTmp->y = pCur->y;
            pCur->x = tmpx;
            pCur->y = tmpy;
            pCur = pCur->pNxt;
         }/* end for() */
         numTries++;
      }/* end while() */  
   }/* end while() */

   //adjust parameter values to reflect new ordering
   if(backProp == false){ LogError(ERR_MISMATCH, "Can't back-propagate reordering");}
   else //parameters are interchangeable, so back-propagate reordering
   {
      pCurAug = m_pInit;
      pCur = m_pFixed;
      while(pCurAug != NULL)
      {
         pCurAug->px->SetEstVal(pCur->x);
         pCurAug->py->SetEstVal(pCur->y);
         pCurAug = pCurAug->pNxt;
         pCur = pCur->pNxt;
      }/* end while() */
   }/* end else() */
   
   //if(done == false){ printf("Re-order failed\n");}
   //else{ printf("Re-order succeeded\n");}
   return done;
}/* end Reorder() */

/******************************************************************************
GeomParamLine3::FixVertices()

Insert vertices wherever the geometry overlaps with pOther. Returns true if 
vertices are corrected successfully (or if correction was not necessary), 
false otherwise.
******************************************************************************/
bool GeomParamLine3::FixVertices(GeomParamABC * pOther)
{
   Segment2D seg;
   VertexList * pCur, * pFix, * pTmp;
   pCur = m_pFixed;
   double za, zb, zmid, xa, ya, xb, yb, xc, yc, d1, d2;

   //let the other geometry check each line segment of this geomettry
   while(pCur->pNxt != NULL)
   {
      xa = seg.p1.x = pCur->x;
      ya = seg.p1.y = pCur->y;
      za = pCur->z;

      xb = seg.p2.x = pCur->pNxt->x;
      yb = seg.p2.y = pCur->pNxt->y;
      zb = pCur->pNxt->z;

      //fix the segment, if necessary
      pFix = pOther->FixVertex(&seg);

      //insert corrections
      if(pFix != NULL)
      { 
         pTmp = pCur->pNxt;
         pCur->pNxt = pFix;
         while(pFix->pNxt != NULL)
         { 
            //interpolate z-values of inserted vertex
            xc = pFix->x;
            yc = pFix->y;
            d1 = sqrt((xc-xa)*(xc-xa)+(yc-ya)*(yc-ya));
            d2 = sqrt((xc-xb)*(xc-xb)+(yc-yb)*(yc-yb));
            zmid = za*(d2/(d1+d2)) + zb*(d1/(d1+d2));
            pFix->z = zmid;

            pFix = pFix->pNxt;
         }
         //interpolate z-values of inserted vertex
         xc = pFix->x;
         yc = pFix->y;
         d1 = sqrt((xc-xa)*(xc-xa)+(yc-ya)*(yc-ya));
         d2 = sqrt((xc-xb)*(xc-xb)+(yc-yb)*(yc-yb));
         zmid = za*(d2/(d1+d2)) + zb*(d1/(d1+d2));
         pFix->z = zmid;

         pFix->pNxt = pTmp;         
         //printf("Fixed geometry\n");
      }/* end if() */
      else
      {
         pCur = pCur->pNxt;
      }
   }/* end while() */
   return true;
}/* end FixVertices() */

/******************************************************************************
GeomParamLine3::FixVertex()

Insert vertices wherever the geometry overlaps with the line segment in pSeg. 
Returns a vertex list containing all inserted vertices, or NULL if no inserions
are needed.
******************************************************************************/
VertexList * GeomParamLine3::FixVertex(Segment2D * pSeg)
{
   Point2D pt;
   Segment2D mySeg;
   VertexList * pCur, * pMyNew, * pRetNew, * pRet, * pTmp;
   double za, zb, zmid, xa, ya, xb, yb, xc, yc, d1, d2;
   int test;

   pCur = m_pFixed;

   pRet = pMyNew = pRetNew = NULL;
   //compare all segments with pVert to see if there is overlap
   for(pCur = m_pFixed; pCur->pNxt != NULL; pCur = pCur->pNxt)      
   {
      xa = mySeg.p1.x = pCur->x;
      ya = mySeg.p1.y = pCur->y;
      za = pCur->z;

      xb = mySeg.p2.x = pCur->pNxt->x;
      yb = mySeg.p2.y = pCur->pNxt->y;
      zb = pCur->pNxt->z;

      test = SegIntersect(&mySeg, pSeg, &pt);

      if(test > NO_SEGS)
      {
         //need to insert a new node in mySeg?
         if((test == BOTHSEG) || (test == LEFTSEG))
         {
            NEW_PRINT("VertexList", 1);
            pMyNew  = new VertexList;         
            MEM_CHECK(pMyNew);

            xc = pMyNew->x = pt.x;
            yc = pMyNew->y = pt.y;
            d1 = sqrt((xc-xa)*(xc-xa)+(yc-ya)*(yc-ya));
            d2 = sqrt((xc-xb)*(xc-xb)+(yc-yb)*(yc-yb));
            zmid = za*(d2/(d1+d2)) + zb*(d1/(d1+d2));
            pMyNew->z = zmid;

            //insert the new vertex
            pTmp = pCur->pNxt;
            pCur->pNxt = pMyNew;
            pMyNew->pNxt = pTmp;
         }/* end if() */

         //need to insert a new node in pSeg?
         if((test == BOTHSEG) || (test == RGHTSEG))
         {         
            if(pRet == NULL)
            { 
               NEW_PRINT("VertexList", 1);
               pRet = new VertexList;
               MEM_CHECK(pRet);
               pRetNew = pRet;
               pRetNew->pNxt = NULL;
            }
            else
            {
               NEW_PRINT("VertexList", 1);
               pRetNew->pNxt = new VertexList;
               MEM_CHECK(pRet);

               pRetNew = pRetNew->pNxt;
               pRetNew->pNxt = NULL;
            }/* end else() */
            pRetNew->x = pt.x;
            pRetNew->y = pt.y;
         }/* end if() */

         //found an intersection, so quit checking
         break;   
      }/* end if() */
   }/* end for() */

   return pRet;
}/* end FixVertex() */

/******************************************************************************
GeomParamLine3::GetValAsStr()

Format the fixed geometry for output to model.
******************************************************************************/
void GeomParamLine3::GetValAsStr(UnmoveableString valStr)
{
   double x, y, px, py;
   double eps = 0.000001;
   char * pTok;
   VertexList * pCur;
   int i;

   pCur = m_pFixed;
   i = 0;
   px = py = -1.00;
   while(pCur != NULL)
   {
      x = pCur->x;
      y = pCur->y;

      if(i == 0){ pTok = &(valStr[0]);}
      else      { pTok = &(valStr[strlen(valStr)]);}

      if((fabs(x-px) < eps) && (fabs(y-py) < eps))
      {
         //no-op
      }
      else if( pCur->pNxt != NULL)
      { 
         sprintf(pTok, "%.6lf  %.6lf  %.6lf\n",pCur->x, pCur->y, pCur->z);
      }
      else //no CR/LF on last line of input
      {
         sprintf(pTok, "%.6lf  %.6lf  %.6lf",pCur->x, pCur->y, pCur->z);
      }
      px = x;
      py = y;
      i++;
      pCur = pCur->pNxt;
   }
}/* end GetValAsStr() */

/******************************************************************************
GeomParamLine3::GetValStrSize()

Detrmine the string size needed to hold the entire geometry.
******************************************************************************/
int  GeomParamLine3::GetValStrSize(void)
{
   VertexList * pCur;
   int i;

   pCur = m_pFixed;
   i = 0;
   while(pCur != NULL){ i+=60; pCur = pCur->pNxt; }

   return i;
}/* end GetValStrSize() */

/******************************************************************************
GeomParamLine3::InsertVertex()

Insert the vertex into the vertex list.
******************************************************************************/
void GeomParamLine3::InsertVertex(AugVertexList * pNew)
{
   AugVertexList * pCur;

   if(m_pInit == NULL){ m_pInit = pNew;}
   else
   {
      pCur = m_pInit;
      while(pCur->pNxt != NULL){ pCur = pCur->pNxt;}
      pCur->pNxt = pNew;
   }
} /* end InsertVertex() */

/******************************************************************************
GeomParamLine3::Write()

Writes formatted output to pFile argument.
******************************************************************************/
void GeomParamLine3::Write(FILE * pFile, int type)
{
   AugVertexList * pCur;
   VertexList conv;
   
   pCur = m_pInit;
   if (type == WRITE_DBG)
   {
      fprintf(pFile, "Name = %s\n", m_pName);
      while(pCur != NULL)
      {
         ConvertAugVertex(pCur, &conv);
         fprintf(pFile, "x-coord = %lf\n", conv.x);
         fprintf(pFile, "y-coord = %lf\n", conv.y);
         fprintf(pFile, "z-value = %lf\n", conv.z);
         pCur = pCur->pNxt;
      }
   }/* end if() */
} /* end GeomParamLine3::Write() */

/******************************************************************************
GeomParamPoly3::Destroy()
******************************************************************************/
void GeomParamPoly3::Destroy(void)
{
   delete [] m_pName;

   //destroy vertices
   DestroyAugVertexList(m_pInit);
   DestroyVertexList(m_pFixed);

   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR (GeomParamPoly3)
******************************************************************************/
GeomParamPoly3::GeomParamPoly3(void)
{
   m_pInit = NULL;
   m_pFixed = NULL;

   m_pName = NULL;
   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
CTOR (GeomParamPoly3)
******************************************************************************/
GeomParamPoly3::GeomParamPoly3(IroncladString name)
{
   int len;

   m_pInit = NULL;
   m_pFixed = NULL;

   len = (int)strlen(name) + 10;
   NEW_PRINT("char", len);
   m_pName = new char[len];
   MEM_CHECK(m_pName);

   strcpy(m_pName, name);

   IncCtorCount();
} /* end CTOR */

/******************************************************************************
GeomParamPoly3::Convert()

Convert augmented vertices to normal vertices. Augmented vertices may contain
parameter names instead of constant values, such as when geometry of a given
element is being calibrated.
******************************************************************************/
void GeomParamPoly3::Convert(void)
{
   AugVertexList * pCur;
   VertexList * pTmp;

   //first time? then allocated space for fixed list
   if((m_pFixed == NULL) && (m_pInit != NULL))
   {
      NEW_PRINT("VertexList", 1);
      m_pFixed = new VertexList;
      m_pFixed->pNxt = NULL;
      MEM_CHECK(m_pFixed);
   }

   //convert each element in list
   pCur = m_pInit;
   pTmp = m_pFixed;
   while(pCur != NULL)
   {
      ConvertAugVertex(pCur, pTmp);
      pCur = pCur->pNxt;
      if(pCur != NULL)          
      {
         if(pTmp->pNxt == NULL)
         {
            NEW_PRINT("VertexList", 1);
            pTmp->pNxt = new VertexList;         
            pTmp = pTmp->pNxt;
            MEM_CHECK(pTmp);
            pTmp->pNxt = NULL;
         }
         else
         {
            pTmp = pTmp->pNxt;
         }
      }
      else
      {
         DestroyVertexList(pTmp->pNxt);
         pTmp->pNxt = NULL;
      }      
   }/* end while() */   
} /* end Convert() */

/******************************************************************************
GeomParamPoly3::Reorder()

Reorder vertices if they overlap one another. Returns true if reorder is 
successful, false otherwise.
******************************************************************************/
bool GeomParamPoly3::Reorder(void)
{
   Point2D pt;
   Segment2D seg1, seg2;
   AugVertexList * pCurAug;   
   VertexList * pCur, * pTmp;
   double xupr, xlwr, yupr, ylwr, tmpx, tmpy;
   bool mustCheck, backProp, done;
   int count, numTries, maxTries, i, r, j;

   pCurAug = m_pInit;
   pCur = m_pFixed;

   yupr = xupr = -1.00;
   ylwr = xlwr = +1.00;
   //first, determine if reordering check is even needed
   mustCheck = false;
   backProp = true;
   count = 0;
   while(pCurAug != NULL)
   {
      count++;
      //only check shapes with parameters attached
      if((pCurAug->px != NULL) || (pCurAug->py != NULL) || 
         (pCurAug->tx != NULL) || (pCurAug->ty != NULL))
      {
         mustCheck = true;
      }
      /* only back-propagate if all vertices are params w/ same bounds */
      if((pCurAug->px == NULL) || (pCurAug->py == NULL))
      {
         backProp = false;
      }
      else //check if same bounds for all
      {
         //first time?
         if(xupr < xlwr){ xupr = pCurAug->px->GetUprBnd(); xlwr = pCurAug->px->GetLwrBnd();}
         if(yupr < ylwr){ yupr = pCurAug->py->GetUprBnd(); ylwr = pCurAug->py->GetLwrBnd();}

         if((pCurAug->px->GetUprBnd() != xupr) || (pCurAug->px->GetLwrBnd() != xlwr)){ backProp = false;}
         if((pCurAug->py->GetUprBnd() != yupr) || (pCurAug->py->GetLwrBnd() != ylwr)){ backProp = false;}
      }
      pCurAug = pCurAug->pNxt;   
   }/* end while() */

   if(mustCheck == false){ return true;}   
   
   done = false;
   numTries = 0;

   maxTries = 1;
   for(i = count; i > 1; i--){ maxTries *= i;}
   maxTries *= 2;

   while((done == false) && (numTries < maxTries))
   {
      done = true;      
      //compare all segments to see if there is overlap
      for(pCur = m_pFixed; pCur->pNxt != NULL; pCur = pCur->pNxt)
      {
         seg1.p1.x = pCur->x;
         seg1.p1.y = pCur->y;
         seg1.p2.x = pCur->pNxt->x;
         seg1.p2.y = pCur->pNxt->y;

         for(pTmp = pCur->pNxt; pTmp != NULL; pTmp = pTmp->pNxt)
         {
            seg2.p1.x = pTmp->x;
            seg2.p1.y = pTmp->y;
            if(pTmp->pNxt != NULL)
            {
               seg2.p2.x = pTmp->pNxt->x;
               seg2.p2.y = pTmp->pNxt->y;
            }
            else //polygons
            {
               seg2.p2.x = m_pFixed->x;
               seg2.p2.y = m_pFixed->y;
            }

            if(SegIntersect(&seg1, &seg2, &pt) == BOTHSEG)
            {
               done = false;
               break;
            }
         }/* end for() */
         if(done == false){ break;}      
      }/* end for() */

      //randomly reorder, if overlap detected      
      if(done == false)
      {
         pCur = m_pFixed->pNxt;
         for(i = 1; i < (count-1); i++)
         {
            r = MyRand()%(count-i);
            pTmp = m_pFixed->pNxt;
            for(j = 0; j < r; j++){ pTmp = pTmp->pNxt;}
            tmpx = pTmp->x;
            tmpy = pTmp->y;
            pTmp->x = pCur->x;
            pTmp->y = pCur->y;
            pCur->x = tmpx;
            pCur->y = tmpy;
            pCur = pCur->pNxt;
         }/* end for() */
         numTries++;
      }/* end while() */  
   }/* end while() */

   //adjust parameter values to reflect new ordering
   if(backProp == false){ LogError(ERR_MISMATCH, "Can't back-propagate reordering");}
   else //parameters are interchangeable, so back-propagate reordering
   {
      pCurAug = m_pInit;
      pCur = m_pFixed;
      while(pCurAug != NULL)
      {
         pCurAug->px->SetEstVal(pCur->x);
         pCurAug->py->SetEstVal(pCur->y);
         pCurAug = pCurAug->pNxt;
         pCur = pCur->pNxt;
      }/* end while() */
   }/* end else() */

   //if(done == false){ printf("Re-order failed\n");}
   //else{ printf("Re-order succeeded\n");}   
   return done;
}/* end Reorder() */

/******************************************************************************
GeomParamPoly3::FixVertices()

Insert vertices wherever the geometry overlaps with pOther. Returns true if 
vertices are corrected successfully (or if correction was not necessary), 
false otherwise.
******************************************************************************/
bool GeomParamPoly3::FixVertices(GeomParamABC * pOther)
{
   Segment2D seg;
   VertexList * pCur, * pFix, *pTmp;
   pCur = m_pFixed;
   double za, zb, zmid, xa, ya, xb, yb, xc, yc, d1, d2;

   //let the other geometry check each line segment of this geomettry
   while(pCur != NULL)
   {
      xa = seg.p1.x = pCur->x;
      ya = seg.p1.y = pCur->y;
      za = pCur->z;

      if(pCur->pNxt != NULL)
      {
         xb = seg.p2.x = pCur->pNxt->x;
         yb = seg.p2.y = pCur->pNxt->y;
         zb = pCur->pNxt->z;
      }
      else //closing segment of polygon
      {
         xb = seg.p2.x = m_pFixed->x;
         yb = seg.p2.y = m_pFixed->y;
         zb = m_pFixed->z;
      }

      //fix the segment, if necessary
      pFix = pOther->FixVertex(&seg);

      //insert corrections
      if(pFix != NULL)
      { 
         pTmp = pCur->pNxt;
         pCur->pNxt = pFix;
         while(pFix->pNxt != NULL)
         { 
            //interpolate z-values of inserted vertex
            xc = pFix->x;
            yc = pFix->y;
            d1 = sqrt((xc-xa)*(xc-xa)+(yc-ya)*(yc-ya));
            d2 = sqrt((xc-xb)*(xc-xb)+(yc-yb)*(yc-yb));
            zmid = za*(d2/(d1+d2)) + zb*(d1/(d1+d2));
            pFix->z = zmid;

            pFix = pFix->pNxt;
         }
         //interpolate z-values of inserted vertex (distance-weighted average)
         xc = pFix->x;
         yc = pFix->y;
         d1 = sqrt((xc-xa)*(xc-xa)+(yc-ya)*(yc-ya));
         d2 = sqrt((xc-xb)*(xc-xb)+(yc-yb)*(yc-yb));
         zmid = za*(d2/(d1+d2)) + zb*(d1/(d1+d2));
         pFix->z = zmid;

         pFix->pNxt = pTmp;
         //printf("Fixed geometry\n");
      }/* end if() */
      else
      {
         pCur = pCur->pNxt;
      }
   }/* end while() */
     
   return true;
}/* end FixVertices() */

/******************************************************************************
GeomParamPoly3::FixVertex()

Insert vertices wherever the geometry overlaps with the line segment in pSeg. 
Returns a vertex list containing all inserted vertices, or NULL if no inserions
are needed.
******************************************************************************/
VertexList * GeomParamPoly3::FixVertex(Segment2D * pSeg)
{
   Point2D pt;
   Segment2D mySeg;
   VertexList * pCur, * pMyNew, * pRetNew, * pRet, * pTmp;
   double za, zb, zmid, xa, ya, xb, yb, xc, yc, d1, d2;
   int test;

   pCur = m_pFixed;

   pRet = pMyNew = pRetNew = NULL;
   //compare all segments with pSeg to see if there is overlap
   for(pCur = m_pFixed; pCur != NULL; pCur = pCur->pNxt)      
   {
      xa = mySeg.p1.x = pCur->x;
      ya = mySeg.p1.y = pCur->y;
      za = pCur->z;

      if(pCur->pNxt != NULL)
      {
         xb = mySeg.p2.x = pCur->pNxt->x;
         yb = mySeg.p2.y = pCur->pNxt->y;
         zb = pCur->pNxt->z;
      }
      else //closing segment of polygon
      {
         xb = mySeg.p2.x = m_pFixed->x;
         yb = mySeg.p2.y = m_pFixed->y;
         zb =  m_pFixed->z;
      }

      test = SegIntersect(&mySeg, pSeg, &pt);
      if(test > NO_SEGS)
      {
         if((test == BOTHSEG) || (test == LEFTSEG))
         {
            NEW_PRINT("VertexList", 1);
            pMyNew  = new VertexList;         
            MEM_CHECK(pMyNew);

            xc = pMyNew->x = pt.x;
            yc = pMyNew->y = pt.y;
            d1 = sqrt((xc-xa)*(xc-xa)+(yc-ya)*(yc-ya));
            d2 = sqrt((xc-xb)*(xc-xb)+(yc-yb)*(yc-yb));
            zmid = za*(d2/(d1+d2)) + zb*(d1/(d1+d2));
            pMyNew->z = zmid;

            //insert the new vertex
            pTmp = pCur->pNxt;
            pCur->pNxt = pMyNew;
            pMyNew->pNxt = pTmp;
         }/* end if() */
         
         if((test == BOTHSEG) || (test == RGHTSEG))
         {
            if(pRet == NULL)
            { 
               NEW_PRINT("VertexList", 1);
               pRet = new VertexList;
               MEM_CHECK(pRet);
               pRetNew = pRet;
               pRetNew->pNxt = NULL;
            }
            else
            {
               NEW_PRINT("VertexList", 1);
               pRetNew->pNxt = new VertexList;
               MEM_CHECK(pRet);

               pRetNew = pRetNew->pNxt;
               pRetNew->pNxt = NULL;
            }/* end else() */
            pRetNew->x = pt.x;
            pRetNew->y = pt.y;
         }/* end if() */

         //found an intersection, so quit checking
         break;   
      }/* end if() */
   }/* end for() */

   return pRet;
}/* end FixVertex() */

/******************************************************************************
GeomParamPoly3::GetValStrSize()

Detrmine the string size needed to hold the entire geometry.
******************************************************************************/
int  GeomParamPoly3::GetValStrSize(void)
{
   VertexList * pCur;
   int i;

   pCur = m_pFixed;
   i = 0;
   while(pCur != NULL){ i+=60; pCur = pCur->pNxt; }

   return i;
}/* end GetValStrSize() */

/******************************************************************************
GeomParamPoly3::GetValAsStr()

Format the fixed geometry for output to model.
******************************************************************************/
void GeomParamPoly3::GetValAsStr(UnmoveableString valStr)
{
   double x, y, px, py; //previous coords
   double eps = 0.000001;
   char * pTok;
   VertexList * pCur;
   int i;

   pCur = m_pFixed;
   px = py = -1.00; //initialize to an unlikely coordinate
   i = 0;
   while(pCur != NULL)
   {
      x = pCur->x;
      y = pCur->y;

      if(i == 0){ pTok = &(valStr[0]);}
      else      { pTok = &(valStr[strlen(valStr)]);}

      if((fabs(x-px) < eps) && (fabs(y-py) < eps))
      {
         //no-op
      }
      else if(pCur->pNxt != NULL)
      {          
         sprintf(pTok, "%.6lf  %.6lf  %.6lf\n",pCur->x, pCur->y, pCur->z);
      }
      else //no CR/LF on last line of input
      {
         sprintf(pTok, "%.6lf  %.6lf  %.6lf",pCur->x, pCur->y, pCur->z);
      }
      py = y;
      px = x;
      i++;
      pCur = pCur->pNxt;
   }
}/* end GetValAsStr() */

/******************************************************************************
GeomParamPoly3::InsertVertex()

Insert the vertex into the vertex list.
******************************************************************************/
void GeomParamPoly3::InsertVertex(AugVertexList * pNew)
{
   AugVertexList * pCur;

   if(m_pInit == NULL){ m_pInit = pNew;}
   else
   {
      pCur = m_pInit;
      while(pCur->pNxt != NULL){ pCur = pCur->pNxt;}
      pCur->pNxt = pNew;
   }
} /* end InsertVertex() */

/******************************************************************************
GeomParamPoly3::Write()

Writes formatted output to pFile argument.
******************************************************************************/
void GeomParamPoly3::Write(FILE * pFile, int type)
{
   AugVertexList * pCur;
   VertexList conv;
   
   pCur = m_pInit;
   if (type == WRITE_DBG)
   {
      fprintf(pFile, "Name = %s\n", m_pName);
      while(pCur != NULL)
      {
         ConvertAugVertex(pCur, &conv);
         fprintf(pFile, "x-coord = %lf\n", conv.x);
         fprintf(pFile, "y-coord = %lf\n", conv.y);
         fprintf(pFile, "z-value = %lf\n", conv.z);
         pCur = pCur->pNxt;
      }
   }/* end if() */
} /* end GeomParamPoly3::Write() */

/******************************************************************************
GeomParamPoly2::Destroy()
******************************************************************************/
void GeomParamPoly2::Destroy(void)
{
   delete [] m_pName;

   //destroy vertices
   DestroyAugVertexList(m_pInit);
   DestroyVertexList(m_pFixed);

   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR (GeomParamPoly2)
******************************************************************************/
GeomParamPoly2::GeomParamPoly2(void)
{
   m_pInit = NULL;
   m_pFixed = NULL;

   m_pName = NULL;
   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
CTOR (GeomParamPoly2)
******************************************************************************/
GeomParamPoly2::GeomParamPoly2(IroncladString name)
{
   int len;

   m_pInit = NULL;
   m_pFixed = NULL;

   len = (int)strlen(name) + 10;
   NEW_PRINT("char", len);
   m_pName = new char[len];
   MEM_CHECK(m_pName);

   strcpy(m_pName, name);

   IncCtorCount();
} /* end CTOR */

/******************************************************************************
GeomParamPoly2::Convert()

Convert augmented vertices to normal vertices. Augmented vertices may contain
parameter names instead of constant values, such as when geometry of a given
element is being calibrated.
******************************************************************************/
void GeomParamPoly2::Convert(void)
{
   AugVertexList * pCur;
   VertexList * pTmp;

   //first time? then allocated space for fixed list
   if((m_pFixed == NULL) && (m_pInit != NULL))
   {
      NEW_PRINT("VertexList", 1);
      m_pFixed = new VertexList;
      m_pFixed->pNxt = NULL;
      MEM_CHECK(m_pFixed);
   }

   //convert each element in list
   pCur = m_pInit;
   pTmp = m_pFixed;
   while(pCur != NULL)
   {
      ConvertAugVertex(pCur, pTmp);
      pCur = pCur->pNxt;
      if(pCur != NULL)          
      {
         if(pTmp->pNxt == NULL)
         {
            NEW_PRINT("VertexList", 1);
            pTmp->pNxt = new VertexList;         
            pTmp = pTmp->pNxt;
            MEM_CHECK(pTmp);
            pTmp->pNxt = NULL;
         }
         else
         {
            pTmp = pTmp->pNxt;
         }
      }
      else
      {
         DestroyVertexList(pTmp->pNxt);
         pTmp->pNxt = NULL;
      }      
   }/* end while() */   
} /* end Convert() */

/******************************************************************************
GeomParamPoly2::Reorder()

Reorder vertices if they overlap one another. Returns true if reorder is 
successful, false otherwise.
******************************************************************************/
bool GeomParamPoly2::Reorder(void)
{
   Point2D pt;
   Segment2D seg1, seg2;
   AugVertexList * pCurAug;   
   VertexList * pCur, * pTmp;
   double xupr, xlwr, yupr, ylwr, tmpx, tmpy;
   bool mustCheck, backProp, done;
   int count, numTries, maxTries, i, r, j;

   pCurAug = m_pInit;
   pCur = m_pFixed;

   yupr = xupr = -1.00;
   ylwr = xlwr = +1.00;
   //first, determine if reordering check is even needed
   mustCheck = false;
   backProp = true;
   count = 0;
   while(pCurAug != NULL)
   {
      count++;
      //only check shapes with parameters attached
      if((pCurAug->px != NULL) || (pCurAug->py != NULL) || 
         (pCurAug->tx != NULL) || (pCurAug->ty != NULL))
      {
         mustCheck = true;
      }
      /* only back-propagate if all vertices are params w/ same bounds */
      if((pCurAug->px == NULL) || (pCurAug->py == NULL))
      {
         backProp = false;
      }
      else //check if same bounds for all
      {
         //first time?
         if(xupr < xlwr){ xupr = pCurAug->px->GetUprBnd(); xlwr = pCurAug->px->GetLwrBnd();}
         if(yupr < ylwr){ yupr = pCurAug->py->GetUprBnd(); ylwr = pCurAug->py->GetLwrBnd();}

         if((pCurAug->px->GetUprBnd() != xupr) || (pCurAug->px->GetLwrBnd() != xlwr)){ backProp = false;}
         if((pCurAug->py->GetUprBnd() != yupr) || (pCurAug->py->GetLwrBnd() != ylwr)){ backProp = false;}
      }
      pCurAug = pCurAug->pNxt;   
   }/* end while() */

   if(mustCheck == false){ return true;}   
   
   done = false;
   numTries = 0;

   maxTries = 1;
   for(i = count; i > 1; i--){ maxTries *= i;}
   maxTries *= 2;
  
   while((done == false) && (numTries < maxTries))
   {
      done = true;      
      //compare all segments to see if there is overlap
      for(pCur = m_pFixed; pCur->pNxt != NULL; pCur = pCur->pNxt)
      {
         seg1.p1.x = pCur->x;
         seg1.p1.y = pCur->y;
         seg1.p2.x = pCur->pNxt->x;
         seg1.p2.y = pCur->pNxt->y;

         for(pTmp = pCur->pNxt; pTmp != NULL; pTmp = pTmp->pNxt)
         {
            seg2.p1.x = pTmp->x;
            seg2.p1.y = pTmp->y;
            if(pTmp->pNxt != NULL)
            {
               seg2.p2.x = pTmp->pNxt->x;
               seg2.p2.y = pTmp->pNxt->y;
            }
            else //polygons only
            {
               seg2.p2.x = m_pFixed->x;
               seg2.p2.y = m_pFixed->y;
            }

            if(SegIntersect(&seg1, &seg2, &pt) == BOTHSEG)
            {
               done = false;
               break;
            }
         }/* end for() */
         if(done == false){ break;}      
      }/* end for() */

      //randomly reorder, if overlap detected      
      if(done == false)
      {
         pCur = m_pFixed->pNxt;
         for(i = 1; i < (count-1); i++)
         {
            r = MyRand()%(count-i);
            pTmp = m_pFixed->pNxt;
            for(j = 0; j < r; j++){ pTmp = pTmp->pNxt;}
            tmpx = pTmp->x;
            tmpy = pTmp->y;
            pTmp->x = pCur->x;
            pTmp->y = pCur->y;
            pCur->x = tmpx;
            pCur->y = tmpy;
            pCur = pCur->pNxt;
         }/* end for() */
         numTries++;
      }/* end if() */  
   }/* end while() */

   //adjust parameter values to reflect new ordering
   if(backProp == false){ LogError(ERR_MISMATCH, "Can't back-propagate reordering");}
   else //parameters are interchangeable, so back-propagate reordering
   {
      pCurAug = m_pInit;
      pCur = m_pFixed;
      while(pCurAug != NULL)
      {
         pCurAug->px->SetEstVal(pCur->x);
         pCurAug->py->SetEstVal(pCur->y);
         pCurAug = pCurAug->pNxt;
         pCur = pCur->pNxt;
      }/* end while() */
   }/* end else() */

   //if(done == false){ printf("Re-order failed\n");}
   //else{ printf("Re-order succeeded\n");}   
   return done;
}/* end Reorder() */

/******************************************************************************
GeomParamPoly2::FixVertices()

Insert vertices wherever the geometry overlaps with pOther. Returns true if 
vertices are corrected successfully (or if correction was not necessary), 
false otherwise.
******************************************************************************/
bool GeomParamPoly2::FixVertices(GeomParamABC * pOther)
{
   Segment2D seg;
   VertexList * pCur, * pFix, * pTmp;
   pCur = m_pFixed;

   //let the other geometry check each line segment of this geomettry
   while(pCur != NULL)
   {
      seg.p1.x = pCur->x;
      seg.p1.y = pCur->y;      

      if(pCur->pNxt != NULL)
      {
         seg.p2.x = pCur->pNxt->x;
         seg.p2.y = pCur->pNxt->y;         
      }
      else //closing segment of polygon
      {
         seg.p2.x = m_pFixed->x;
         seg.p2.y = m_pFixed->y;         
      }

      //fix the segment, if necessary
      pFix = pOther->FixVertex(&seg);

      //insert corrections
      if(pFix != NULL)
      { 
         pTmp = pCur->pNxt;
         pCur->pNxt = pFix;
         while(pFix->pNxt != NULL)
         {                         
            pFix->z = 0.00;
            pFix = pFix->pNxt;
         }
         pFix->z = 0.00;
         pFix->pNxt = pTmp;         

         //printf("Fixed geometry\n");
      }/* end if() */
      else
      {
         pCur = pCur->pNxt;      
      }
   }/* end while() */
   return true;
}/* end FixVertices() */

/******************************************************************************
GeomParamPoly2::FixVertex()

Insert vertices wherever the geometry overlaps with the line segment in pSeg. 
Returns a vertex list containing all inserted vertices, or NULL if no inserions
are needed.
******************************************************************************/
VertexList * GeomParamPoly2::FixVertex(Segment2D * pSeg)
{
   int test;
   Point2D pt;
   Segment2D mySeg;
   VertexList * pCur, * pMyNew, * pRetNew, * pRet, * pTmp;

   pCur = m_pFixed;

   pRet = pMyNew = pRetNew = NULL;
   //compare all segments with pSeg to see if there is overlap
   for(pCur = m_pFixed; pCur != NULL; pCur = pCur->pNxt)      
   {
      mySeg.p1.x = pCur->x;
      mySeg.p1.y = pCur->y;      

      if(pCur->pNxt != NULL)
      {
         mySeg.p2.x = pCur->pNxt->x;
         mySeg.p2.y = pCur->pNxt->y;         
      }
      else //closing segment of polygon
      {
         mySeg.p2.x = m_pFixed->x;
         mySeg.p2.y = m_pFixed->y;         
      }

      test = SegIntersect(&mySeg, pSeg, &pt);
      if(test > NO_SEGS)
      {
         if((test == BOTHSEG) || (test == LEFTSEG))
         {
            NEW_PRINT("VertexList", 1);
            pMyNew  = new VertexList;         
            MEM_CHECK(pMyNew);

            pMyNew->x = pt.x;
            pMyNew->y = pt.y;
            pMyNew->z = 0.00;

            //insert the new vertex
            pTmp = pCur->pNxt;
            pCur->pNxt = pMyNew;
            pMyNew->pNxt = pTmp;
         }/* end if() */

         if((test == BOTHSEG) || (test == RGHTSEG))
         {        
            if(pRet == NULL)
            { 
               NEW_PRINT("VertexList", 1);
               pRet = new VertexList;
               MEM_CHECK(pRet);
               pRetNew = pRet;
               pRetNew->pNxt = NULL;
            }
            else
            {
               NEW_PRINT("VertexList", 1);
               pRetNew->pNxt = new VertexList;
               MEM_CHECK(pRet);

               pRetNew = pRetNew->pNxt;
               pRetNew->pNxt = NULL;
            }/* end else() */
            pRetNew->x = pt.x;
            pRetNew->y = pt.y;            
         }/* end if() */

         //found an intersection, so quit
         break;
      }/* end if() */
   }/* end for() */

   return pRet;
}/* end FixVertex() */

/******************************************************************************
GeomParamPoly2::GetValStrSize()

Detrmine the string size needed to hold the entire geometry.
******************************************************************************/
int  GeomParamPoly2::GetValStrSize(void)
{
   VertexList * pCur;
   int i;

   pCur = m_pFixed;
   i = 0;
   while(pCur != NULL){ i+=40; pCur = pCur->pNxt; }

   return i;
}/* end GetValStrSize() */

/******************************************************************************
GeomParamPoly2::GetValAsStr()

Format the fixed geometry for output to model.
******************************************************************************/
void GeomParamPoly2::GetValAsStr(UnmoveableString valStr)
{
   double x, y, px, py;
   double eps = 0.000001;
   char * pTok;
   VertexList * pCur;
   int i;

   pCur = m_pFixed;
   i = 0;
   px = py = -1.00;
   while(pCur != NULL)
   {
      if(i == 0){ pTok = &(valStr[0]);}
      else      { pTok = &(valStr[strlen(valStr)]);}

      x = pCur->x;
      y = pCur->y;
      if((fabs(x-px) < eps) && (fabs(y-py) < eps))
      {
         //no-op
      }
      else if( pCur->pNxt != NULL)
      { 
         sprintf(pTok, "%.6lf  %.6lf\n",pCur->x, pCur->y);
      }
      else //no CR/LF on last line of input
      {
         sprintf(pTok, "%.6lf  %.6lf",pCur->x, pCur->y);
      }
      px = x;
      py = y;
      i++;
      pCur = pCur->pNxt;
   }
}/* end GetValAsStr() */

/******************************************************************************
GeomParamPoly2::InsertVertex()

Insert the vertex into the vertex list.
******************************************************************************/
void GeomParamPoly2::InsertVertex(AugVertexList * pNew)
{
   AugVertexList * pCur;

   if(m_pInit == NULL){ m_pInit = pNew;}
   else
   {
      pCur = m_pInit;
      while(pCur->pNxt != NULL){ pCur = pCur->pNxt;}
      pCur->pNxt = pNew;
   }
} /* end InsertVertex() */

/******************************************************************************
GeomParamPoly2::Write()

Writes formatted output to pFile argument.
******************************************************************************/
void GeomParamPoly2::Write(FILE * pFile, int type)
{
   AugVertexList * pCur;
   VertexList conv;
   
   pCur = m_pInit;
   if (type == WRITE_DBG)
   {
      fprintf(pFile, "Name = %s\n", m_pName);
      while(pCur != NULL)
      {
         ConvertAugVertex(pCur, &conv);
         fprintf(pFile, "x-coord = %lf\n", conv.x);
         fprintf(pFile, "y-coord = %lf\n", conv.y);
         pCur = pCur->pNxt;
      }
   }/* end if() */
} /* end GeomParamPoly2::Write() */

/******************************************************************************
GeomParamCirc4::Destroy()
******************************************************************************/
void GeomParamCirc4::Destroy(void)
{
   delete [] m_pName;
   delete m_pInit;
   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR (GeomParamCirc4)
******************************************************************************/
GeomParamCirc4::GeomParamCirc4(void)
{
   m_pName = NULL;
   m_pInit = NULL;
   m_Fixed.r = 0.00;
   m_Fixed.x = 0.00;
   m_Fixed.y = 0.00;
   m_Zcur    = 0.00;
   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
CTOR (GeomParamCirc4)
******************************************************************************/
GeomParamCirc4::GeomParamCirc4(IroncladString name, AugCircle * pData)
{
   int len;

   len = (int)strlen(name) + 10;
   NEW_PRINT("char", len);
   m_pName = new char[len];
   MEM_CHECK(m_pName);

   strcpy(m_pName, name);

   m_pInit = pData;
   m_Fixed.r = 0.00;
   m_Fixed.x = 0.00;
   m_Fixed.y = 0.00;
   m_Zcur    = 0.00;

   IncCtorCount();
} /* end CTOR */

/******************************************************************************
GeomParamCirc4::Convert()

Convert augmented circle data to a normal circle. Augmented circles may contain
parameter names instead of constant values, such as when geometry of a given
element is being calibrated.
******************************************************************************/
void GeomParamCirc4::Convert(void)
{
   ConvertAugCircle(m_pInit, &m_Fixed, &m_Zcur);
} /* end Convert() */

/******************************************************************************
GeomParamCirc4::GetVertexList()

Convert (x,y,r) components of the circle into a VertexList structure.
******************************************************************************/
VertexList * GeomParamCirc4::GetVertexList(int * type)
{
   *type = MY_CIRCLE_TYPE;
   return((VertexList *)&m_Fixed);
}/* end GetVertexList() */

/******************************************************************************
GeomParamCirc4::FixVertices()

Reduce circle radius so that it does not intersect with vertices of other 
geometries. Returns true if radius is corrected successfully (or if correction 
was not necessary), false otherwise.
******************************************************************************/
bool GeomParamCirc4::FixVertices(GeomParamABC * pOther)
{
   int type;
   VertexList * pList, * pCur;
   Segment2D seg;
   double rmin, x1, y1, r1, x2, y2, r2, d;
   double eps = 0.000001;

   pList = pOther->GetVertexList(&type);

   //segments tested depends on the pOther type 
   switch(type)
   {
      //must check closing segment of polygons
      case(MY_POLYGON_TYPE) :
      {         
         for(pCur = pList; pCur != NULL; pCur = pCur->pNxt)
         {
            seg.p1.x = pCur->x;
            seg.p1.y = pCur->y;
            if(pCur->pNxt != NULL)
            {
               seg.p2.x = pCur->pNxt->x;
               seg.p2.y = pCur->pNxt->y;
            }
            else
            {
               seg.p2.x = pList->x;
               seg.p2.y = pList->y;
            }
            //check segment-circle intersection
            if(CircSegIntersect(&m_Fixed, &seg, &rmin) == true)
            {
               m_Fixed.r = rmin;
               //back-propagate
               if(m_pInit->pr != NULL){ m_pInit->pr->SetEstVal(m_Fixed.r);}
            }/* end if() */
         }/* end for() */
         break;
      }/* end case() */
      case(MY_LINE_TYPE) :
      {
         for(pCur = pList; pCur->pNxt != NULL; pCur = pCur->pNxt)
         {
            seg.p1.x = pCur->x;
            seg.p1.y = pCur->y;
            seg.p2.x = pCur->pNxt->x;
            seg.p2.y = pCur->pNxt->y;

            //check segment-circle intersection
            if(CircSegIntersect(&m_Fixed, &seg, &rmin) == true)
            {
               m_Fixed.r = rmin;
               //back-propagate
               if(m_pInit->pr != NULL){ m_pInit->pr->SetEstVal(m_Fixed.r);}
            }/* end if() */
         }/* end for() */
         break;
      }/* end case() */
      //compare distance between the two center-points with radius
      case(MY_CIRCLE_TYPE) :
      {
         x1 = m_Fixed.x;
         y1 = m_Fixed.y;
         r1 = m_Fixed.r;
         x2 = pList->x;
         y2 = pList->y;
         r2 = pList->z;
         d = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));

         if(d > (r1+r2+eps)){ return true;}
         m_Fixed.r = MyMax(eps, d - r2 - eps);
         //back-propagate
         if(m_pInit->pr != NULL){ m_pInit->pr->SetEstVal(m_Fixed.r);}
         break;
      }/* end case() */
      default :
      {
         break;
      }
   }/* end switch() */

   return true;
}/* end FixVertices() */

/******************************************************************************
GeomParamCirc4::FixVertex()

Adjust radius if the circle intersects with the given segment
******************************************************************************/
VertexList * GeomParamCirc4::FixVertex(Segment2D * pSeg)
{
   double rmin;

   //check segment-circle intersection
   if(CircSegIntersect(&m_Fixed, pSeg, &rmin) == true)
   {
      m_Fixed.r = rmin;
      //back-propagate
      if(m_pInit->pr != NULL){ m_pInit->pr->SetEstVal(m_Fixed.r);}
   }/* end if() */

   return NULL;
}/* end FixVertex() */

/******************************************************************************
GeomParamCirc4::GetValAsStr()

Format the fixed geometry for output to model.
******************************************************************************/
void GeomParamCirc4::GetValAsStr(UnmoveableString valStr)
{
   sprintf(valStr, "%.6lf  %.6lf  %.6lf  %.6lf",m_Fixed.x, m_Fixed.y, m_Zcur, m_Fixed.r);
}/* end GetValAsStr() */

/******************************************************************************
GeomParamCirc4::Write()

Writes formatted output to pFile argument.
******************************************************************************/
void GeomParamCirc4::Write(FILE * pFile, int type)
{    
   if (type == WRITE_DBG)
   {
      fprintf(pFile, "Name = %s\n", m_pName);
      ConvertAugCircle(m_pInit, &m_Fixed, &m_Zcur);
      fprintf(pFile, "x-ctr  = %lf\n", m_Fixed.x);
      fprintf(pFile, "y-ctr  = %lf\n", m_Fixed.y);
      fprintf(pFile, "z-val  = %lf\n", m_Zcur);
      fprintf(pFile, "radius = %lf\n", m_Fixed.r);
   }/* end if() */
} /* end GeomParamCirc4::Write() */
