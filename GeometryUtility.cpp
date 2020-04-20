/******************************************************************************
File      : GemoetryUtility.cpp
Author    : L. Shawn Matott
Copyright : 2004, L. Shawn Matott

This file contains a bunch of useful c-style computational geometry routines, 
ranging from point-in-polygon determination to point-to-polygon distance 
calculations.

Version History
05-12-04    lsm   created
12-02-04    lsm   Added intersection tests and additional distance calclations
******************************************************************************/
#include <math.h>

#include "GeometryUtility.h"

#include "MyComplex.h"
#include "Utility.h"
#include "Exception.h"

/******************************************************************************
SegIntersect()

Checks to see if two line segments intersect. 

Return value depends on the relationship between the two segments. There are four
possibilities:
   SHARSEG (-1) : segments share a node
   NO_SEGS (0)  : neither segment intersects
   LEFTSEG (1)  : left segment (seg1) intersects a node of the right segment
   RGHTSEG (2)  : right segment (seg2) intersects a node of the left segment
   BOTHSEG (3)  : both segments intersect (true intersection)
******************************************************************************/
int SegIntersect(Segment2D * seg1, Segment2D * seg2, Point2D * pNew)
{
   double m1, m2, b1, b2, x, eps;
   eps = 0.000001;

   //check to see if either segment is actually a point
   if((fabs(seg1->p1.x - seg1->p2.x) < eps) && (fabs(seg1->p1.y - seg1->p2.y) < eps)){ return NO_SEGS;}
   if((fabs(seg2->p1.x - seg2->p2.x) < eps) && (fabs(seg2->p1.y - seg1->p2.y) < eps)){ return NO_SEGS;}

   //check if segments share a node
   if((fabs(seg1->p1.x - seg2->p2.x) < eps) && (fabs(seg1->p1.y - seg2->p2.y) < eps))
   { 
      pNew->x = seg1->p1.x;
      pNew->y = seg1->p1.y;
      return SHARSEG;
   }
   if((fabs(seg1->p1.x - seg2->p1.x) < eps) && (fabs(seg1->p1.y - seg2->p1.y) < eps))
   { 
      pNew->x = seg1->p1.x;
      pNew->y = seg1->p1.y;
      return SHARSEG;
   }
   if((fabs(seg1->p2.x - seg2->p2.x) < eps) && (fabs(seg1->p2.y - seg2->p2.y) < eps))
   { 
      pNew->x = seg1->p2.x;
      pNew->y = seg1->p2.y;
      return SHARSEG;
   }
   if((fabs(seg1->p2.x - seg2->p1.x) < eps) && (fabs(seg1->p2.y - seg2->p1.y) < eps))
   { 
      pNew->x = seg1->p2.x;
      pNew->y = seg1->p2.y;
      return SHARSEG;
   }

   //compute slopes
   m1 = (seg1->p2.y - seg1->p1.y)/(seg1->p2.x - seg1->p1.x);
   m2 = (seg2->p2.y - seg2->p1.y)/(seg2->p2.x - seg2->p1.x);

   //compute intercepts
   b1 = seg1->p2.y - (m1 * seg1->p2.x);
   b2 = seg2->p2.y - (m2 * seg2->p2.x);

   if(m1 == m2) //parallel
   {
      if(b1 == b2)//unlikely, but troublesome....
      { 
         LogError(ERR_MISMATCH, "Parallel geometries with identical intercepts");
      }
      return NO_SEGS;
   }/* end if() */
   else
   {
      pNew->x = x = (b2 - b1)/(m1 - m2);
      pNew->y = m1*x + b1;
      //test for true intersection
      if((x < MyMax(seg1->p1.x, seg1->p2.x) - eps) && 
         (x > MyMin(seg1->p1.x, seg1->p2.x) + eps))
      {
         if((x < MyMax(seg2->p1.x, seg2->p2.x) - eps) && 
            (x > MyMin(seg2->p1.x, seg2->p2.x) + eps))
            { 
               return BOTHSEG;
            }/* end if() */
      }/* end if() */

      //test for intersection between left segement and right node
      if((fabs(seg2->p1.x - x) < eps) || (fabs(seg2->p2.x - x) < eps))
      {
         if((x < MyMax(seg1->p1.x, seg1->p2.x) - eps) && 
            (x > MyMin(seg1->p1.x, seg1->p2.x) + eps))
            { 
               return LEFTSEG;
            }/* end if() */        
      }
      //test for intersection between right segement and left node
      if((fabs(seg1->p1.x - x) < eps) || (fabs(seg1->p2.x - x) < eps))
      {
         if((x < MyMax(seg2->p1.x, seg2->p2.x) - eps) && 
            (x > MyMin(seg2->p1.x, seg2->p2.x) + eps))
            { 
               return RGHTSEG;
            }/* end if() */
      }/* end else() */

      return NO_SEGS;
   }/* end else() */
}/* end SegIntersect() */

/******************************************************************************
PointInPoly()

Checks to see if a point is inside a polygon. Returns true if inside or on 
boundary of the polygon, false otherwise.

Code taken from James Craig via Paul Bourke.
******************************************************************************/
bool PointInPoly(Point2D pt, Point2D * poly, int n)
{
   bool  in = false;
   int   j, i;
   int NLines = n;
   for (i=0, j=NLines-1; i<NLines; j=i++) { 
    if ((((poly[i].y<=pt.y) && (pt.y<poly[j].y)) ||
         ((poly[j].y<=pt.y) && (pt.y<poly[i].y))) &&
          (pt.x<(poly[j].x-poly[i].x)*(pt.y-poly[i].y)/
            (poly[j].y-poly[i].y)+poly[i].x)){
            in=!in;
            }
   }/* end for() */
   return in;
}/* end PointInPoly() */

/******************************************************************************
DistToLine()

Computes the closest distance from a point (p) to a line (a-b).
******************************************************************************/
double DistToLine(Point2D p, Point2D a, Point2D b)
{
   cmp z1(a.x,a.y);
   cmp z2(b.x,b.y);
   cmp z (p.x,p.y);
   
   //if line seg is really a point....
   if(z1 == z2){ return abs(z1-z);}

   cmp Z = (z-0.5*(z1+z2))/(0.5*(z2-z1));

   //test the three possibilities
   if(Z.RE < -1.00)    { return abs(z1-z);} //left of line
   else if(Z.RE > 1.00){ return abs(z2-z);} //right of line
   else                { return fabs(Z.IM)*0.5*abs(z2-z1);} //over line   
}/* end DistToLine() */

/******************************************************************************
DistToSegment()

Computes the closest distance from a point (p) to a line segment (seg).
******************************************************************************/
double DistToSegment(Point2D pt, Segment2D * seg)
{
   double d;
   Point2D a, b;
   a.x = seg->p1.x;
   a.y = seg->p1.y;
   b.x = seg->p2.x;
   b.y = seg->p2.y;
   d = DistToLine(pt, a, b);
   return (d);
}/* end DistToSegment() */

/******************************************************************************
CircSegIntersect()

Determines whether or not the given circle and line segments intersect. Returns
the largest circle radius that would eliminate the intersection.
******************************************************************************/
bool   CircSegIntersect(Circle2D * circ, Segment2D * seg, double * rmin)
{
   double dist, d1, d2, x1, x2, y1, y2, xc, yc;
   Point2D ctr;
   double eps = 0.000001;

   //check to see if the segment is actually a point
   if((fabs(seg->p1.x - seg->p2.x) < eps) && (fabs(seg->p1.y - seg->p2.y) < eps)){ *rmin = circ->r; return false;}

   xc = ctr.x = circ->x;
   yc = ctr.y = circ->y;
   /*
   if the distance to the line segment is farther 
   than the circle radius, then no intersection
   */
   dist = DistToSegment(ctr, seg);
   if(dist > (circ->r + eps))
   {
      *rmin = circ->r;
      //printf("Segment outside circle\n");
      return false;
   }
   /*
   if the distance to each end point is less than 
   the circle radius, the entire line is inside the
   circle and no intersection
   */
   x1 = seg->p1.x;
   y1 = seg->p1.y;
   d1 = sqrt((xc-x1)*(xc-x1) + (yc-y1)*(yc-y1));

   x2 = seg->p2.x;
   y2 = seg->p2.y;
   d2 = sqrt((xc-x2)*(xc-x2) + (yc-y2)*(yc-y2));
   if((d1 < (circ->r - eps)) && (d2 < (circ->r - eps)))
   {
      //printf("Segment within circle\n");
      *rmin = circ->r;
      return false;
   }
   
   /* geometries intersect, compute max 'acceptable' radius */
   //printf("Segment intersects circle\n");
   *rmin = MyMax(eps, dist - eps);
   return true;
}/* end CircSegIntersect() */

/******************************************************************************
DistToPoly()

Computes the closest distance from a point to a polygon.
******************************************************************************/
double DistToPoly(Point2D p, Point2D * poly, int n)
{
   int i, nxt;
   double dmin, tmp;

   dmin = DistToLine(p, poly[0], poly[1]);
   //compute the closest distance for each segment and select min of these
   for(i = 0; i < n; i++)
   {
      nxt = (i+1)%n;
      tmp = DistToLine(p, poly[i], poly[nxt]);
      if(tmp < dmin){ dmin = tmp;}
   }/* end for() */

   //printf("Particle (%lf,%lf) is %lf from plume\n", p.x,p.y,dmin);
   return dmin;
}/* end DistToPoly() */
