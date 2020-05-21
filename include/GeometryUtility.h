/******************************************************************************
File      : GemoetryUtility.h
Author    : L. Shawn Matott
Copyright : 2004, L. Shawn Matott

This file contains a bunch of useful c-style computational geometry routines, 
ranging from point-in-polygon determination to point-to-polygon distance 
calculations.

Version History
05-12-04    lsm   created
12-02-04    lsm   Added intersection tests and additional distance calclations
******************************************************************************/
#ifndef GEOMETRY_UTILITY_H
#define GEOMETRY_UTILITY_H

#define SHARSEG (-1) //segments share a node
#define NO_SEGS (0)  //neither segment intersects
#define LEFTSEG (1)  //left segment intersects a node of the left segment
#define RGHTSEG (2)  //right segment intersects a node of the left segment
#define BOTHSEG (3)  //both segments intersect (true)

extern "C" {
   typedef struct POINT_2D_STRUCT{double x,y;}Point2D;
   typedef struct SEGMENT_2D_STRUCT{Point2D p1,p2;}Segment2D;
   typedef struct CIRCLE_2D_STRUCT{double x,y,r;}Circle2D;

   bool   PointInPoly(Point2D pt, Point2D * poly, int n);
   int    SegIntersect(Segment2D * seg1, Segment2D * seg2, Point2D * pNew);
   bool   CircSegIntersect(Circle2D * circ, Segment2D * seg, double * rmin);
   double DistToPoly(Point2D p, Point2D * poly, int n);
   double DistToLine(Point2D p, Point2D a, Point2D b);
   double DistToSeg(Point2D pt, Segment2D * seg);
}
#endif /* GEOMETRY_UTILITY_H */

