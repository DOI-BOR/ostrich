/******************************************************************************
File      : TelescopingBounds.h
Author    : L. Shawn Matott
Copyright : 2012, L. Shawn Matott

Encapsulates a group of functions for dynamically adjusting parameter bounds.

Version History
12-12-12    lsm   added copyright information and initial comments.
******************************************************************************/
#ifndef TELESCOPING_BOUNDS_H
#define TELESCOPING_BOUNDS_H

extern "C" {
  double telescope_parameter(double xmin, double xmax, double xbest, double a, double xnew, double (*f)(double, double));
  void revise_bounds(double xmin, double xmax, double xbest, double a, double * xmin_new, double * xmax_new, double (*f)(double, double));
  double flin(double a, double b);
  double fpvx(double a, double b);
  double fdcv(double a, double b);
  double fcve(double a, double b);
  double fvex(double a, double b);
}

#endif
