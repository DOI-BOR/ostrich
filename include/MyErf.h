/******************************************************************************
File      : MyErf.h
Author    : L. Shawn Matott
Copyright : 2002, L. Shawn Matott

Contains a function declaration for the complimentray error function and 
a macro implementation for the error function.

Version History
11-18-02    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
******************************************************************************/
#ifndef MY_ERF_H
#define MY_ERF_H

double MyErfc(double x);

#define MyErf(x) (1 - MyErfc(x))

#endif /* MY_ERF_H */

