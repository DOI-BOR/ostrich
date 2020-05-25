/******************************************************************************
File     : BoxCoxModel.h
Author   : L. Shawn Matott
Copyright: 2013, L. Shawn Matott

Computes an objective function value for a given BoxCox transformation. The obj.
function value measures the degree to which the transformation incurs normality
on the transformed and weighted residuals.

Version History
01-09-13    lsm   Created
******************************************************************************/
#ifndef BOX_COX_MODEL_H
#define BOX_COX_MODEL_H

#define BOX_IN_FILE   "BoxCoxIn.txt"
#define BOX_OUT_FILE  "BoxCoxOut.txt"

extern "C" {
int BoxCoxModel(void);
}

#endif /* BOX_COX_MODEL_H */
