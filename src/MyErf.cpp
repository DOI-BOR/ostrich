/******************************************************************************
File      : MyErf.cpp
Author    : L. Shawn Matott
Copyright : 2002, L. Shawn Matott

Contains a function definition for the complimentray error function.

Version History
11-18-02    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
******************************************************************************/
#include <math.h>

#include "MyErf.h"

#include "MyTypes.h"

/******************************************************************************
MyErfc()

Complimentary Error Function, courtesy of James Craig.

This only computes positive quadrant of the erfc and takes advantage of 
symmetry of erfc to compute values in the negative quadrant.
******************************************************************************/
double MyErfc(double x)
{
	double tmp(fabs(x)); //take abs so that we are always in positive quadrant.
	double fun;
	double f1;
	double tmp2;
	double tmp3;

	/*
	if |x| > 7, just set at asymptotic limit if you need 
	less than 12 digits of precision
	if(tmp > 7.0)
	{
		fun = 0.0;
	}
	*/
	//less precise approximation is needed for |x| > 3
	if(tmp > 3.0)
	{
		/*
		Traditional Taylor series polynomial expansion of the error function
		f = 1 - 1/(2|x|^2) + 3/(4|x|^4) -5/(6|x|^6)
		*/
		f1  = (1.0 - 1.0/(2.0 * tmp * tmp) 
			       + 3.0/(4.0 * pow(tmp,4)) 
			       - 5.0/(6.0 * pow(tmp,6)));

		fun = f1 * exp(-tmp * tmp) / (tmp * sqrt(MY_PI));
	} /* end if() */
	/*
	Need good approximation in the range |x| < 3, this is where 
	all the fun stuff happens.
	*/
	else
	{
		tmp2 = 1.0 / (1.0 + (0.3275911 * tmp));

		//5th order polynomial interpolation
		tmp3 =   0.254829592  * tmp2 
			   - (0.284496736 * tmp2 * tmp2) 
			   + (1.421413741 * pow(tmp2,3))
			   - (1.453152027 * pow(tmp2,4))
			   + (1.061405429 * pow(tmp2,5));

		fun = tmp3 * exp(-tmp * tmp);
	} /* end else() */

	//convert result according to quadrant.
	if (tmp == x) 
	{
		//pos. quadrant, no change needed
		return fun;
	}
	else
	{
		//neg. quadrant, need to invert
		return (2-fun);
	}
} /* end MyErfc() */


