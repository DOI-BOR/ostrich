/******************************************************************************
File     : MyComlpex.h
Author   : James Craig
Copyright: 2004, James Craig

A complex number class, courtesy James Craig.

Version History
08-20-04    lsm   added to project
01-28-05    lsm   added support for IRIX/GMAKE environment
******************************************************************************/
#ifndef MY_COMPLEX_H
#define MY_COMPLEX_H

#include "MyHeaderInc.h"

/************************************************************************
 *  Class cmp
 *  Complex number Data Abstraction
 ***********************************************************************/
class cmp 
{
   public:
   	double RE;
      double IM;

      cmp()                                 {RE =        IM = 0.00;  IncCtorCount();}
      cmp(const cmp &in)                    {RE = in.RE; IM = in.IM; IncCtorCount();}
      cmp(const double &d1,const double &d2){RE = d1;    IM = d2;    IncCtorCount();}
      cmp(const double &d1)                 {RE = d1;    IM = 0.00;  IncCtorCount();}
      ~cmp(){IncDtorCount();}

      double real(void) {return RE;}
      double imag(void) {return IM;}
		cmp& operator=(const cmp &z)   {RE=z.RE;  IM=z.IM; return *this;}
		cmp& operator=(const double &d){RE=d;     IM=0.0;  return *this;}
};/* end class cmp */

inline double abs(const cmp &z){ return pow(z.RE*z.RE + z.IM*z.IM,0.5);}
inline cmp conj(const cmp &z)  {	return cmp(z.RE,-z.IM);}
inline double arg(const cmp &z)
{
   if ((z.RE!=0.0) || (z.IM!=0.0)) {return atan2(z.IM,z.RE);}
	else                            {return 0.0;}
}
 
inline cmp log(const cmp &z){	return cmp(log(abs(z)),arg(z));}

//addition-----------------------------------------------------
inline cmp operator+(const cmp &z1, const cmp &z2){
  return cmp(z1.RE + z2.RE,z1.IM + z2.IM);
}
inline cmp operator+(const double &d1,const cmp &z2){
  return cmp(z2.RE + d1, z2.IM);
}
inline cmp operator+(const cmp &z1, const double &d2){
  return cmp(z1.RE + d2, z1.IM);
}
inline void operator+=(cmp &z1, const cmp &z2){
  z1.RE=z1.RE+z2.RE;
  z1.IM=z1.IM+z2.IM;
}
inline void operator+=(cmp &z1, const double &d2){
  z1.RE=z1.RE+d2;
}

//subtraction-----------------------------------------------------
inline cmp operator-(const cmp &z1, const cmp &z2){
  return cmp(z1.RE - z2.RE,z1.IM - z2.IM);
}
inline cmp operator-(const double &d1, const cmp &z2){
  return cmp(d1 - z2.RE, -z2.IM);
}
inline cmp operator-(const cmp &z1, const double &d2){
  return cmp(z1.RE - d2, z1.IM);
}
inline void operator-=(cmp &z1,const cmp &z2){
  z1.RE=z1.RE-z2.RE;
  z1.IM=z1.IM-z2.IM;
}
inline void operator-=(cmp &z1, const double &d2){
  z1.RE=z1.RE-d2;
}

//multiplication-----------------------------------------------------
inline cmp operator*(const cmp &z1, const cmp &z2){
  return cmp(z1.RE * z2.RE - z1.IM * z2.IM,
				   z1.RE * z2.IM + z1.IM * z2.RE);
}
inline cmp operator*(const double &d1, const cmp &z2){
  return cmp(d1 * z2.RE, d1 * z2.IM);
}
inline cmp operator*(const cmp &z1, const double &d2){
  return cmp(d2 * z1.RE, d2 * z1.IM);
}
inline void operator*=(cmp &z1, const cmp &z2){
  z1= z1*z2;
}
inline void operator*=(cmp &z1, const double &d2){
  z1.RE=z1.RE*d2;
	z1.IM=z1.IM*d2;
}

//division-----------------------------------------------------
inline cmp operator/(const cmp &z1, const cmp &z2){
	double den=z2.RE*z2.RE+z2.IM*z2.IM;
  return cmp((z1.RE * z2.RE + z1.IM * z2.IM)/den,
				   (z2.RE * z1.IM - z1.RE * z2.IM)/den);
}
inline cmp operator/(const double &d1, const cmp &z2){
  double den=z2.RE*z2.RE+z2.IM*z2.IM;
  return cmp(( d1 * z2.RE )/den,
				   (-d1 * z2.IM )/den);
}
inline cmp operator/(const cmp &z1, const double &d2){
  return cmp(z1.RE / d2, z1.IM / d2);
}
inline void operator/=(cmp &z1, const cmp &z2){
  z1=z1/z2;
}
inline void operator/=(cmp &z1, const double &d2){
  z1.RE=z1.RE/d2;
	z1.IM=z1.IM/d2;
}
//unary operators-----------------------------------------------------
//inline cmp operator+(const cmp &z){
//}
inline cmp operator-(const cmp &z){
	return cmp(-z.RE,-z.IM);
}

//equality operators-----------------------------------------------------
inline bool operator ==(const cmp  &z1,const cmp &z2){
return (z1.RE == z2.RE) && (z1.IM == z2.IM );
}
inline bool operator ==(const cmp  & z1,const double &d2){
return (z1.RE == d2) && (z1.IM == 0.0 );
}
inline bool operator ==(const double  &d1,const cmp &z2){
return (z2.RE == d1) && (z2.IM == 0.0 );
}

//inequality operators-----------------------------------------------------
inline bool operator !=(const cmp  &z1, const cmp  &z2 ){
return (z1.RE != z2.RE) || (z1.IM != z2.IM );
}
inline bool operator !=(const cmp  & z1,const double d2 ){
return (z1.IM != 0.0 ) || (z1.RE != d2);
}
inline bool operator !=(const double  &d1,const cmp &z2){
return (z2.IM != 0.0 ) || (z2.RE != d1);
}

#endif /* MY_COMPLEX_H */
