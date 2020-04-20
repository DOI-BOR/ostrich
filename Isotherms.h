/******************************************************************************
File     : Isotherms.h
Author   : L. Shawn Matott
Copyright: 2004, L. Shawn Matott

Isotherm classes are used to compute isotherms using the parameters and output 
concentrations specified in the IsothermIn.txt input file.

Version History
03-16-04    lsm   Created
03-04-05    lsm   Added PolanyiPartition isotherm
06-09-05    lsm   Added LangmuirPartition isotherm
01-01-07    lsm   Added additional isotherms, the current list of options are:
                     1.  BET Isotherm
                     2.  Freundlich Isotherm
                     3.  Freundlich-Partition Isotherm
                     4.  Linear Isotherm
                     5.  Langmuir Isotherm
                     6.  Generalized Langmuir-Freundlich Isotherm
                     7.  Langmuir-Partition Isotherm
                     8.  Polanyi Isotherm
                     9.  Polanyi-Partition Isotherm
                     10. Toth Isotherm
06-22-07    lsm      11. Dual-Langmuir Isotherm
07-30-07    lsm      12. Orear Isotherm (for testing purposes only)
                     13. McCammon Isotherm (for testing purposes only)
******************************************************************************/
#ifndef ISOTHERMS_H
#define ISOTHERMS_H

#include "MyHeaderInc.h"

//forward decs
class ParameterGroup;
class ObservationGroup;

/******************************************************************************
class IsothermABC

Abstract base class (ABC) defintion for Isotherms. Defines interface used by
concrete instances of isotherms.
********************************************************************************/
class IsothermABC
{ 
   public :
      virtual ~IsothermABC(void){ DBG_PRINT("IsothermABC::DTOR"); }
      virtual void Destroy(void)    = 0;
      virtual bool Initialize(char * pStr) = 0;
      virtual double q(double c)    = 0;
      virtual double dqdc(double c) = 0;
      virtual void Compute(void)    = 0;
      virtual double * GetPtrToC(int * n) = 0;
      virtual char * GetPtrToOutFile(void) = 0;
      virtual bool Initialize(ParameterGroup * pgroup) = 0;
      virtual void Compute(ObservationGroup * ogroup) = 0;
   private:
}; /* end class IsothermABC */

/******************************************************************************
class LinearIsotherm

Implements a Linear Isotherm (q = K C).
********************************************************************************/
class LinearIsotherm : public IsothermABC
{ 
   public :
      LinearIsotherm(void){m_NumOut=0; m_pC=m_pq=NULL; m_Kd=0.00;}
      ~LinearIsotherm(void){ DBG_PRINT("LinearIsotherm::DTOR"); Destroy(); }

      bool Initialize(char * pStr);
      void Destroy(void);
      double q(double c) {return (m_Kd * c);}
      double dqdc(double c) {return m_Kd;}
      void Compute(void);
      double * GetPtrToC(int * n) { *n = m_NumOut; return m_pC;}
      char * GetPtrToOutFile(void) { return m_OutFile;}
      bool Initialize(ParameterGroup * pgroup);
      void Compute(ObservationGroup * ogroup);

   private:
      int    m_NumOut; //number of C and q outputs
      double * m_pC;   //array of concentrations
      double * m_pq;   //array of q
      double m_Kd;     //sorption distribution coeff.
      char   m_OutFile[DEF_STR_SZ];      
}; /* end class LinearIsotherm */

/******************************************************************************
class LangmuirIsotherm

Implements a Langmuir Isotherm (q = (Q0*b C)/(1 + b C)).
********************************************************************************/
class LangmuirIsotherm : public IsothermABC
{ 
   public :
      LangmuirIsotherm(void){m_NumOut=0; m_pC=m_pq=NULL; m_Q0=m_b=0.00;}
      ~LangmuirIsotherm(void){ DBG_PRINT("LangmuirIsotherm::DTOR"); Destroy(); }

      bool Initialize(char * pStr);
      void Destroy(void);
      double q(double c) { return ((m_Q0*m_b*c)/(1 + (m_b*c)));}
      double dqdc(double c) {return ((m_Q0*m_b)/((1 + (m_b*c))*(1 + (m_b*c))));}   
      void Compute(void);
      double * GetPtrToC(int * n) { *n = m_NumOut; return m_pC;}
      char * GetPtrToOutFile(void) { return m_OutFile;}
      bool Initialize(ParameterGroup * pgroup);
      void Compute(ObservationGroup * ogroup);

   private:
      int    m_NumOut; //number of C and q outputs
      double * m_pC;   //array of concentrations
      double * m_pq;   //array of q
      double m_Q0;     //max. sorption density
      double m_b;      //affinity of adsorbent to the adsorbate
      char   m_OutFile[DEF_STR_SZ];      
}; /* end class LangmuirIsotherm */

/******************************************************************************
class DualLangmuirIsotherm

Implements a Dual-Langmuir Isotherm (q = (Q01*b1 C)/(1 + b1 C) + (Q02*b2 C)/(1 + b2 C)).
********************************************************************************/
class DualLangmuirIsotherm : public IsothermABC
{ 
   public :
      DualLangmuirIsotherm(void){m_NumOut=0; m_pC=m_pq=NULL; m_Q01=m_b1=m_Q02=m_b2=0.00;}
      ~DualLangmuirIsotherm(void){ DBG_PRINT("DualLangmuirIsotherm::DTOR"); Destroy(); }

      bool Initialize(char * pStr);
      void Destroy(void);
      double q(double c) { return ( ((m_Q01*m_b1*c)/(1 + (m_b1*c))) + 
                                    ((m_Q02*m_b2*c)/(1 + (m_b2*c))));}

      double dqdc(double c) {return (((m_Q01*m_b1)/((1 + (m_b1*c))*(1 + (m_b1*c)))) +
                                     ((m_Q02*m_b2)/((1 + (m_b2*c))*(1 + (m_b2*c)))));}   

      void Compute(void);
      double * GetPtrToC(int * n) { *n = m_NumOut; return m_pC;}
      char * GetPtrToOutFile(void) { return m_OutFile;}
      bool Initialize(ParameterGroup * pgroup);
      void Compute(ObservationGroup * ogroup);

   private:
      int    m_NumOut; //number of C and q outputs
      double * m_pC;   //array of concentrations
      double * m_pq;   //array of q
      //first sorption type
      double m_Q01;     //max. sorption density
      double m_b1;      //affinity of adsorbent to the adsorbate
      //second sorption type
      double m_Q02;     //max. sorption density
      double m_b2;      //affinity of adsorbent to the adsorbate

      char   m_OutFile[DEF_STR_SZ];      
}; /* end class DualLangmuirIsotherm */

/******************************************************************************
class FreundlichIsotherm

Implements a Freundlich Isotherm (q = Kf C^Nf).
********************************************************************************/
class FreundlichIsotherm : public IsothermABC
{ 
   public :
      FreundlichIsotherm(void){m_NumOut=0; m_pC=m_pq=NULL; m_Kf=m_Nf=0.00;}
      ~FreundlichIsotherm(void){ DBG_PRINT("FreundlichIsotherm::DTOR"); Destroy(); }

      bool Initialize(char * pStr);
      void Destroy(void);
      double q(double c) { return (m_Kf*pow(c, m_Nf));}
      double dqdc(double c) {return (m_Kf*m_Nf*pow(c, (m_Nf-1.00)));}
      void Compute(void);
      double * GetPtrToC(int * n) { *n = m_NumOut; return m_pC;}
      char * GetPtrToOutFile(void) { return m_OutFile;}
      bool Initialize(ParameterGroup * pgroup);
      void Compute(ObservationGroup * ogroup);

   private:
      int    m_NumOut; //number of C and q outputs
      double * m_pC;   //array of concentrations
      double * m_pq;   //array of q
      double m_Kf;     
      double m_Nf;   //(1/n)
      char   m_OutFile[DEF_STR_SZ];      
}; /* end class FreundlichIsotherm */

/******************************************************************************
class PolanyiPartitionIsotherm

Implements a Polanyi-Partition Isotherm (q = Q0 * 10^(-a[log(Sw/C)]^b) + Kp C).
********************************************************************************/
class PolanyiPartitionIsotherm : public IsothermABC
{ 
   public :
      PolanyiPartitionIsotherm(void){m_NumOut=0; m_pC=m_pq=NULL; m_Q0=m_a=m_b=m_Sw=m_Kp=0.00;}
      ~PolanyiPartitionIsotherm(void){ DBG_PRINT("PolanyiPartitionIsotherm::DTOR"); Destroy(); }

      bool Initialize(char * pStr);
      void Destroy(void);
      double q(double c) { return (m_Q0*(pow(10.00, -m_a*pow(log10(m_Sw/c), m_b))) + m_Kp*c);}

      double dqdc(double c) {return ((((-m_Q0*-m_a*m_b)/c)*
                                    (pow(10.00, -m_a*pow(log10(m_Sw/c), m_b)))*
                                    (pow(log10(m_Sw/c), (m_b-1.00)))) + m_Kp);}

      void Compute(void);
      double * GetPtrToC(int * n) { *n = m_NumOut; return m_pC;}
      char * GetPtrToOutFile(void) { return m_OutFile;}
     bool Initialize(ParameterGroup * pgroup);
      void Compute(ObservationGroup * ogroup);

   private:
      int    m_NumOut; //number of C and q outputs
      double * m_pC;   //array of concentrations
      double * m_pq;   //array of q
      double m_Kp; //partition coefficient
      double m_Sw; //Solubility
      double m_Q0; //capacity parameter
      double m_a, m_b; //fitting parameters
      char   m_OutFile[DEF_STR_SZ];      
}; /* end class PolanyiPartitionIsotherm */

/******************************************************************************
class LangmuirPartitionIsotherm

Implements a Langmuir isotherm with partitioning:
   q = [(Q0 * b * C)/(1.00 + b * C)] + Kp * C
********************************************************************************/
class LangmuirPartitionIsotherm : public IsothermABC
{ 
   public :
      LangmuirPartitionIsotherm(void){m_NumOut=0; m_pC=m_pq=NULL; m_Q0=m_b=m_Kp=0.00;}
      ~LangmuirPartitionIsotherm(void){ DBG_PRINT("LangmuirPartitionIsotherm::DTOR"); Destroy(); }

      bool Initialize(char * pStr);
      void Destroy(void);
      double q(double c) { return (((m_Q0 * m_b * c)/(1.00 + m_b * c)) + m_Kp * c);}
      double dqdc(double c) {return (((m_Q0*m_b)/((1 + (m_b*c))*(1 + (m_b*c)))) + m_Kp);}      
      void Compute(void);
      double * GetPtrToC(int * n) { *n = m_NumOut; return m_pC;}
      char * GetPtrToOutFile(void) { return m_OutFile;}
      bool Initialize(ParameterGroup * pgroup);
      void Compute(ObservationGroup * ogroup);

   private:
      int    m_NumOut; //number of C and q outputs
      double * m_pC;   //array of concentrations
      double * m_pq;   //array of q
      double m_Kp; //partition coefficient
      double m_Q0; //capacity parameter
      double m_b; //fitting parameter
      char   m_OutFile[DEF_STR_SZ];      
}; /* end class LangmuirPartitionIsotherm */

/******************************************************************************
class BET_Isotherm

Implements a BET isotherm:
   q = [(Q0*b * C)/[(Sw - C)(1.00 + [(b-1.00)(C/Sw)])]]
********************************************************************************/
class BET_Isotherm : public IsothermABC
{ 
   public :
      BET_Isotherm(void){m_NumOut=0; m_pC=m_pq=NULL; m_Sw=m_Q0=m_b=0.00;}
      ~BET_Isotherm(void){ DBG_PRINT("BET_Isotherm::DTOR"); Destroy(); }

      bool Initialize(char * pStr);
      void Destroy(void);
      double q(double c) { return ((m_Q0 * m_b * c)/((m_Sw - c)*(1.00 + ((m_b - 1.00) * (c / m_Sw)))));}

      double dqdc(double c) { return (((m_Q0 * m_b)*(m_Sw+((m_b - 1.00)*(c*c/m_Sw))))/
                                  (((m_Sw - c)*(m_Sw - c))*(1.00+((m_b - 1.00)*(c/m_Sw)))*
                                  (1.00+((m_b - 1.00)*(c/m_Sw)))));}
      void Compute(void);
      double * GetPtrToC(int * n) { *n = m_NumOut; return m_pC;}
      char * GetPtrToOutFile(void) { return m_OutFile;}
      bool Initialize(ParameterGroup * pgroup);
      void Compute(ObservationGroup * ogroup);

   private:
      int    m_NumOut; //number of C and q outputs
      double * m_pC;   //array of concentrations
      double * m_pq;   //array of q
      double m_Sw; //solubility constant 
      double m_Q0; //capacity parameter
      double m_b; //fitting parameter
      char   m_OutFile[DEF_STR_SZ];      
}; /* end class BET_Isotherm */

/******************************************************************************
class TothIsotherm

Implements a Toth isotherm:
   q = [(Q0 * b * C)/[1.00 + (b*C)^nT]^(1/nT)]
********************************************************************************/
class TothIsotherm : public IsothermABC
{ 
   public :
      TothIsotherm(void){m_NumOut=0; m_pC=m_pq=NULL; m_nT=m_Q0=m_b=0.00;}
      ~TothIsotherm(void){ DBG_PRINT("TothIsotherm::DTOR"); Destroy(); }

      bool Initialize(char * pStr);
      void Destroy(void);
      double q(double c) { return ((m_Q0 * m_b * c)/pow((1.00 + pow(m_b * c, m_nT)),(1.00 / m_nT)));}

      double dqdc(double c) { return (((m_Q0 * m_b)/pow((1.00 + pow(m_b * c, m_nT)),(1.00 / m_nT)))*
                                     (1.00- (pow((m_b*c), m_nT)/ (1.00+ pow((m_b*c), m_nT)))));}

      void Compute(void);
      double * GetPtrToC(int * n) { *n = m_NumOut; return m_pC;}
      char * GetPtrToOutFile(void) { return m_OutFile;}
     bool Initialize(ParameterGroup * pgroup);
      void Compute(ObservationGroup * ogroup);

   private:
      int    m_NumOut; //number of C and q outputs
      double * m_pC;   //array of concentrations
      double * m_pq;   //array of q
      double m_nT; //exponent
      double m_Q0; //capacity parameter
      double m_b; //fitting parameter
      char   m_OutFile[DEF_STR_SZ];      
}; /* end class TothIsotherm */

/******************************************************************************
class LangmuirFreundlichIsotherm

Implements a general Langmuir-Freundlich isotherm:
   q = [(Q0 * (b * C)^nG)/[1.00 + (b*C)^nG]
********************************************************************************/
class LangmuirFreundlichIsotherm : public IsothermABC
{ 
   public :
      LangmuirFreundlichIsotherm(void){m_NumOut=0; m_pC=m_pq=NULL; m_nG=m_Q0=m_b=0.00;}
      ~LangmuirFreundlichIsotherm(void){ DBG_PRINT("LangmuirFreundlichIsotherm::DTOR"); Destroy(); }

      bool Initialize(char * pStr);
      void Destroy(void);
      double q(double c) { return ((m_Q0 * pow((m_b * c), m_nG)) / (1.00 + pow((m_b * c), m_nG)));}

      double dqdc(double c) { return ((m_nG * m_Q0 * m_b* pow((m_b * c), (m_nG - 1.00))) / 
                                     ((1.00 + pow((m_b * c), m_nG)) * (1.00 + pow((m_b * c), m_nG))));}

      void Compute(void);
      double * GetPtrToC(int * n) { *n = m_NumOut; return m_pC;}
      char * GetPtrToOutFile(void) { return m_OutFile;}
      bool Initialize(ParameterGroup * pgroup);
      void Compute(ObservationGroup * ogroup);

   private:
      int    m_NumOut; //number of C and q outputs
      double * m_pC;   //array of concentrations
      double * m_pq;   //array of q
      double m_nG; //exponent
      double m_Q0; //capacity parameter
      double m_b; //fitting parameter
      char   m_OutFile[DEF_STR_SZ];      
}; /* end class LangmuirFreundlichIsotherm */

/******************************************************************************
class PolanyiIsotherm

Implements a Polanyi isotherm:
 q = Q0 * 10^(-a[log(Sw/C)]^b)
********************************************************************************/
class PolanyiIsotherm : public IsothermABC
{ 
   public :
      PolanyiIsotherm(void){m_NumOut=0; m_pC=m_pq=NULL; m_a=m_b=m_Sw=m_Q0=0.00;}
      ~PolanyiIsotherm(void){ DBG_PRINT("PolanyiIsotherm::DTOR"); Destroy(); }

      bool Initialize(char * pStr);
      void Destroy(void);
      double q(double c) { return (m_Q0*(pow(10.00, -m_a*pow(log10(m_Sw/c), m_b))));}

      double dqdc(double c) {return (((-m_Q0*-m_a*m_b)/c)*
                                    (pow(10.00, -m_a*pow(log10(m_Sw/c), m_b)))*
                                    (pow(log10(m_Sw/c), (m_b-1.00))));}

      void Compute(void);
      double * GetPtrToC(int * n) { *n = m_NumOut; return m_pC;}
      char * GetPtrToOutFile(void) { return m_OutFile;}
      bool Initialize(ParameterGroup * pgroup);
      void Compute(ObservationGroup * ogroup);

   private:
      int    m_NumOut; //number of C and q outputs
      double * m_pC;   //array of concentrations
      double * m_pq;   //array of q
      double m_Sw; //solubility
      double m_a; //exponent
      double m_Q0; //capacity parameter
      double m_b; //fitting parameter
      char   m_OutFile[DEF_STR_SZ];      
}; /* end class PolanyiIsotherm */

/******************************************************************************
class FreundlichPartitionIsotherm

Implements a freundlich isotherm with partitioning:
   q = Kf * C^nf + Kp * C
********************************************************************************/
class FreundlichPartitionIsotherm : public IsothermABC
{ 
   public :
      FreundlichPartitionIsotherm(void){m_NumOut=0; m_pC=m_pq=NULL; m_Kf=m_Nf=m_Kp=0.00;}
      ~FreundlichPartitionIsotherm(void){ DBG_PRINT("FreundlichPartitionIsotherm::DTOR"); Destroy(); }

      bool Initialize(char * pStr);
      void Destroy(void);
      double q(double c) { return (m_Kf*pow(c, m_Nf) + m_Kp*c);}
      double dqdc(double c) {return ((m_Kf*m_Nf*pow(c, (m_Nf-1.00))) + m_Kp);}
      void Compute(void);
      double * GetPtrToC(int * n) { *n = m_NumOut; return m_pC;}
      char * GetPtrToOutFile(void) { return m_OutFile;}
      bool Initialize(ParameterGroup * pgroup);
      void Compute(ObservationGroup * ogroup);

   private:
      int    m_NumOut; //number of C and q outputs
      double * m_pC;   //array of concentrations
      double * m_pq;   //array of q
      double m_Kp; //partition coefficient
      double m_Kf; //capacity parameter
      double m_Nf; //exponent
      char   m_OutFile[DEF_STR_SZ];      
}; /* end class FreundlichPartitionIsotherm */

/******************************************************************************
class OrearIsotherm

Implements an Orear "isotherm" (q = a*C - b/C).
********************************************************************************/
class OrearIsotherm : public IsothermABC
{ 
   public :
      OrearIsotherm(void){m_NumOut=0; m_pC=m_pq=NULL; m_a=m_b=0.00;}
      ~OrearIsotherm(void){ DBG_PRINT("OrearIsotherm::DTOR"); Destroy(); }

      bool Initialize(char * pStr);
      void Destroy(void);
      double q(double c) {return ((m_a*c) - (m_b/c));}
      double dqdc(double c) {return (m_a + (m_b/(c*c)));}
      void Compute(void);
      double * GetPtrToC(int * n) { *n = m_NumOut; return m_pC;}
      char * GetPtrToOutFile(void) { return m_OutFile;}
      bool Initialize(ParameterGroup * pgroup){ return true;}
      void Compute(ObservationGroup * ogroup){ return;}

   private:
      int    m_NumOut; //number of C and q outputs
      double * m_pC;   //array of concentrations
      double * m_pq;   //array of q
      double m_a;
      double m_b;
      char   m_OutFile[DEF_STR_SZ];      
}; /* end class OrearIsotherm */

/******************************************************************************
class McCammonIsotherm

Implements a McCammon "isotherm" (A*q^2 +B*q*c + C*c*c + D*q + E*c + F = 0).
********************************************************************************/
class McCammonIsotherm : public IsothermABC
{ 
   public :
      McCammonIsotherm(void){m_NumOut=0; m_pC=m_pq=NULL; m_D=-0.20; 
                             m_F=1.3534; m_A=m_B=m_C=m_E=0.00;}

      ~McCammonIsotherm(void){ DBG_PRINT("McCammonIsotherm::DTOR"); Destroy(); }

      bool Initialize(char * pStr);
      void Destroy(void);
      double q(double c);
      double dqdc(double c);
      void Compute(void);
      double * GetPtrToC(int * n) { *n = m_NumOut; return m_pC;}
      char * GetPtrToOutFile(void) { return m_OutFile;}
      bool Initialize(ParameterGroup * pgroup){ return true;}
      void Compute(ObservationGroup * ogroup){ return;}

   private:
      int    m_NumOut; //number of C and q outputs
      double * m_pC;   //array of concentrations
      double * m_pq;   //array of q
      double m_A;
      double m_B;
      double m_C;
      double m_D;
      double m_E;
      double m_F;
      char   m_OutFile[DEF_STR_SZ];      
}; /* end class OrearIsotherm */

#endif /* ISOTHERMS_H */
