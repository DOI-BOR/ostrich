/******************************************************************************
File      : PumpAndTreat.h
Author    : L. Shawn Matott
Copyright : 2004, L. Shawn Matott

Defines a pump-and-treat optimization extension to the ObjectiveFunction class.

This class will support the following Pump-and-Treat objectives:
   Minimze the pumping rate
   Minimize the cost of pumping
   Minimize the cost of installation and pumping
   Minimize the cost of installation, pumping and treatment

Additionally, this class instantiates a set of constraint classes which can be 
added as a penalty to the objective function using a user-defined method 
(additive penalty, multiplicative penalty, etc.). The following constraints are
supported:
   Hydraulic gradient constraint which contain the plume
   Drawdown constraints which limit pumping rates
   Particle capture constraints which ensure plume cature
   Treatment capacity constraints which limit total pumping rate
   
Version History
05-07-04    lsm   created
01-10-05    lsm   Generalized the PatoConstraintABC class and modified to interface 
                  with abstract response variables (RespVarABC)
02-25-05    lsm   Added support for Mayer cost formulation
******************************************************************************/
#ifndef PUMP_AND_TREAT_H
#define PUMP_AND_TREAT_H

#include "MyHeaderInc.h"

// parent class
#include "ObjectiveFunction.h"

// forward decs
class ConstraintABC;
class ParameterGroup;
class RespVarABC;
class GeometryUtility;
class ParameterABC;
class ResponseVarGroup;

/* define a structure to contain plume vertices */
typedef struct PLUME_2D_STRUCT
{
	StringType name; // a name assigned to the plume
	Point2D * poly; //polygon vertices
	int nv; //number of vertices
}Plume2D;

/* define a structure to encapsulate well information */
typedef struct WELL_STRUCT
{
	StringType name;
	//design variables
	ParameterABC * pQ;
	ParameterABC * pXloc;
	ParameterABC * pYloc;
	//head at the well (or very near the well)
	RespVarABC * pHead;
	//surface topography at the well (can be a response variable or a constant)
	RespVarABC * pTopo;
	double Topo;
	//base of aquifer at the (can be a response variable or a constant)
	RespVarABC * pBase;
	double Base;
	//cost information
	double Cdrill, Cpump, Cnrg, Ctot;
}WellStruct;

/* define a structure that maps pump sizes */
typedef struct PUMP_LKUP_TABLE_STRUCT
{
	double Qmin, Qmax, Lmin, Lmax, cost;
}PumpLkupTableStruct;

/******************************************************************************
class PATO (pumo-and-treat optimization)
******************************************************************************/
class PATO : public ObjectiveFunction
{
   public :
      ~PATO(void){ DBG_PRINT("PATO::DTOR"); Destroy(); }
      void Destroy(void);
      PATO(ParameterGroup * pParamGroup);      
      double CalcObjFunc(void);
      int CalcMultiObjFunc(double * pF, int nObj){ return -1; }
      void WriteSetupToFile(FILE * pFile);
      void WriteWells(FILE * pFile, int type);
      void WriteCost(FILE * pFile, int type);
      void WriteConstraints(FILE * pFile, int type);
      ConstraintABC * GetConstraintPtr(IroncladString pName);
	  void * GetResponseVarGroup(void);

   private :
      double CalcPumpingRate(void);
      double CalcOperationCost(void);
      double CalcCapitalCost(void);
      double CalcMayerCost(void);
      double CalcTreatmentCost(void);
      double LookupPumpCost(double rate, double lift);
      void InitFromFile(void);
      void InitResponseVars(void);
      void InitConstraints(void);
      void InitWells(void);
      void InitLookupTable(void);      
      void InitPlumes(void);
      void ResizePlume(double x, double y, Plume2D * pPlume);

      PatoObjType m_ObjType;
      LmtPenType  m_PenType;

      double m_RateThresh;   //implicitly defines whether or not a well is active

      //cost factors for TOTQ formulation
      double m_ExtRateCF;    //extraction rate cost factor
      double m_InjRateCF;    //injection rate cost factor

      //captial cost factors for Mayer's formulation
      double m_MayerPumpCF;  //pump cost factor
      double m_MayerDrillCF; //drill cost factor

      //captial cost factors for RS Means formulation
      double m_FixWellCF;    //fixed well installation cost factor
      double m_VarWellCF;    //depth-dependent well installation cost factor

      //unit conversion factors for RS Means formulation
      double m_RateUCF;      //lookup table rate unit conversion factor
      double m_LiftUCF;      //lookup table lift unit conversion factor

      //cost factors for operational costs
      double m_TimeFrame;    //remediation time frame (years)
      double m_IntRate;      //interest rate
      double m_LaborRate;    //operational cost labor rate
      double m_ExtEnergyRate; //energy cost rate for extraction
      double m_InjEnergyRate; //energy cost rate for injection
      double m_AnalyticRate;  //analysis cost ($/sample)
      double m_SampleFreq;    //sample events per year
      double m_DisposalRate;  //disposal cost rate
      double m_MaintFactor;   //maintenance factor

      //cost factors for treatment costs (both operational and capital)
      double m_TreatCapCoeff;
      double m_TreatCapExpon;
      double m_TreatOpCoeff;
      double m_TreatOpExpon;

      double m_Costs[11];

      WellStruct * m_pWells;
      int m_MaxNumWells;

      PumpLkupTableStruct * m_pTbl;
      int m_TblSize;

      ParameterGroup * m_pParamGroup; // design variables
      ResponseVarGroup * m_pRespGroup; // response variables

      ConstraintABC * m_pConstraints; //linked-list of constraints

      Plume2D * m_pPlumes; //plume vertices
      int m_NumPlumes; //number of plumes
}; /* end class PATO */

#endif /* PUMP_AND_TREAT_H */
