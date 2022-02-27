/******************************************************************************
File      : ParamInitializerABC.h
Author    : L. Shawn Matott 
Copyright : 2022, L. Shawn Matott

This is an abstract base class that supports methods for assigning sets of
initial parameter values.

Version History
02-27-2022   lsm   added copyright information and initial comments.
******************************************************************************/
#ifndef PARAM_INITIALIZER_ABC_H
#define PARAM_INITIALIZER_ABC_H

#include "MyHeaderInc.h"

//forward decs
class ParameterGroup;
class ParameterABC;

/******************************************************************************
class ParamInitializerABC (parameter initializer) base class
******************************************************************************/
class ParamInitializerABC
{
    public :
    	virtual ~ParamInitializerABC(void){ DBG_PRINT("ParamInitializerABC::DTOR"); }
		virtual void Destroy(void)=0;
		virtual void GetParameterSets(double ** pVals, int start)=0;
		virtual void Write(FILE * pFile, int type)=0; 
		virtual int GetNumParameterSets(void)=0;
}; /* end class ParamInitializerABC */

/******************************************************************************
class HamedParamInitializer

Initialize sets of parameter values using Hamed's method. This is thought to be 
useful for capturing larger sections of the pareto front in multi-objective problems.
******************************************************************************/
class HamedParamInitializer : public ParamInitializerABC
{
    public:
	~HamedParamInitializer(void){ DBG_PRINT("HamedParamInitializer::DTOR"); Destroy(); }
	void Destroy(void);
	HamedParamInitializer(ParameterGroup * pParamGroup, FILE * pInFile);
	void GetParameterSets(double ** pVals, int start);
	int GetNumParameterSets(void){ return m_NumSets; }
	void Write(FILE * pFile, int type);

private:
	ParameterGroup * m_pParams;
	int m_HamedOffset;
	int m_NumSets;
	int m_NumParams;
}; /* end class HamedParamInitializer */

/******************************************************************************
class KMeansParamInitializer

Initialize parameter values using K-Means clustering. This is thought to be useful
for grouping spatially related parameters.
******************************************************************************/
class KMeansParamInitializer : public ParamInitializerABC
{
    public:
	~KMeansParamInitializer(void){ DBG_PRINT("KMeansParamInitializer::DTOR"); Destroy(); }
	void Destroy(void);
	KMeansParamInitializer(ParameterGroup * pParamGroup, FILE * pInFIle);
	void GetParameterSets(double ** pVals, int start);
	void Write(FILE * pFile, int type);
	int GetNumParameterSets(void){ return m_NumSets; }

private:
	ParameterGroup * m_pParams;
	int m_HamedOffset;
    int m_NumSets;
	int m_NumParams;
}; /* end class KMeansParamInitializer */

#endif /* PARAM_INITIALIZER_ABC_H */






