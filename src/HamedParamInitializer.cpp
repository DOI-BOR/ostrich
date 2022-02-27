/******************************************************************************
File      : HamedParamInitializer.cpp
Author    : L. Shawn Matott
Copyright : 2022, L. Shawn Matott

Encapsulates the Hamed parameter initialization strategy.

Version History
02-27-2022    lsm   created
******************************************************************************/
#include "ParamInitializerABC.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "Exception.h"

/******************************************************************************
GetParameterSets()

Apply the algorithm and populate the "pVals" matrix with initial parameter values.
******************************************************************************/
void HamedParamInitializer::GetParameterSets(double ** pVals, int start)
{
   int i, j, end;

   end = start + m_NumSets;
   for(i = start; i < end; i++)
   {
      for(j = 0; j < m_NumParams; j++)
      {
         pVals[i][j] = m_pParams->GetParamPtr(j)->ConvertInVal(atof("0.0"));
      }
   }
	return; 
} /* end GetParameterSets() */

/******************************************************************************
CTOR
  
Assign member variables.
******************************************************************************/
HamedParamInitializer::HamedParamInitializer(ParameterGroup * pParamGroup, FILE * pInFIle)
{
   m_HamedOffset = 397;
   m_NumSets = 0;
   m_NumParams = 0;
   m_pParams = pParamGroup;
   IncCtorCount();
} /* end CTOR */

/******************************************************************************
Destroy()

Free up memory of the initializer.
******************************************************************************/
void HamedParamInitializer::Destroy(void)
{
   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
Write()

Writes formatted output to pFile.
******************************************************************************/
void HamedParamInitializer::Write(FILE * pFile, int type)
{
   fprintf(pFile, "****** Parameter Initialization ******\n");
   fprintf(pFile, "Name       : Hamed's Method\n");
   fprintf(pFile, "Offset     : %d\n", m_HamedOffset);
   fprintf(pFile, "Num Sets   : %d\n", m_NumSets);
} /* end Write() */

