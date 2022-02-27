/******************************************************************************
File      : KMeansParamInitializer.cpp
Author    : L. Shawn Matott
Copyright : 2022, L. Shawn Matott

Encapsulates the KMeans clustering parameter initialization strategy.

Version History
02-27-2022    lsm   created
******************************************************************************/
#include "ParamInitializerABC.h"
#include "ParameterABC.h"
#include "ParameterGroup.h"
#include "Exception.h"

/******************************************************************************
GetParameterSets()

Apply the algorithm and populate the "pVals" matrix with initial parameter sets.
******************************************************************************/
void KMeansParamInitializer::GetParameterSets(double ** pVals, int start)
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
KMeansParamInitializer::KMeansParamInitializer(ParameterGroup * pParamGroup, FILE * pInFile)
{
   m_pParams = pParamGroup;
   m_NumSets = 0;
   m_NumParams = 0;
   IncCtorCount();
} /* end CTOR */

/******************************************************************************
Destroy()

Free up memory of the initializer.
******************************************************************************/
void KMeansParamInitializer::Destroy(void)
{
   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
Write()

Writes formatted output to pFile.
******************************************************************************/
void KMeansParamInitializer::Write(FILE * pFile, int type)
{
   fprintf(pFile, "****** Parameter Initialization ******\n");
   fprintf(pFile, "Name       : KMeans Clustering\n");
   fprintf(pFile, "Num Sets   : %d\n", m_NumSets);
} /* end Write() */
