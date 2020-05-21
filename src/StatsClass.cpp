/******************************************************************************
File     : StatsClass.cpp
Author   : L. Shawn Matott
Copyright: 2003, L. Shawn Matott

The StatsClass is used to compute statistical measures, following a successful
calibration. Many of the statistics require the Jacobian matrix, evaluated at 
the calibrated minimum. If the native calibration algorithm has the Jacobian
available, this can be passed into the StatsClass, otherwise the StatsClass 
will compute the Jacobian internally.

The following statistical measures are avaiable:

   variance/standard deviation   
   covariance/standard error
   correlation coefficients
   Beale's Nonlinearity Measure
   Linssen's Nonlinearity Measure
   Cook's D Observation Influence Measure
   DFBETAS Observation Influence Measure

Version History
06-12-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
11-19-03    lsm   Parameter output uses WRITE_TX_BNR definition.
11-25-03    lsm   Added sensitivity measures, corrected FSTAT in Beale/Linssen
03-24-04    lsm   Added GetLinearity() to support PSO-LevMar hybrid
                  Added CalcOptimalStepSize() option
07-08-04    lsm   Added Jacobian sensitivity checks
01-11-05    lsm   Added metrics
03-21-05    lsm   Added support for parameter-specific relative increments
01-01-07    lsm   StatsClass now uses abstract model base class (ModelABC). 
                  The Jacobian can now be calculated in parallel. Added support
                  for the detection of temporarily insensitive parameters and
                  observations; when detected such insensitive terms are "held",
                  i.e. removed from consideration.
01-01-07    lsm   Additional types for finte difference increments were defined.
                  User selects these types by including the following lines in the
                  MathAndStats section of the input file:
                     DiffIncType <type>, where type is one of:
                        range-relative
                        value-relative
                        absolute
                        optimal
                  Unless optimal is selected, user must also specify the step 
                  size by including the following line in the MathAndStats 
                  section of the input file:
                     DiffIncrement <val1 val2 val3 ... valn>
                        where val1 is the step size for the first parameter, 
                        val2 for the second parameter, etc. If only one value
                        is given, it is applied to all parameters.
01-01-07    lsm   Added Durbin-Watson, Runs Test and MMRI statistics. The use 
                  can select these statistics by adding the following lines to 
                  the MathAndStats section of the input file:
                     RunsTest
                     Durbin-Watson
                     MMRI
07-18-07    lsm   Added support for SuperMUSE
******************************************************************************/
#include "mpi_stub.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "StatsClass.h"
#include "Model.h"
#include "ModelBackup.h"
#include "ObservationGroup.h"
#include "Observation.h"
#include "ResponseVarGroup.h"
#include "RespVarABC.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "ObjectiveFunction.h"
#include "SuperMUSE.h"

#include "Utility.h"
#include "Exception.h"
#include "StatUtility.h"
#include "WriteUtility.h"
#include "SuperMuseUtility.h"

int EvalInitParamsParallel(int np, int id, double ** pList, int size, 
                           ModelABC * pModel, int * num_left);

/******************************************************************************
AdjustObjFunc()

Adjust the given objective function to remove the influence of any insensitive
observations.
******************************************************************************/
double StatsClass::AdjustObjFunc(double val)
{
   double tmp;
   int i;
   //restore 'true' residuals
   CalcResiduals();

   for(i = 0; i < m_NumObs; i++)
   {      
      if(m_bHoldObs[i] == true)
      {
         tmp = m_pResid[i] * m_pResid[i];
         val -= tmp;
      }            
   }/* end for() */

   //restore adjusted residuals
   AdjustResiduals();   

   return val;
}/* end AdjustObjFunc() */

/******************************************************************************
AdjustJacobian()

Adjust the Jacobian matrix, if needed, to eliminate insensitive observations
and/or parameters.
******************************************************************************/
void StatsClass::AdjustJacobian(void)
{
   int i, j, row, col;

   m_bAdjustedJac = true;

   //remove rows of insensitive observations
   if(m_NumHeldObs > 0)
   {
      i = 0;
      for(j = 0; j < m_NumObs; j++)
      {
         if(m_bHoldObs[j] == false)
         {
            i++;
         }
         else
         {
            for(col = 0; col < m_NumParams; col++)
            {
               for(row = i; row < m_NumObs - 1; row++)
               {
                  m_pJacob[row][col] = m_pJacob[row+1][col];
                  m_pJacobUW[row][col] = m_pJacobUW[row+1][col];
               }
               m_pJacob[row][col] = 0.00;
               m_pJacobUW[row][col] = 0.00;
            }/* end for() */
         }/* end else() */
      }/* end for() */
   }/* end if() */

   if(m_NumHeldParams > 0)
   {
      i = 0;
      for(j = 0; j < m_NumParams; j++)
      {
         if(m_bHoldParam[j] == false)
         {
            i++;
         }
         else
         {
            for(row = 0; row < m_NumObs; row++)
            {
               for(col = i; col < m_NumParams - 1; col++)
               {
                  m_pJacob[row][col] = m_pJacob[row][col+1];
                  m_pJacobUW[row][col] = m_pJacobUW[row][col+1];
               }
               m_pJacob[row][col] = 0.00;
               m_pJacobUW[row][col] = 0.00;
            }/* end for() */
         }/* end else() */
      }/* end for() */
   }/* end if() */

   //re-compute the transpose
   for(row = 0; row < m_NumObs; row++)
   {
      for(col = 0; col < m_NumParams; col++)
      {
         m_pJacobT[col][row] = m_pJacob[row][col];
      }
   }/* end for() */   
}/* end AdjustJacobian() */

/******************************************************************************
AdjustResiduals()

Adjust the residuals vector, if needed, to eliminate insensitive observations.
******************************************************************************/
void StatsClass::AdjustResiduals(void)
{ 
   int i, j, k;

   if(m_NumHeldObs > 0)
   {
      i = 0;
      for(j = 0; j < m_NumObs; j++)
      {
         if(m_bHoldObs[j] == false)
         {
            i++;
         }
         else
         {
            for(k = i; k < m_NumObs - 1; k++)
            {
               m_pResid[k] = m_pResid[k+1];
            }
            m_pResid[k] = 0.00;
         }/* end else() */
      }/* end for() */
   }/* end if() */
}/* end AdjustResiduals() */

/******************************************************************************
AdjustVector()

Adjust the given vector, if needed, to INSERT zeroes where insensitive 
observations or parameters are normally located, but have been shifted out
due to previous adjustment calls.
******************************************************************************/
void StatsClass::AdjustVector(double * vec, bool obs)
{
   int i, j, k;

   if((obs == true) && (m_NumHeldObs > 0))
   {
      i = 0;
      for(j = 0; j < m_NumObs; j++)
      {
         if(m_bHoldObs[j] == false)
         {
            i++;
         }
         else
         {
            for(k = m_NumObs - 1; k > i; k--)
            {
               vec[k] = vec[k-1];
            }
            vec[k] = 0.00;
            i++;
         }/* end else() */
      }/* end for() */
   }/* end if() */
   else if((obs == false) && (m_NumHeldParams > 0))
   {
      i = 0;
      for(j = 0; j < m_NumParams; j++)
      {
         if(m_bHoldParam[j] == false)
         {
            i++;
         }
         else
         {
            for(k = m_NumParams - 1; k > i; k--)
            {
               vec[k] = vec[k-1];
            }
            vec[k] = 0.00;
            i++;
         }/* end else() */
      }/* end for() */
   }
}/* end AdjustVector() */

/******************************************************************************
CTOR

Init. everything based on user configuration file.
******************************************************************************/
StatsClass::StatsClass(ModelABC * pModel)
{
   int i;
   ParameterGroup * pParamGroup;
   ObservationGroup * pObsGroup;

   m_Phi = 0.00;
   m_CIpct = 95.00;

   m_DiffType = FD_FORWARD;
   m_MinInc = NEARLY_ZERO;
         
   m_bNoStats = false;
   m_StdDevFlag = false;
   m_StdErrFlag = false;
   m_CorrCoefFlag = false;
   m_NormPlotFlag = false;
   m_BealeFlag = false;
   m_LinssenFlag = false;
   m_CooksFlag = false;
   m_DfbetasFlag = false;     
   m_MatricesFlag = false;   
   m_CIflag = false;
   m_SensFlag = false;
   m_RunsTestFlag = false;
   m_AutorunFunctionFlag = false;
   m_MMRI_Flag = false;
   m_bOkToHoldParams = true;
   m_bOkToHoldObs = true;
   m_BestBoxCoxFlag = false;
   m_bWriteIterationResiduals = false;
      
   //allocate parabolic central diff. matrices
   NEW_PRINT("double *", 3);
   m_ParaMat = new double *[3];

   NEW_PRINT("double *", 3);
   m_ParaInv = new double *[3];
   MEM_CHECK(m_ParaInv);

   for(i = 0; i < 3; i++)
   {   
      NEW_PRINT("double", 3);
      m_ParaMat[i] = new double[3];

      NEW_PRINT("double", 3);
      m_ParaInv[i] = new double[3];
   }
   MEM_CHECK(m_ParaInv[i-1]);

   /* allocate those matrices which are a function of the model setup */
   m_pModel = pModel;
   pParamGroup = m_pModel->GetParamGroupPtr();
   pObsGroup = m_pModel->GetObsGroupPtr();

   m_NumObs = pObsGroup->GetNumObs();
   m_NumParams = pParamGroup->GetNumParams();

   m_pPredictions = NULL;
   m_pJacPred = NULL; //partial derivatives for predictions
   m_Pred = NULL;
   m_PredSD = NULL;
   m_PredCIlwr = NULL;
   m_PredCIupr = NULL;

   m_NumHeldParams = 0;
   m_NumHeldObs = 0;
   m_bAdjustedJac = false;

   NEW_PRINT("bool", m_NumParams);
   m_bHoldParam =  new bool[m_NumParams];
   MEM_CHECK(m_bHoldParam);

   NEW_PRINT("bool", m_NumObs);
   m_bHoldObs = new bool[m_NumObs];
   MEM_CHECK(m_bHoldObs);

   NEW_PRINT("double", m_NumParams + m_NumObs + 1);
   m_pMinJac = new double[m_NumParams + m_NumObs + 1];
   MEM_CHECK(m_pMinJac);

   NEW_PRINT("FiniteDiffType", m_NumParams);
   m_pDType = new FiniteDiffType[m_NumParams];
   MEM_CHECK(m_pDType);

   NEW_PRINT("double", m_NumParams);
   m_pDx = new double[m_NumParams];
   MEM_CHECK(m_pDx)

   NEW_PRINT("double", m_NumParams);
   m_pMid = new double[m_NumParams];
   MEM_CHECK(m_pMid);

   NEW_PRINT("double", m_NumParams);
   m_pHi = new double[m_NumParams];
   MEM_CHECK(m_pHi);

   NEW_PRINT("double", m_NumParams);
   m_pLow = new double[m_NumParams];
   MEM_CHECK(m_pLow);

   NEW_PRINT("double", m_NumParams);
   m_pDiffInc  = new double[m_NumParams];
   MEM_CHECK(m_pDiffInc);

   m_DiffIncType = FD_RANGE_REL;
   for(i = 0; i < m_NumParams; i++){ m_pDiffInc[i] = 0.001;}

   NEW_PRINT("double", m_NumObs);
   m_pCooksD = new double[m_NumObs];

   NEW_PRINT("double", m_NumObs);
   m_pResid = new double[m_NumObs];

   NEW_PRINT("double", m_NumObs);
   m_pOrdResid = new double[m_NumObs];

   NEW_PRINT("double", m_NumObs);
   m_pExpResid = new double[m_NumObs];

   NEW_PRINT("double", m_NumParams);
   m_pCIupr = new double[m_NumParams];   

   NEW_PRINT("double", m_NumParams);
   m_pCIlwr = new double[m_NumParams];   

   NEW_PRINT("double *", m_NumObs);
   m_pJacob = new double*[m_NumObs];

   NEW_PRINT("double *", m_NumObs);
   m_pJacobUW = new double*[m_NumObs];

   NEW_PRINT("double *", m_NumObs);
   m_pDFBETAS = new double*[m_NumObs];

   NEW_PRINT("double *", m_NumParams);
   m_pJacobT = new double*[m_NumParams];   

   NEW_PRINT("double *", m_NumParams);
   m_pNormal = new double*[m_NumParams];   

   NEW_PRINT("double *", m_NumParams);
   m_pInvNormal = new double*[m_NumParams];   

   NEW_PRINT("double *", m_NumParams);
   m_pPbyO1 = new double*[m_NumParams];   

   NEW_PRINT("double *", m_NumObs);
   m_pHat = new double*[m_NumObs];

   NEW_PRINT("double *", m_NumParams);
   m_pChange = new double*[m_NumParams];
   
   NEW_PRINT("double *", m_NumObs);
   m_pScaledSens = new double*[m_NumObs];

   NEW_PRINT("double *", m_NumObs);
   m_pPctScaledSens = new double*[m_NumObs];

   NEW_PRINT("double", m_NumParams);
   m_pCompScaledSens = new double[m_NumParams];
   
   NEW_PRINT("double *", m_NumParams);
   m_pCovar = new double*[m_NumParams];
   MEM_CHECK(m_pCompScaledSens);

   for(i = 0; i < m_NumObs; i++)
   {
      NEW_PRINT("double", m_NumParams);
      m_pJacob[i] = new double[m_NumParams];
   }   
   for(i = 0; i < m_NumObs; i++)
   {
      NEW_PRINT("double", m_NumParams);
      m_pJacobUW[i] = new double[m_NumParams];
   }   
   for(i = 0; i < m_NumObs; i++)
   {
      NEW_PRINT("double", m_NumParams);
      m_pDFBETAS[i] = new double[m_NumParams];
   }   
   for(i = 0; i < m_NumParams; i++)
   {
      NEW_PRINT("double", m_NumObs);
      m_pJacobT[i] = new double[m_NumObs];
   }   
   for(i = 0; i < m_NumParams; i++)
   {
      NEW_PRINT("double", m_NumParams);
      m_pNormal[i] = new double[m_NumParams];
   }
   for(i = 0; i < m_NumParams; i++)
   {
      NEW_PRINT("double", m_NumParams);
      m_pInvNormal[i] = new double[m_NumParams];
   }
   for(i = 0; i < m_NumParams; i++)
   {
      NEW_PRINT("double", m_NumObs);
      m_pPbyO1[i] = new double[m_NumObs];
   }
   for(i = 0; i < m_NumObs; i++)
   {
      NEW_PRINT("double", m_NumObs);
      m_pHat[i] = new double[m_NumObs];
   }
   for(i = 0; i < m_NumParams; i++)
   {
      NEW_PRINT("double", m_NumObs);
      m_pChange[i] = new double[m_NumObs];
   }
   for(i = 0; i < m_NumObs; i++)
   {
      NEW_PRINT("double", m_NumParams);
      m_pScaledSens[i] = new double[m_NumParams];
   }
   for(i = 0; i < m_NumParams; i++)
   {
      NEW_PRINT("double", m_NumParams);
      m_pCovar[i] = new double[m_NumParams];
   }
   for(i = 0; i < m_NumObs; i++)
   {
      NEW_PRINT("double", m_NumParams);
      m_pPctScaledSens[i] = new double[m_NumParams];
   }
   MEM_CHECK(m_pPctScaledSens[i-1]);
    
   //setup model backups
   NEW_PRINT("ModelBackup", 1);
   m_pMidBkup = new ModelBackup(m_pModel);

   NEW_PRINT("ModelBackup", 1);
   m_pLowBkup = new ModelBackup(m_pModel);

   NEW_PRINT("ModelBackup", 1);
   m_pHiBkup = new ModelBackup(m_pModel);
   MEM_CHECK(m_pHiBkup);

   m_pBuf = NULL;
   m_pMyBuf = NULL;

   //configuration file can override certain defaults
   InitFromFile(GetInFileName());

   IncCtorCount();
} /* end CTOR */

/******************************************************************************
InitFromFile()

Read configuration information from the given filename.
******************************************************************************/
void StatsClass::InitFromFile(IroncladString pStatsFileName)
{
   FILE * pFile;
   char * line;
   char tmp[DEF_STR_SZ];
   char tmp2[DEF_STR_SZ];
   char * pTmp;
   int i,j;

   m_DiffCount = 0;
   m_StepCount = 0;
   m_StatsCount = 0;

   //read in statistics configuration
   pFile = fopen(pStatsFileName, "r");
   if(pFile == NULL) 
   {
      //couldn't open file, use defaults and log the error.
      LogError(ERR_FILE_IO, "Couldn't open stats. config. file. Using Defaults");      
      return;
   }/* end if() */   

   //make sure correct tokens are present
   if(CheckToken(pFile, "BeginMathAndStats", pStatsFileName) == true)
   {
      FindToken(pFile, "EndMathAndStats", pStatsFileName);
      rewind(pFile);

      /*
      User has specified config file, clear all stat flags as they will be
      overridden.
      */
      m_bNoStats = false;
      m_StdDevFlag = false;
      m_StdErrFlag = false;
      m_CorrCoefFlag = false;
      m_NormPlotFlag = false;
      m_BealeFlag = false;
      m_LinssenFlag = false;
      m_CooksFlag = false;
      m_DfbetasFlag = false;
      m_MatricesFlag = false;    
      m_CIflag = false;
      m_SensFlag = false;
      m_RunsTestFlag = false;
      m_AutorunFunctionFlag = false;
      m_MMRI_Flag = false;
      m_bWriteIterationResiduals = false;

      FindToken(pFile, "BeginMathAndStats", pStatsFileName);
      line = GetNxtDataLine(pFile, pStatsFileName);
      while(strstr(line, "EndMathAndStats") == NULL)
      {
         if(strstr(line, "DiffType") != NULL)
         {
            sscanf(line, "%s %s", tmp, tmp2);
            strcpy(line, tmp2);
            MyStrLwr(line);
            if(strstr(line, "forward") != NULL) {m_DiffType = FD_FORWARD;}
            else if(strstr(line, "outside") != NULL) {m_DiffType = FD_OUT_CEN;}
            else if(strstr(line, "parabolic") != NULL) {m_DiffType = FD_PAR_CEN;}
            else if(strstr(line, "best-fit") != NULL) {m_DiffType = FD_FIT_CEN;}
         }/*end if() */
         else if(strstr(line, "DiffIncType") != NULL)
         {
            sscanf(line, "%s %s", tmp, tmp2);
            strcpy(line, tmp2);
            MyStrLwr(line);
            if     (strstr(line, "range-relative") != NULL) {m_DiffIncType = FD_RANGE_REL;}
            else if(strstr(line, "value-relative") != NULL) {m_DiffIncType = FD_VALUE_REL;}
            else if(strstr(line, "absolute")       != NULL) {m_DiffIncType = FD_ABSOLUTE;}
            else if(strstr(line, "optimal")        != NULL) {m_DiffIncType = FD_OPTIMAL;}
         }/*end if() */
         else if(strstr(line, "DiffRelIncrement") != NULL)
         {
            MyStrLwr(line);
            pTmp = line + (int)strlen("DiffRelIncrement");

            for(i = 0; i < m_NumParams; i++)
            {
                j = ExtractString(pTmp, tmp);
                m_pDiffInc[i] = atof(tmp);
                if(j == -1){ break;}
                pTmp += j;               
            }

            for(i = i; i < m_NumParams; i++)
            { 
               m_pDiffInc[i] = m_pDiffInc[0];
            }

            /* this keyword (DiffRelIncrement) is range-relative */
            m_DiffIncType = FD_RANGE_REL;
         }/*end else if() */         
         else if(strstr(line, "DiffIncrement") != NULL)
         {
            MyStrLwr(line);
            pTmp = line + (int)strlen("DiffIncrement");

            for(i = 0; i < m_NumParams; i++)
            {
                j = ExtractString(pTmp, tmp);                
                m_pDiffInc[i] = atof(tmp);
                if(j == -1){ break;}
                pTmp += j;               
            }

            for(i = i; i < m_NumParams; i++)
            { 
               m_pDiffInc[i] = m_pDiffInc[0];
            }
         }/*end else if() */
         else if(strstr(line, "DiffMinIncrement") != NULL)
         {
            sscanf(line, "DiffMinIncrement %lf", &m_MinInc);
         }/*end else if() */
         else if(strstr(line, "CI_Pct") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_CIpct); 
            if((m_CIpct < 0.00) || (m_CIpct > 100.00)){ m_CIpct = 95.00;}            
         }
         else if(strstr(line, "Default") != NULL)
         {
            m_bNoStats = false;
            m_StdDevFlag = true;
            m_StdErrFlag = true;
            m_CorrCoefFlag = true;
            m_NormPlotFlag = false;
            m_BealeFlag = false;
            m_LinssenFlag = false;
            m_CooksFlag = false;
            m_DfbetasFlag = false;
            m_MatricesFlag = false;
            m_CIflag = false;
            m_SensFlag = false;
            m_RunsTestFlag = false;
            m_AutorunFunctionFlag = false;
            m_MMRI_Flag = false;
            m_bOkToHoldParams = true;
            m_bOkToHoldObs = true;
            m_BestBoxCoxFlag = false;
         }
         else if(strstr(line, "AllStats") != NULL)
         {
            m_bNoStats = false;
            m_StdDevFlag = true;
            m_StdErrFlag = true;
            m_CorrCoefFlag = true;
            m_NormPlotFlag = true;
            m_BealeFlag = true;
            m_LinssenFlag = true;
            m_CooksFlag = true;
            m_DfbetasFlag = true;
            m_MatricesFlag = true;
            m_CIflag = true;
            m_SensFlag = true;
            m_RunsTestFlag = true;
            m_AutorunFunctionFlag = true;
            m_MMRI_Flag = true;
            m_bOkToHoldParams = true;
            m_bOkToHoldObs = true;
            m_BestBoxCoxFlag = true;
         }
         else if(strstr(line, "NoStats") != NULL)
         {
            m_bNoStats = true;
            m_StdDevFlag = false;
            m_StdErrFlag = false;
            m_CorrCoefFlag = false;
            m_NormPlotFlag = false;
            m_BealeFlag = false;
            m_LinssenFlag = false;
            m_CooksFlag = false;
            m_DfbetasFlag = false;
            m_MatricesFlag = false;
            m_CIflag = false;
            m_SensFlag = false;
            m_RunsTestFlag = false;
            m_AutorunFunctionFlag = false;
            m_MMRI_Flag = false;
            m_bOkToHoldParams = false;
            m_bOkToHoldObs = false;
            m_BestBoxCoxFlag = false; 
         } 
         else if(strstr(line, "BestBoxCox") != NULL){m_bNoStats = false; m_BestBoxCoxFlag = true;}
         else if(strstr(line, "StdDev") != NULL){m_bNoStats = false; m_StdDevFlag = true;}
         else if(strstr(line, "StdErr") != NULL){m_bNoStats = false; m_StdErrFlag = true;}
         else if(strstr(line, "CorrCoeff") != NULL){m_bNoStats = false; m_CorrCoefFlag = true;}
         else if(strstr(line, "NormPlot") != NULL){m_bNoStats = false; m_NormPlotFlag = true;}
         else if(strstr(line, "Beale") != NULL){m_bNoStats = false; m_BealeFlag = true;}       
         else if(strstr(line, "Linssen") != NULL){m_bNoStats = false; m_LinssenFlag = true;}
         else if(strstr(line, "CooksD") != NULL){m_bNoStats = false; m_CooksFlag = true;}
         else if(strstr(line, "DFBETAS") != NULL){m_bNoStats = false; m_DfbetasFlag = true;}
         else if(strstr(line, "Matrices") != NULL){m_bNoStats = false; m_MatricesFlag = true;}
         else if(strstr(line, "Confidence") != NULL){m_bNoStats = false; m_CIflag = true;}
         else if(strstr(line, "Sensitivity") != NULL){m_bNoStats = false; m_SensFlag = true;}
         else if(strstr(line, "RunsTest") != NULL){m_bNoStats = false; m_RunsTestFlag = true;}
         else if(strstr(line, "AutorunFunction") != NULL){m_bNoStats = false; m_AutorunFunctionFlag = true;}
         else if(strstr(line, "MMRI") != NULL){m_bNoStats = false; m_MMRI_Flag = true;}
         else if(strstr(line, "ExcludeInsensitiveParameters")   != NULL){ m_bOkToHoldParams = true;}
         else if(strstr(line, "IncludeInsensitiveParameters")   != NULL){ m_bOkToHoldParams = false;}
         else if(strstr(line, "ExcludeInsensitiveObservations") != NULL){ m_bOkToHoldObs = true;}
         else if(strstr(line, "IncludeInsensitiveObservations") != NULL){ m_bOkToHoldObs = false;}
         else if(strstr(line, "WriteResidualsEachIteration") != NULL){ m_bWriteIterationResiduals = true; }         
         else
         {
            sprintf(tmp, "Unknown token: %s", line);
            LogError(ERR_FILE_IO, tmp);
         }/* end else() */
         line = GetNxtDataLine(pFile, pStatsFileName);
      } /* end while() */
   }/* end if() */   

   //make sure correct tokens are present
   if(CheckToken(pFile, "BeginMoreStats", pStatsFileName) == true)
   {
      FindToken(pFile, "EndMoreStats", pStatsFileName);
      rewind(pFile);

      FindToken(pFile, "BeginMoreStats", pStatsFileName);
      line = GetNxtDataLine(pFile, pStatsFileName);
      while(strstr(line, "EndMoreStats") == NULL)
      {
         if(strstr(line, "DiffType") != NULL)
         {
            sscanf(line, "%s %s", tmp, tmp2);
            strcpy(line, tmp2);
            MyStrLwr(line);
            if(strstr(line, "forward") != NULL) {m_DiffType = FD_FORWARD;}
            else if(strstr(line, "outside") != NULL) {m_DiffType = FD_OUT_CEN;}
            else if(strstr(line, "parabolic") != NULL) {m_DiffType = FD_PAR_CEN;}
            else if(strstr(line, "best-fit") != NULL) {m_DiffType = FD_FIT_CEN;}
         }/*end if() */
         else if(strstr(line, "DiffIncType") != NULL)
         {
            sscanf(line, "%s %s", tmp, tmp2);
            strcpy(line, tmp2);
            MyStrLwr(line);
            if     (strstr(line, "range-relative") != NULL) {m_DiffIncType = FD_RANGE_REL;}
            else if(strstr(line, "value-relative") != NULL) {m_DiffIncType = FD_VALUE_REL;}
            else if(strstr(line, "absolute")       != NULL) {m_DiffIncType = FD_ABSOLUTE;}
            else if(strstr(line, "optimal")        != NULL) {m_DiffIncType = FD_OPTIMAL;}
         }/*end if() */
         else if(strstr(line, "DiffRelIncrement") != NULL)
         {
            MyStrLwr(line);
            pTmp = line + (int)strlen("DiffRelIncrement");

            for(i = 0; i < m_NumParams; i++)
            {
                j = ExtractString(pTmp, tmp);
                m_pDiffInc[i] = atof(tmp);
                if(j == -1){ break;}
                pTmp += j;               
            }

            for(i = i; i < m_NumParams; i++)
            { 
               m_pDiffInc[i] = m_pDiffInc[0];
            }

            /* this keyword (DiffRelIncrement) is range-relative */
            m_DiffIncType = FD_RANGE_REL;
         }/*end else if() */         
         else if(strstr(line, "DiffIncrement") != NULL)
         {
            MyStrLwr(line);
            pTmp = line + (int)strlen("DiffIncrement");

            for(i = 0; i < m_NumParams; i++)
            {
                j = ExtractString(pTmp, tmp);                
                m_pDiffInc[i] = atof(tmp);
                if(j == -1){ break;}
                pTmp += j;               
            }

            for(i = i; i < m_NumParams; i++)
            { 
               m_pDiffInc[i] = m_pDiffInc[0];
            }
         }/*end else if() */
         else if(strstr(line, "DiffMinIncrement") != NULL)
         {
            sscanf(line, "DiffMinIncrement %lf", &m_MinInc);
         }/*end else if() */
         else if(strstr(line, "CI_Pct") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_CIpct); 
            if((m_CIpct < 0.00) || (m_CIpct > 100.00)){ m_CIpct = 95.00;}            
         }
         else if(strstr(line, "Default") != NULL)
         {
            m_bNoStats = false;
            m_StdDevFlag = true;
            m_StdErrFlag = true;
            m_CorrCoefFlag = true;
            m_NormPlotFlag = false;
            m_BealeFlag = false;
            m_LinssenFlag = false;
            m_CooksFlag = false;
            m_DfbetasFlag = false;
            m_MatricesFlag = false;
            m_CIflag = false;
            m_SensFlag = false;
            m_RunsTestFlag = false;
            m_AutorunFunctionFlag = false;
            m_MMRI_Flag = false;
            m_bOkToHoldParams = true;
            m_bOkToHoldObs = true;
            m_BestBoxCoxFlag = false;
         }
         else if(strstr(line, "AllStats") != NULL)
         {
            m_bNoStats = false;
            m_StdDevFlag = true;
            m_StdErrFlag = true;
            m_CorrCoefFlag = true;
            m_NormPlotFlag = true;
            m_BealeFlag = true;
            m_LinssenFlag = true;
            m_CooksFlag = true;
            m_DfbetasFlag = true;
            m_MatricesFlag = true;
            m_CIflag = true;
            m_SensFlag = true;
            m_RunsTestFlag = true;
            m_AutorunFunctionFlag = true;
            m_MMRI_Flag = true;
            m_bOkToHoldParams = true;
            m_bOkToHoldObs = true;
            m_BestBoxCoxFlag = true;
         }
         else if(strstr(line, "NoStats") != NULL)
         {
            m_bNoStats = true;
            m_StdDevFlag = false;
            m_StdErrFlag = false;
            m_CorrCoefFlag = false;
            m_NormPlotFlag = false;
            m_BealeFlag = false;
            m_LinssenFlag = false;
            m_CooksFlag = false;
            m_DfbetasFlag = false;
            m_MatricesFlag = false;
            m_CIflag = false;
            m_SensFlag = false;
            m_RunsTestFlag = false;
            m_AutorunFunctionFlag = false;
            m_MMRI_Flag = false;
            m_bOkToHoldParams = false;
            m_bOkToHoldObs = false;
            m_BestBoxCoxFlag = false; 
         } 
         else if(strstr(line, "BestBoxCox") != NULL){m_BestBoxCoxFlag = true;}
         else if(strstr(line, "StdDev") != NULL){m_StdDevFlag = true;}
         else if(strstr(line, "StdErr") != NULL){m_StdErrFlag = true;}
         else if(strstr(line, "CorrCoeff") != NULL){m_CorrCoefFlag = true;}
         else if(strstr(line, "NormPlot") != NULL){m_NormPlotFlag = true;}
         else if(strstr(line, "Beale") != NULL){m_BealeFlag = true;}       
         else if(strstr(line, "Linssen") != NULL){m_LinssenFlag = true;}
         else if(strstr(line, "CooksD") != NULL){m_CooksFlag = true;}
         else if(strstr(line, "DFBETAS") != NULL){m_DfbetasFlag = true;}
         else if(strstr(line, "Matrices") != NULL){m_MatricesFlag = true;}
         else if(strstr(line, "Confidence") != NULL){m_CIflag = true;}
         else if(strstr(line, "Sensitivity") != NULL){m_SensFlag = true;}
         else if(strstr(line, "RunsTest") != NULL){m_RunsTestFlag = true;}
         else if(strstr(line, "AutorunFunction") != NULL){m_AutorunFunctionFlag = true;}
         else if(strstr(line, "MMRI") != NULL){m_MMRI_Flag = true;}
         else if(strstr(line, "ExcludeInsensitiveParameters")   != NULL){ m_bOkToHoldParams = true;}
         else if(strstr(line, "IncludeInsensitiveParameters")   != NULL){ m_bOkToHoldParams = false;}
         else if(strstr(line, "ExcludeInsensitiveObservations") != NULL){ m_bOkToHoldObs = true;}
         else if(strstr(line, "IncludeInsensitiveObservations") != NULL){ m_bOkToHoldObs = false;}
         else
         {
            sprintf(tmp, "Unknown token: %s", line);
            LogError(ERR_FILE_IO, tmp);
         }/* end else() */
         line = GetNxtDataLine(pFile, pStatsFileName);
      } /* end while() */
   }/* end if() */   

   fclose(pFile);
} /* end InitFromFile() */

/******************************************************************************
Destroy()

Free up arrays.
******************************************************************************/
void StatsClass::Destroy(void)
{
   int i;
   int nrv = 0;   

   if(m_pPredictions != NULL)
   {
      nrv = m_pPredictions->GetNumRespVars();
      delete m_pPredictions;
   }

   delete [] m_pMinJac;
   delete [] m_pBuf;
   delete [] m_pMyBuf;
   delete [] m_pCooksD;
   delete [] m_pResid;
   delete [] m_pOrdResid;
   delete [] m_pExpResid;
   delete [] m_bHoldParam;
   delete [] m_bHoldObs;

   for(i = 0; i < m_NumObs; i++){delete [] m_pDFBETAS[i];}
   delete [] m_pDFBETAS;

   for(i = 0; i < m_NumObs; i++){delete [] m_pJacob[i];}
   delete [] m_pJacob;

   for(i = 0; i < m_NumObs; i++){delete [] m_pJacobUW[i];}
   delete [] m_pJacobUW;

   for(i = 0; i < nrv; i++){delete [] m_pJacPred[i];}
   delete [] m_pJacPred;

   for(i = 0; i < m_NumObs; i++){delete [] m_pScaledSens[i];}
   delete [] m_pScaledSens;

   for(i = 0; i < m_NumObs; i++){delete [] m_pPctScaledSens[i];}
   delete [] m_pPctScaledSens;

   delete [] m_pCompScaledSens;

   for(i = 0; i < m_NumParams; i++){delete [] m_pJacobT[i];}
   delete [] m_pJacobT;

   for(i = 0; i < m_NumParams; i++){delete [] m_pPbyO1[i];}
   delete [] m_pPbyO1;

   for(i = 0; i < m_NumObs; i++){delete [] m_pHat[i];}
   delete [] m_pHat;

   for(i = 0; i < m_NumParams; i++){delete [] m_pChange[i];}
   delete [] m_pChange;

   for(i = 0; i < m_NumParams; i++){delete [] m_pNormal[i];}
   delete [] m_pNormal;

   for(i = 0; i < m_NumParams; i++){delete [] m_pInvNormal[i];}
   delete [] m_pInvNormal;

   for(i = 0; i < m_NumParams; i++){delete [] m_pCovar[i];}
   delete [] m_pCovar;

   //free up parabolic central diff. matrices
   for(i = 0; i < 3; i++){   
      delete m_ParaMat[i];
      delete m_ParaInv[i];}
   delete [] m_ParaMat;
   delete [] m_ParaInv;
   delete m_pMidBkup;
   delete m_pLowBkup;
   delete m_pHiBkup;
   delete [] m_pCIupr;
   delete [] m_pCIlwr;

   delete [] m_pDiffInc;
   delete [] m_pDType;
   delete [] m_pDx;
   delete [] m_pMid;
   delete [] m_pHi;
   delete [] m_pLow;

   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
CalcResiduals()

Compute the differences between current model-computed obsertvation values and
the experimental observation values read in from the model output file.
******************************************************************************/
Unchangeable1DArray StatsClass::CalcResiduals(void)
{   
   int i;
   ObservationGroup * pObsGroup;
   Observation * pObs;

   pObsGroup = m_pModel->GetObsGroupPtr();
   //calculate vector of residuals
   for(i = 0; i < m_NumObs; i++)
   {
      pObs = pObsGroup->GetObsPtr(i);
      m_pResid[i] = pObs->CalcResidual(true, true);
   }/* end for() */

   return m_pResid;
}/* end CalcResiduals() */

/******************************************************************************
CalcJacobian()

Calculate the Jacobian matrix (matrix of dObs/dParams). Four finite difference
methods can be used: forward, outside central, parabolic central and best-fit
central. The default is forward, but this can be changed in the input file 
along  with the relative increments.
******************************************************************************/
Unchangeable2DArray StatsClass::CalcJacobian(double * pBestSavedF)
{   
   Unchangeable2DArray rVal = CalcJacobian(true, true, pBestSavedF);
   //compute unweighted Jacobian
   for(int i = 0; i < m_NumObs; i++)
   {
      double wt = GetObsWeight(m_pModel->GetObsGroupPtr()->GetObsPtr(i));
      for(int j = 0; j < m_NumParams; j++)
      {
         m_pJacobUW[i][j] = UnWeightJacobian(m_pJacob[i][j], wt);
      }/* end for() */
   }/* end for() */
   return rVal;
}/* end CalcJacobian(); */

/******************************************************************************
CalcJacobian()

Calculate the Jacobian matrix (matrix of dObs/dParams). Four finite difference
methods can be used: forward, outside central, parabolic central and best-fit
central. The default is forward, but this can be changed in the input file 
along  with the relative increments.

Takes 2 boolean arguments that indicate whether or not observations and 
parameters should be 'held' if they are insensitive.
******************************************************************************/
Unchangeable2DArray StatsClass::CalcJacobian(bool bOkToHoldParams, bool bOkToHoldObs, double * pBestSavedF)
{   
   int n, i, id, j;
   char msg[DEF_STR_SZ];
   FILE * pFile;
   double sumCol, sumRow, sumTot;
   ObservationGroup * pObsGroup = m_pModel->GetObsGroupPtr();
   Observation * pObs;
   ParameterGroup * pParamGroup = m_pModel->GetParamGroupPtr();
   ParameterABC * pParam;

   MPI_Comm_size(MPI_COMM_WORLD, &n);
   MPI_Comm_rank(MPI_COMM_WORLD, &id);

   m_bAdjustedJac = false;
   m_NumHeldParams = 0;
   m_NumHeldObs = 0;
   for(i = 0; i < m_NumObs; i++){ m_bHoldObs[i] = false;}
   for(i = 0; i < m_NumParams; i++){ m_bHoldParam[i] = false;}
   
   if(n == 1) //serial execution
   {
      if(IsSuperMUSE() == false)
      {
         EvalJacSerial(pBestSavedF);
      }
      else
      {
         EvalJacSuperMUSE();
      }
   }/* end if() */
   else /* parallel execution */
   {
      BcastJacobian();
      EvalJacParallel();
   }/* end else() */

   //Perform sensitivity checks
   if(id == 0)
   {
      sumTot = 0.00;
      for(i = 0; i < m_NumObs; i++)
      {
         for(j = 0; j < m_NumParams; j++)
         {
            sumTot += fabs(m_pJacob[i][j]);
         }
      }
      if(sumTot <= NEARLY_ZERO)
      {
         sprintf(msg, "Jacobian matrix is completely insensitive");
         LogError(ERR_JACOBIAN, msg);
#ifndef ISOFIT_BUILD
         if(GetProgramType() != JACOBIAN_PROGRAM)
           goto exit_program;
#endif
      }

      //check for observation insensitivity
      for(i = 0; i < m_NumObs; i++)
      {
         sumRow = 0.00;
         for(j = 0; j < m_NumParams; j++)
         {
            sumRow += fabs(m_pJacob[i][j]);
            if(sumRow > NEARLY_ZERO) break;
         }
         if(sumRow <= NEARLY_ZERO)
         {
            if(bOkToHoldObs == true)
            {
               m_bHoldObs[i] = true;
               m_NumHeldObs++;
            }
            else
            {
               pObs = pObsGroup->GetObsPtr(i);
               sprintf(msg, "%s appears to be insensitive", pObs->GetName());
               LogError(ERR_INS_OBS, msg);
            }
         }/* end if() */
      }/* end for() */

      if(m_NumHeldObs > 0)
      {
         sprintf(msg, "Jacobian has %d insensitive observations", m_NumHeldObs);
         LogError(ERR_JACOBIAN, msg);
      }

      //check for parameter insensitivity
      for(j = 0; j < m_NumParams; j++)
      {
         sumCol = 0.00;
         for(i = 0; i < m_NumObs; i++)
         {
            sumCol += fabs(m_pJacob[i][j]);
            if(sumCol > NEARLY_ZERO) break;
         }
         if(sumCol <= NEARLY_ZERO)
         {
            if(bOkToHoldParams == true)
            {
               m_bHoldParam[j] = true;
               m_NumHeldParams++;
            }
            else
            {
               pParam = pParamGroup->GetParamPtr(j);
               sprintf(msg, "%s appears to be insensitive", pParam->GetName());
               LogError(ERR_INS_PARM, msg);
            }
         }/* end if() */
      }/* end for() */

      if(m_NumHeldParams > 0)
      {
         sprintf(msg, "Jacobian has %d insensitive parameters", m_NumHeldParams);
         LogError(ERR_JACOBIAN, msg);
      }
   }/* end if(master processor) */

   return(m_pJacob);

#ifndef ISOFIT_BUILD
exit_program:
#endif
   m_StdDevFlag = m_StdErrFlag = m_CorrCoefFlag = false;
   m_NormPlotFlag = m_BealeFlag = false;
   m_LinssenFlag = m_CooksFlag = m_DfbetasFlag = false;
   m_MatricesFlag = m_CIflag = m_SensFlag = false;

   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   sprintf(msg, "OstOutput%d.txt", id);
   pFile = fopen(msg, "a");   
   WriteStats(pFile);
   fclose(pFile);

   ExitProgram(1);
   return NULL;
} /* end CalcJacobian() */

/******************************************************************************
BcastJacobian()

When in parallel, BcastJacobian() routine is called upon to broadcast the 
current parameter, objective function and observation values from the master 
processor to all of the slave processors.
******************************************************************************/
void StatsClass::BcastJacobian(void)
{   
   int buf_size;
   int num_procs, id;

   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   MPI_Comm_rank(MPI_COMM_WORLD, &id);

   //size the flattened variable matrix   
   buf_size = 1 + m_NumParams + m_NumObs;
   if(m_pBuf == NULL)
   {
      NEW_PRINT("double", buf_size);
      m_pBuf = new double[buf_size];   
      MEM_CHECK(m_pBuf);
   }

   //first element is objective function value
   m_pBuf[0] = m_pModel->GetObjFuncVal();
   //next elements are the parameter settings
   m_pModel->GetParamGroupPtr()->ReadParams(&(m_pBuf[1]));
   //final elements are the simulated observation values
   m_pModel->GetObsGroupPtr()->ReadObservations(&(m_pBuf[1+m_NumParams]));

   //broadcast the flattened matrix
   MPI_Bcast(m_pBuf, buf_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

   //use the flattened matrix to adjust model and parameter and observation groups
   m_pModel->SetObjFuncVal(m_pBuf[0]);
   m_pModel->GetParamGroupPtr()->WriteParams(&(m_pBuf[1]));
   m_pModel->GetObsGroupPtr()->WriteObservations(&(m_pBuf[1+m_NumParams]));
}/* end BcastJacobian() */

/******************************************************************************
BcastMinJac()

When in parallel, BcastMinJac() routine is called upon to collect the minimum
configuration discovered by the most recent Jacobian evaluation.
******************************************************************************/
void StatsClass::BcastMinJac(void)
{   
   int i, proc;
   int buf_size;
   int num_procs, id;

   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   MPI_Comm_rank(MPI_COMM_WORLD, &id);

   //size the flattened variable matrix   
   buf_size = 1 + m_NumParams + m_NumObs;
   if(m_pBuf == NULL)
   {
      NEW_PRINT("double", buf_size);
      m_pBuf = new double[buf_size];   
      MEM_CHECK(m_pBuf);
   }

   //initializing the broadcast buffer
   for(i = 0; i < buf_size; i++){ m_pBuf[i] = m_pMinJac[i];}

   for(proc = 0; proc < num_procs; proc++)
   {
      //broadcast the min Jac data to other processors
      MPI_Bcast(m_pBuf, buf_size, MPI_DOUBLE, proc, MPI_COMM_WORLD);

      //revise min Jac, if the broadcasted processor data is better
      if(m_pBuf[0] < m_pMinJac[0])
      {
         for(i = 0; i < buf_size; i++){ m_pMinJac[i] = m_pBuf[i];}
      }
      //restore broadcast buffer
      for(i = 0; i < buf_size; i++){ m_pBuf[i] = m_pMinJac[i];}
   }/* end for() */
}/* end BcastMinJac() */

/******************************************************************************
EvalJacParallel()

Compute the Jocobian matrix in parallel. Each processor evaluates a pre-
determined number of parameter sets, based on their processor id.
******************************************************************************/
void StatsClass::EvalJacParallel(void)
{    
   int i ,j, num_procs, id, bufsize, idx;
   bool flipSign;
   double cur, next, F;
   double upr;    //upper bound of parameter
   double lwr;    //lower bound of parameter
   double dx;     //parameter perterbation
   double diff; //derivative approximation, to be stored in Jacobian
   double dObs;  //forward and outside central diff. numerator
   double lowParam; //central diff. LHS
   double hiParam; //central diff. RHS and forward diff. RHS
   double midParam; //central diff middle parameter value, forward diff. LHS
   double lowObs; //central diff. LHS
   double hiObs; //central diff. RHS and forward diff.
   double midObs; //current observation value
   double paraCof[3];    //coeff. of parabolic solution
   double paraObs[3];    //obs. vals. of parabolic solution
   double Sxy, Sx, Sy, Sxx; //central diff. best-fit parameters   
   double total_diff;
   ParameterGroup * pParamGroup;
   ParameterABC * pParam;
   double * pPoint;
   FiniteDiffType dType;
   FiniteDiffIncType dIncType;
   
   //setup processor id and number of processors
   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

   /* ------------------------------------------
   Intialize the min. Jacobian
   ------------------------------------------ */
   //first element is objective function value
   m_pMinJac[0] = m_pModel->GetObjFuncVal();
   //next elements are the parameter settings
   m_pModel->GetParamGroupPtr()->ReadParams(&(m_pMinJac[1]));
   //final elements are the simulated observation values
   m_pModel->GetObsGroupPtr()->ReadObservations(&(m_pMinJac[1+m_NumParams]));

   //allocate space for intermediate buffers
   bufsize = m_NumParams*m_NumObs;
   if(m_pMyBuf == NULL)
   {   
      NEW_PRINT("double", bufsize);
      m_pMyBuf = new double[bufsize];
      MEM_CHECK(m_pMyBuf);
   }
   for(i = 0; i < bufsize; i++) m_pMyBuf[i] = 0.00;

   pParamGroup = m_pModel->GetParamGroupPtr();
   m_pMidBkup->Store();
   flipSign = false;

   /*------------------------------------------------------
   Perterb each parameter, one at a time, and compute 
   resulting change in observations and store in a flattened 
   Jacobian matrix.
   -------------------------------------------------------*/
   for(j = id; j < m_NumParams; j+=num_procs)
   {  
      dType = m_DiffType;
      dIncType = m_DiffIncType;

//if FD calculation is ~0.00, retry on same processor using an alternative increment type
retry:
      pParam = pParamGroup->GetParamPtr(j);

      cur = pParam->GetEstVal();
      upr = pParam->GetUprBnd();
      lwr = pParam->GetLwrBnd();
      if(dIncType == FD_OPTIMAL)
      {
         NEW_PRINT("double", m_NumParams);
         pPoint = new double[m_NumParams];
         MEM_CHECK(pPoint);

         pParamGroup->ReadParams(pPoint);
         dx = CalcOptimalStepSize(j, pPoint);
         delete [] pPoint;      
      }
      else if(dIncType == FD_RANGE_REL)
      {
         dx = fabs(m_pDiffInc[j]*(upr-lwr));
      }
      else if(dIncType == FD_VALUE_REL)
      {
         dx = MyMax(fabs(m_pDiffInc[j]*cur), m_MinInc);
      }
      else if(dIncType == FD_ABSOLUTE)
      {
         dx = fabs(m_pDiffInc[j]);
      }
      else //default to range-relative increment
      {
         dx = fabs(m_pDiffInc[j]*(upr-lwr));
      }
      //trick from NR in C
      next = cur + dx;
      dx = next - cur;

      //perterb parameter      
      midParam = pParam->GetEstVal();
      hiParam = midParam  + dx;
      lowParam = midParam - dx;

      //avoid exceeding parameter limits
      if(hiParam > upr) 
      { 
         hiParam = lowParam; //move in opposite direction
         flipSign = true;     //must reverse sign when done
         dType = FD_FORWARD;  //only use forward difference
      } /* end if() */
      if(lowParam < lwr) 
      { 
         dType = FD_FORWARD; //only use forward difference
      }/* end if() */

      //adjust change in parameter
      switch(dType)
      {
         case(FD_FIT_CEN) :         
         case(FD_OUT_CEN) :
            hiParam  = midParam + (0.5 * dx);
            lowParam = midParam - (0.5 * dx);
            break;
         case(FD_PAR_CEN) :
            hiParam  = midParam + (0.5 * dx);
            lowParam = midParam - (0.5 * dx);

            //prepare parabolic matrix            
            m_ParaMat[0][2] = 1.00;
            m_ParaMat[1][2] = 1.00;
            m_ParaMat[2][2] = 1.00;
            m_ParaMat[0][0] = lowParam * lowParam;
            m_ParaMat[0][1] = lowParam;
            m_ParaMat[1][0] = midParam * midParam;
            m_ParaMat[1][1] = midParam;
            m_ParaMat[2][0] = hiParam * hiParam;
            m_ParaMat[2][1] = hiParam;
            //compute inverse
            MatInv(m_ParaMat, m_ParaInv, 3);
            break;         
         case(FD_FORWARD) :            
         default:
            if(flipSign == true) 
            { 
               flipSign = false; //clear flag, for next iteration
               dx *= -1.0;
            }/* end if() */
            break;         
      }/* end switch() */

      /*---------------------------------------
      Perform required model executions, saving 
      results into the appropriate model backup.
      -----------------------------------------*/      
      pParam->SetEstVal(hiParam);
      F = m_pModel->Execute();
      m_DiffCount++;
      m_pHiBkup->Store();

      /* ------------------------------------------
      Update the min. Jacobian, if necessary
      ------------------------------------------ */
      if(F < m_pMinJac[0])
      {
         //first element is objective function value
         m_pMinJac[0] = F;
         //next elements are the parameter settings
         m_pModel->GetParamGroupPtr()->ReadParams(&(m_pMinJac[1]));
         //final elements are the simulated observation values
         m_pModel->GetObsGroupPtr()->ReadObservations(&(m_pMinJac[1+m_NumParams]));
      }

      if(dType != FD_FORWARD)
      {
         pParam->SetEstVal(lowParam);
         F = m_pModel->Execute();
         m_DiffCount++;
         m_pLowBkup->Store();

         /* ------------------------------------------
         Update the min. Jacobian, if necessary
         ------------------------------------------ */
         if(F < m_pMinJac[0])
         {
            //first element is objective function value
            m_pMinJac[0] = F;
            //next elements are the parameter settings
            m_pModel->GetParamGroupPtr()->ReadParams(&(m_pMinJac[1]));
            //final elements are the simulated observation values
            m_pModel->GetObsGroupPtr()->ReadObservations(&(m_pMinJac[1+m_NumParams]));
         }/* end if() */
      }/* end if() */

      //compute change of each observation
      total_diff = 0.00;
      for(i = 0; i < m_NumObs; i++)
      {
         midObs = m_pMidBkup->GetObs(i, true, true);
         hiObs = m_pHiBkup->GetObs(i, true, true);

         if(dType != FD_FORWARD)
         {
            lowObs = m_pLowBkup->GetObs(i, true, true);
         }/* end if() */
         
         //computation method for derivative depends on dType
         switch(dType)
         {
            case(FD_OUT_CEN) : //outside central
               dObs = (hiObs - lowObs);
               diff = (dObs / dx);
               break;
            case(FD_PAR_CEN) : //parabolic central
               //fill obs. vector
               paraObs[0] = lowObs;
               paraObs[1] = midObs;
               paraObs[2] = hiObs;
               //compute coefficients
               VectMult(m_ParaInv, paraObs, paraCof, 3, 3);
               //compute derivative approx. (dy/dx = 2ax + b)
               diff = (2.00 * paraCof[0] * midParam) + paraCof[1];
               break;
            case(FD_FIT_CEN) : //Least-Squares best-fit central
               /*-------------------------------------------------------               
               The derivative of (y = bx + a) using Least-Sqaures is: 
                  dy/dx = b = (SSxy - SxSy)/(SSxx - (Sx)^2)
                  where, 
                  S = #points (3),
                  Sxy = (x1y1 + x2y2 + x3y3),
                  Sx = (x1 + x2 + x3),
                  Sy = (y1 + y2 + y3),
                  Sxx = (x1^2 + x2^2 + x3^2),
                  x1,x2,x3 -> parameters, 
                  y1,y2,y3 -> observations,
                  and equal variance has been given to each data point.
               -------------------------------------------------------*/       
               Sxy = ((lowObs*lowParam)+(midObs*midParam)+(hiObs*hiParam));
               Sx = (lowParam + midParam + hiParam);
               Sy = (lowObs + midObs + hiObs);
               Sxx=(lowParam*lowParam)+(midParam*midParam)+(hiParam*hiParam);
               diff = ((3.00 * Sxy) - (Sx * Sy)) / ((3.00 * Sxx) - (Sx * Sx));
               break;
            case(FD_FORWARD) : //forward difference
            default:
               dObs = hiObs - midObs;
               diff = (dObs / dx);
               break;         
         }/* end switch() */

         idx = (i*m_NumParams) + j;
         m_pMyBuf[idx] = diff;

         m_pJacob[i][j] = diff;
         m_pJacobT[j][i] = diff;
         total_diff += fabs(diff);
      }/* end for(observations) */
      m_pMidBkup->SemiRestore();

      //if sensitivity of all observations is ~0.00, retry using an alternative increment type
      if((total_diff <= NEARLY_ZERO) && (dIncType != FD_RANGE_REL) && (GetProgramType() != JACOBIAN_PROGRAM))
      {
         dIncType = FD_RANGE_REL;
         goto retry;
      }
   } /* end for(parameters) */   

   //gather results
   double * pTmp = new double[bufsize];
   for(i = 0; i < bufsize; i++) pTmp[i] = 0.00;
   MPI_Reduce(m_pMyBuf, pTmp, bufsize, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   //store the processors results into the Jacobian matrix
   for(j = 0; j < m_NumParams; j ++)
   {
      for(i = 0; i < m_NumObs; i++)
      {
         idx = (i*m_NumParams) + j;
         m_pJacob[i][j] = pTmp[idx];
         m_pJacobT[j][i] = pTmp[idx];         
      }
   }
   delete [] pTmp;

   //collect Minimum Jacobian data
   BcastMinJac();
}/* end EvalJacParallel() */

/******************************************************************************
EvalJacSuperMUSE()

Compute the Jacobian entries using SuperMUSE. This routine interfaces with the 
RepeatTasker SuperMUSE program, which assigns model evaluations to SuperMUSE 
clients on a first-come-first-served basis.
******************************************************************************/
void StatsClass::EvalJacSuperMUSE(void)
{  
   bool flipSign, bOk;
   double upr;    //upper bound of parameter
   double lwr;    //lower bound of parameter
   double F, dx, cur, next;     //parameter perterbation
   int i;
   int j;
   int task;
   //double Flwr, Fupr, range;
   double diff; //derivative approximation, to be stored in Jacobian
   double dObs;  //forward and outside central diff. numerator
   double lowParam; //central diff. LHS
   double hiParam; //central diff. RHS and forward diff. RHS
   double midParam; //central diff middle parameter value, forward diff. LHS
   double lowObs; //central diff. LHS
   double hiObs; //central diff. RHS and forward diff.
   double midObs; //current observation value
   double paraCof[3];    //coeff. of parabolic solution
   double paraObs[3];    //obs. vals. of parabolic solution
   double Sxy, Sx, Sy, Sxx; //central diff. best-fit parameters   
   ParameterGroup * pParamGroup;
   ParameterABC * pParam;
   double * pPoint;
   FiniteDiffType dType;
   FiniteDiffIncType dIncType;
   double total_diff;
   double negOne = -1.00;
   SuperMUSE * pSMUSE = GetSuperMusePtr();
   
   pParamGroup = m_pModel->GetParamGroupPtr();

   m_pMidBkup->Store();
   flipSign = false;

   /* ------------------------------------------
   Intialize the min. Jacobian
   ------------------------------------------ */
   //first element is objective function value
   m_pMinJac[0] = m_pModel->GetObjFuncVal();
   //next elements are the parameter settings
   m_pModel->GetParamGroupPtr()->ReadParams(&(m_pMinJac[1]));
   //final elements are the simulated observation values
   m_pModel->GetObsGroupPtr()->ReadObservations(&(m_pMinJac[1+m_NumParams]));

   /*------------------------------------------------------
   Perterb each parameter, one at a time, and assemble a 
   list of the required model evaluations.
   -------------------------------------------------------*/
   for(j = 0; j < m_NumParams; j++)
   {  
      dType = m_DiffType;
      dIncType = m_DiffIncType;
      pParam = pParamGroup->GetParamPtr(j);

      cur = pParam->GetEstVal();
      upr = pParam->GetUprBnd();
      lwr = pParam->GetLwrBnd();
      if(dIncType == FD_OPTIMAL)
      {
         NEW_PRINT("double", m_NumParams);
         pPoint = new double[m_NumParams];
         MEM_CHECK(pPoint);

         pParamGroup->ReadParams(pPoint);
         dx = CalcOptimalStepSize(j, pPoint);
         delete [] pPoint;      
      }
      else if(dIncType == FD_RANGE_REL)
      {
         dx = fabs(m_pDiffInc[j]*(upr-lwr));
      }
      else if(dIncType == FD_VALUE_REL)
      {
         dx = MyMax(fabs(m_pDiffInc[j]*cur), m_MinInc);
      }
      else if(dIncType == FD_ABSOLUTE)
      {
         dx = fabs(m_pDiffInc[j]);
      }
      else //default to range-relative increment
      {
         dx = fabs(m_pDiffInc[j]*(upr-lwr));
      }
      //trick from NR in C
      next = cur + dx;
      dx = next - cur;

      //perterb parameter      
      midParam = pParam->GetEstVal();
      hiParam = midParam  + dx;
      lowParam = midParam - dx;
 
      //avoid exceeding parameter limits
      if(hiParam > upr) 
      { 
         hiParam = lowParam; //move in opposite direction
         flipSign = true;     //must reverse sign when done
         dType = FD_FORWARD;  //only use forward difference
      } /* end if() */
      if(lowParam < lwr) 
      { 
         dType = FD_FORWARD; //only use forward difference
      }/* end if() */

      //adjust change in parameter
      switch(dType)
      {
         case(FD_FIT_CEN) :         
         case(FD_OUT_CEN) :
         case(FD_PAR_CEN) :
            hiParam  = midParam + (0.5 * dx);
            lowParam = midParam - (0.5 * dx);
            break;
         case(FD_FORWARD) :            
         default:
            if(flipSign == true) 
            { 
               flipSign = false; //clear flag, for next iteration
               dx *= -1.0;
            }/* end if() */
            break;         
      }/* end switch() */

      /* --------------------------------------
      Save parameter vars (will need them again
      when SuperMUSE finishes its work.
      -------------------------------------- */
      m_pDType[j] = dType;
      m_pDx[j]    = dx;
      m_pMid[j]   = midParam;
      m_pHi[j]    = hiParam;
      m_pLow[j]   = lowParam;

      /*---------------------------------------
      Store required model executions as items
      in the SuperMUSE task list. The convention
      is the high parameter first and the low
      parameter (if needed) second.
      -----------------------------------------*/      
      pParam->SetEstVal(hiParam);
      //pass group to supermuse
      pSMUSE->WriteTask(pParamGroup);

      if(dType != FD_FORWARD)
      {
         pParam->SetEstVal(lowParam);
        //pass group to supermuse
        pSMUSE->WriteTask(pParamGroup);
      }
      
      m_pMidBkup->SemiRestore();
   }/* end for() */

   /* -------------------------------------------------
   Now that a list of required evaluations has been 
   assembled, signal SuperMUSE to begin processing the
   job and then wait for SuperMUSE to complete.
   ------------------------------------------------- */
   //Finish task file (this will cause RepeatTasker to begin processing the job)
   pSMUSE->FinishTaskFile();

   //wait for SuperMUSE to report back (via the success or error files)
   bOk = pSMUSE->WaitForTasker();

   if(bOk == false) //SuperMUSE failed
   {
      LogError(ERR_SMUSE, "Reverting to serial execution.");
      DisableSuperMUSE();
      CalcJacobian(&negOne);
   }
   else //SuperMUSE was successful
   {     
      /* -------------------------------------------------
      Now that SuperMUSE has completed all of the required
      evaluations, compute the Jacobian entries.
      ------------------------------------------------- */
      task = 0;
      for(j = 0; j < m_NumParams; j++)
      {  
         pParam = pParamGroup->GetParamPtr(j);

         //retrieved stored vars for the given parameter
         dType    = m_pDType[j];
         dx       = m_pDx[j];
         midParam = m_pMid[j];
         hiParam  = m_pHi[j];
         lowParam = m_pLow[j];
 
         //prepare parabolic matrix, if needed
         if(dType == FD_PAR_CEN)
         {
            m_ParaMat[0][2] = 1.00;
            m_ParaMat[1][2] = 1.00;
            m_ParaMat[2][2] = 1.00;
            m_ParaMat[0][0] = lowParam * lowParam;
            m_ParaMat[0][1] = lowParam;
            m_ParaMat[1][0] = midParam * midParam;
            m_ParaMat[1][1] = midParam;
            m_ParaMat[2][0] = hiParam * hiParam;
            m_ParaMat[2][1] = hiParam;
            //compute inverse
            MatInv(m_ParaMat, m_ParaInv, 3);
         }/* end if() */

         /*---------------------------------------
         Extract results from the model executions, 
         saving results into the appropriate model 
         backup.
         -----------------------------------------*/      
         pParam->SetEstVal(hiParam);
         F = pSMUSE->GatherResult(task);
         task++;
         m_DiffCount++;
         m_pHiBkup->Store();

         /* ------------------------------------------
         Update the min. Jacobian, if necessary
         ------------------------------------------ */
         if(F < m_pMinJac[0])
         {
            //first element is objective function value
            m_pMinJac[0] = F;
            //next elements are the parameter settings
            m_pModel->GetParamGroupPtr()->ReadParams(&(m_pMinJac[1]));
            //final elements are the simulated observation values
            m_pModel->GetObsGroupPtr()->ReadObservations(&(m_pMinJac[1+m_NumParams]));
         }/* end if() */

         if(dType != FD_FORWARD)
         {
            pParam->SetEstVal(lowParam);
            F = pSMUSE->GatherResult(task);
            task++;
            m_DiffCount++;
            m_pLowBkup->Store();

            /* ------------------------------------------
            Update the min. Jacobian, if necessary
            ------------------------------------------ */
            if(F < m_pMinJac[0])
            {
               //first element is objective function value
               m_pMinJac[0] = F;
               //next elements are the parameter settings
               m_pModel->GetParamGroupPtr()->ReadParams(&(m_pMinJac[1]));
               //final elements are the simulated observation values
               m_pModel->GetObsGroupPtr()->ReadObservations(&(m_pMinJac[1+m_NumParams]));
            }/* end if() */
         }/* end if() */

         //compute change of each observation
         total_diff = 0.00;
         for(i = 0; i < m_NumObs; i++)
         {
            midObs = m_pMidBkup->GetObs(i, true, true);
            hiObs = m_pHiBkup->GetObs(i, true, true);

            if(dType != FD_FORWARD)
            {
               lowObs = m_pLowBkup->GetObs(i, true, true);
            }/* end if() */
         
            //computation method for derivative depends on dType
            switch(dType)
            {
               case(FD_OUT_CEN) : //outside central
                 dObs = (hiObs - lowObs);
                  diff = (dObs / dx);
                  break;
               case(FD_PAR_CEN) : //parabolic central
                  //fill obs. vector
                  paraObs[0] = lowObs;
                  paraObs[1] = midObs;
                  paraObs[2] = hiObs;
                  //compute coefficients
                  VectMult(m_ParaInv, paraObs, paraCof, 3, 3);
                  //compute derivative approx. (dy/dx = 2ax + b)
                  diff = (2.00 * paraCof[0] * midParam) + paraCof[1];
                  break;
               case(FD_FIT_CEN) : //Least-Squares best-fit central
                  /*-------------------------------------------------------               
                  The derivative of (y = bx + a) using Least-Sqaures is: 
                     dy/dx = b = (SSxy - SxSy)/(SSxx - (Sx)^2)
                     where, 
                     S = #points (3),
                     Sxy = (x1y1 + x2y2 + x3y3),
                     Sx = (x1 + x2 + x3),
                     Sy = (y1 + y2 + y3),
                     Sxx = (x1^2 + x2^2 + x3^2),
                     x1,x2,x3 -> parameters, 
                     y1,y2,y3 -> observations,
                     and equal variance has been given to each data point.
                  -------------------------------------------------------*/       
                  Sxy = ((lowObs*lowParam)+(midObs*midParam)+(hiObs*hiParam));
                  Sx = (lowParam + midParam + hiParam);
                  Sy = (lowObs + midObs + hiObs);
                  Sxx=(lowParam*lowParam)+(midParam*midParam)+(hiParam*hiParam);
                  diff = ((3.00 * Sxy) - (Sx * Sy)) / ((3.00 * Sxx) - (Sx * Sx));
                  break;
               case(FD_FORWARD) : //forward difference
               default:
                  dObs = hiObs - midObs;
                  diff = (dObs / dx);
                  break;         
            }/* end switch() */
         
            m_pJacob[i][j] = diff;
            m_pJacobT[j][i] = diff;         
            total_diff += fabs(diff);
         }/* end for(observations) */
      
         m_pMidBkup->SemiRestore();
      } /* end for(parameters) */   
   }/* end if() */
}/* end EvalJacSuperMUSE() */

/******************************************************************************
EvalJacSerial()

Evaluate the Jacobian in serial.
******************************************************************************/
void StatsClass::EvalJacSerial(double * pBestSavedF)
{
   bool flipSign;
   double upr;    //upper bound of parameter
   double lwr;    //lower bound of parameter
   double F, dx, cur, next;     //parameter perterbation
   int i;
   int j;
   //double Flwr, Fupr, range;
   double diff; //derivative approximation, to be stored in Jacobian
   double total_diff; 
   double dObs;  //forward and outside central diff. numerator
   double lowParam; //central diff. LHS
   double hiParam; //central diff. RHS and forward diff. RHS
   double midParam; //central diff middle parameter value, forward diff. LHS
   double lowObs; //central diff. LHS
   double hiObs; //central diff. RHS and forward diff.
   double midObs; //current observation value
   double paraCof[3];    //coeff. of parabolic solution
   double paraObs[3];    //obs. vals. of parabolic solution
   double Sxy, Sx, Sy, Sxx; //central diff. best-fit parameters   
   ParameterGroup * pParamGroup;
   ParameterABC * pParam;
   double * pPoint;
   FiniteDiffType dType;
   FiniteDiffIncType dIncType;
      
   pParamGroup = m_pModel->GetParamGroupPtr();

   m_pMidBkup->Store();
   flipSign = false;

   /* ------------------------------------------
   Intialize the min. Jacobian
   ------------------------------------------ */
   //first element is objective function value
   m_pMinJac[0] = m_pModel->GetObjFuncVal();
   //next elements are the parameter settings
   m_pModel->GetParamGroupPtr()->ReadParams(&(m_pMinJac[1]));
   //final elements are the simulated observation values
   m_pModel->GetObsGroupPtr()->ReadObservations(&(m_pMinJac[1+m_NumParams]));

   /*------------------------------------------------------
   Perterb each parameter, one at a time, and compute 
   resulting change in observations and store in the 
   Jacobian matrix.
   -------------------------------------------------------*/
   for(j = 0; j < m_NumParams; j++)
   {  
      dType = m_DiffType;
      dIncType = m_DiffIncType;

//if FD calculation is ~0.00, retry using an alternative increment type
retry:
      pParam = pParamGroup->GetParamPtr(j);

      cur = pParam->GetEstVal();
      upr = pParam->GetUprBnd();
      lwr = pParam->GetLwrBnd();
      if(dIncType == FD_OPTIMAL)
      {
         NEW_PRINT("double", m_NumParams);
         pPoint = new double[m_NumParams];
         MEM_CHECK(pPoint);

         pParamGroup->ReadParams(pPoint);
         dx = CalcOptimalStepSize(j, pPoint);
         delete [] pPoint;      
      }
      else if(dIncType == FD_RANGE_REL)
      {
         dx = fabs(m_pDiffInc[j]*(upr-lwr));
      }
      else if(dIncType == FD_VALUE_REL)
      {
         dx = MyMax(fabs(m_pDiffInc[j]*cur), m_MinInc);
      }
      else if(dIncType == FD_ABSOLUTE)
      {
         dx = fabs(m_pDiffInc[j]);
      }
      else //default to range-relative increment
      {
         dx = fabs(m_pDiffInc[j]*(upr-lwr));
      }

      //trick from NR in C
      next = cur + dx;
      dx = next - cur;

      //perterb parameter      
      midParam = pParam->GetEstVal();
      hiParam = midParam  + dx;
      lowParam = midParam - dx;
 
      //avoid exceeding parameter limits
      if(hiParam > upr) 
      { 
         hiParam = lowParam; //move in opposite direction
         flipSign = true;     //must reverse sign when done
         dType = FD_FORWARD;  //only use forward difference
      } /* end if() */
      if(lowParam < lwr) 
      { 
         dType = FD_FORWARD; //only use forward difference
      }/* end if() */

      //adjust change in parameter
      switch(dType)
      {
         case(FD_FIT_CEN) :         
         case(FD_OUT_CEN) :
            hiParam  = midParam + (0.5 * dx);
            lowParam = midParam - (0.5 * dx);
            break;
         case(FD_PAR_CEN) :
            hiParam  = midParam + (0.5 * dx);
            lowParam = midParam - (0.5 * dx);

            //prepare parabolic matrix            
            m_ParaMat[0][2] = 1.00;
            m_ParaMat[1][2] = 1.00;
            m_ParaMat[2][2] = 1.00;
            m_ParaMat[0][0] = lowParam * lowParam;
            m_ParaMat[0][1] = lowParam;
            m_ParaMat[1][0] = midParam * midParam;
            m_ParaMat[1][1] = midParam;
            m_ParaMat[2][0] = hiParam * hiParam;
            m_ParaMat[2][1] = hiParam;
            //compute inverse
            MatInv(m_ParaMat, m_ParaInv, 3);
            break;         
         case(FD_FORWARD) :            
         default:
            if(flipSign == true) 
            { 
               flipSign = false; //clear flag, for next iteration
               dx *= -1.0;
            }/* end if() */
            break;         
      }/* end switch() */

      /*---------------------------------------
      Perform required model executions, saving 
      results into the appropriate model backup.
      -----------------------------------------*/      
      pParam->SetEstVal(hiParam);
      F = m_pModel->Execute();
      if(F < (*pBestSavedF))
      {
         m_pModel->SaveBest(0);
         (*pBestSavedF) = F;
      }
      m_DiffCount++;
      m_pHiBkup->Store();

      /* ------------------------------------------
      Update the min. Jacobian, if necessary
      ------------------------------------------ */
      if(F < m_pMinJac[0])
      {
         //first element is objective function value
         m_pMinJac[0] = F;
         //next elements are the parameter settings
         m_pModel->GetParamGroupPtr()->ReadParams(&(m_pMinJac[1]));
         //final elements are the simulated observation values
         m_pModel->GetObsGroupPtr()->ReadObservations(&(m_pMinJac[1+m_NumParams]));
      }/* end if() */

      if(dType != FD_FORWARD)
      {
         pParam->SetEstVal(lowParam);
         F = m_pModel->Execute();
         if(F < (*pBestSavedF))
         {
            m_pModel->SaveBest(0);
            (*pBestSavedF) = F;
         }

         m_DiffCount++;
         m_pLowBkup->Store();

         /* ------------------------------------------
         Update the min. Jacobian, if necessary
         ------------------------------------------ */
         if(F < m_pMinJac[0])
         {
            //first element is objective function value
            m_pMinJac[0] = F;
            //next elements are the parameter settings
            m_pModel->GetParamGroupPtr()->ReadParams(&(m_pMinJac[1]));
            //final elements are the simulated observation values
            m_pModel->GetObsGroupPtr()->ReadObservations(&(m_pMinJac[1+m_NumParams]));
         }/* end if() */
      }/* end if() */

      //compute change of each observation
      total_diff = 0.00;
      for(i = 0; i < m_NumObs; i++)
      {
         midObs = m_pMidBkup->GetObs(i, true, true);
         hiObs = m_pHiBkup->GetObs(i, true, true);

         if(dType != FD_FORWARD)
         {
            lowObs = m_pLowBkup->GetObs(i, true, true);
         }/* end if() */
         
         //computation method for derivative depends on dType
         switch(dType)
         {
            case(FD_OUT_CEN) : //outside central
               dObs = (hiObs - lowObs);
               diff = (dObs / dx);
               break;
            case(FD_PAR_CEN) : //parabolic central
               //fill obs. vector
               paraObs[0] = lowObs;
               paraObs[1] = midObs;
               paraObs[2] = hiObs;
               //compute coefficients
               VectMult(m_ParaInv, paraObs, paraCof, 3, 3);
               //compute derivative approx. (dy/dx = 2ax + b)
               diff = (2.00 * paraCof[0] * midParam) + paraCof[1];
               break;
            case(FD_FIT_CEN) : //Least-Squares best-fit central
               /*-------------------------------------------------------               
               The derivative of (y = bx + a) using Least-Sqaures is: 
                  dy/dx = b = (SSxy - SxSy)/(SSxx - (Sx)^2)
                  where, 
                  S = #points (3),
                  Sxy = (x1y1 + x2y2 + x3y3),
                  Sx = (x1 + x2 + x3),
                  Sy = (y1 + y2 + y3),
                  Sxx = (x1^2 + x2^2 + x3^2),
                  x1,x2,x3 -> parameters, 
                  y1,y2,y3 -> observations,
                  and equal variance has been given to each data point.
               -------------------------------------------------------*/       
               Sxy = ((lowObs*lowParam)+(midObs*midParam)+(hiObs*hiParam));
               Sx = (lowParam + midParam + hiParam);
               Sy = (lowObs + midObs + hiObs);
               Sxx=(lowParam*lowParam)+(midParam*midParam)+(hiParam*hiParam);
               diff = ((3.00 * Sxy) - (Sx * Sy)) / ((3.00 * Sxx) - (Sx * Sx));
               break;
            case(FD_FORWARD) : //forward difference
            default:
               dObs = hiObs - midObs;
               diff = (dObs / dx);
               //printf("%E,%E,%E,%E\n", midObs, hiObs, midParam, hiParam);
               break;         
         }/* end switch() */
         
         m_pJacob[i][j] = diff;
         m_pJacobT[j][i] = diff;
         total_diff += fabs(diff);
      }/* end for(observations) */      

      //compute change of each prediction, if applicable
      int nrv = 0;
      if (m_pPredictions != NULL) nrv = m_pPredictions->GetNumRespVars();
      for(i = 0; i < nrv; i++)
      {
         midObs = m_pMidBkup->GetPred(i);
         hiObs = m_pHiBkup->GetPred(i);

         if(dType != FD_FORWARD)
         {
            lowObs = m_pLowBkup->GetPred(i);
         }/* end if() */
         
         //computation method for derivative depends on dType
         switch(dType)
         {
            case(FD_OUT_CEN) : //outside central
               dObs = (hiObs - lowObs);
               diff = (dObs / dx);
               break;
            case(FD_PAR_CEN) : //parabolic central
               //fill obs. vector
               paraObs[0] = lowObs;
               paraObs[1] = midObs;
               paraObs[2] = hiObs;
               //compute coefficients
               VectMult(m_ParaInv, paraObs, paraCof, 3, 3);
               //compute derivative approx. (dy/dx = 2ax + b)
               diff = (2.00 * paraCof[0] * midParam) + paraCof[1];
               break;
            case(FD_FIT_CEN) : //Least-Squares best-fit central
               /*-------------------------------------------------------               
               The derivative of (y = bx + a) using Least-Sqaures is: 
                  dy/dx = b = (SSxy - SxSy)/(SSxx - (Sx)^2)
                  where, 
                  S = #points (3),
                  Sxy = (x1y1 + x2y2 + x3y3),
                  Sx = (x1 + x2 + x3),
                  Sy = (y1 + y2 + y3),
                  Sxx = (x1^2 + x2^2 + x3^2),
                  x1,x2,x3 -> parameters, 
                  y1,y2,y3 -> observations,
                  and equal variance has been given to each data point.
               -------------------------------------------------------*/       
               Sxy = ((lowObs*lowParam)+(midObs*midParam)+(hiObs*hiParam));
               Sx = (lowParam + midParam + hiParam);
               Sy = (lowObs + midObs + hiObs);
               Sxx=(lowParam*lowParam)+(midParam*midParam)+(hiParam*hiParam);
               diff = ((3.00 * Sxy) - (Sx * Sy)) / ((3.00 * Sxx) - (Sx * Sx));
               break;
            case(FD_FORWARD) : //forward difference
            default:
               dObs = hiObs - midObs;
               diff = (dObs / dx);
               //printf("%E,%E,%E,%E\n", midObs, hiObs, midParam, hiParam);
               break;         
         }/* end switch() */
         
         m_pJacPred[i][j] = diff;
      }/* end for(predictions) */

      m_pMidBkup->SemiRestore();

      //if sensivity of all observations ~0.00, retry using an alternative increment type
      if((total_diff <= NEARLY_ZERO) && (dIncType != FD_RANGE_REL) && (GetProgramType() != JACOBIAN_PROGRAM))
      {
         dIncType = FD_RANGE_REL;
         goto retry;
      }
   } /* end for(parameters) */   
}/* end EvalJacSerial() */

/******************************************************************************
CalcOptimalStepSize()

Calculate the optimal step size using the equations (4) and (5) from Yager, 2004.
"Effects of Model Sensitivty and Nonlinearity on Nonlinear Regression of Ground-
Water Flow".
******************************************************************************/
double StatsClass::CalcOptimalStepSize(int idx, double * params)
{
   int numTries, maxTries;
   double delta = 1.00;
   double db, old_db, eps, tmp;
   double Fmid, Fupr, Flwr;
   double bMid;
   double Sjj;

   bMid = params[idx];
   Fmid = m_pModel->Execute();
   m_StepCount++;

   eps = 1E-6;
   db = old_db = 2.00*sqrt(eps)*fabs(bMid);
   numTries = 0;
   maxTries = 5;
   /*-----------------------------------------------------
   Try and iterate until convergence; if anything unusual 
   happens, just use the original guess.
   ------------------------------------------------------*/   
   while(delta > eps)
   {
      if(numTries >= maxTries)
      {
         db = 2.00*sqrt(eps)*fabs(bMid); 
         break;
      }
      numTries++;

      params[idx] = bMid + db;
      m_pModel->GetParamGroupPtr()->WriteParams(params);
      Fupr = m_pModel->Execute();
      m_StepCount++;

      params[idx] = bMid - db;
      m_pModel->GetParamGroupPtr()->WriteParams(params);
      Flwr = m_pModel->Execute();
      m_StepCount++;

      Sjj = (Fupr - 2.00*Fmid + Flwr) / (db*db);
      if(Sjj == 0.00)
      { 
         db = 2.00*sqrt(eps)*fabs(bMid); 
         break;
      }
      tmp = (4.00*eps*Fmid)/Sjj;
      if(tmp <= 0.00)
      { 
         db = 2.00*sqrt(eps)*fabs(bMid); 
         break;
      }
      db = sqrt(fabs(tmp));      
      delta = fabs(db-old_db);
      old_db = db;
   }/* end if() */

   params[idx] = bMid;
   m_pModel->GetParamGroupPtr()->WriteParams(params);

   return db;
}/* end CalcOptimalStepSize() */

/******************************************************************************
GetJacobT()

Retrieve transpose of the Jacobian.
******************************************************************************/
Unchangeable2DArray StatsClass::GetJacobT(void)
{
   return(m_pJacobT);
}/* end GetJacobT() */

/******************************************************************************
GetJacobUW()

Retrieve the unweighted Jacobian.
******************************************************************************/
Unchangeable2DArray StatsClass::GetJacobUW(void)
{
   return(m_pJacobUW);
}/* end GetJacobUW() */

/******************************************************************************
CalcNormal()

Calculate what is known as the 'normal' regression matrix, also referred to as 
the Fisher information matrix:
   (J^T)*Q*J.
******************************************************************************/
Unchangeable2DArray StatsClass::CalcNormal(void)
{
   int n, p;

   n = m_NumObs;
   p = m_NumParams;
   if(m_bAdjustedJac == true)
   {
      n -= m_NumHeldObs;
      p -= m_NumHeldParams;
   }
   MatMult(m_pJacobT, m_pJacob, m_pNormal, p, n, p);

   return m_pNormal;
}/* end CalcNormal() */

/******************************************************************************
CalcStats()

Calculate the statistics requested by the user.
******************************************************************************/
void StatsClass::CalcStats(void)
{
   if(m_bNoStats) return;

   int i, j, id, n, p, nrv;
   double negOne = -1.00;

   //read in predictions and store as response variables
   nrv = 0;
   if(m_pPredictions != NULL)
   {
      nrv = m_pPredictions->GetNumRespVars();
      for(i = 0; i < nrv; i++){delete [] m_pJacPred[i];}
      delete [] m_pJacPred;
      delete m_pPredictions;
   }
   NEW_PRINT("ResponseVarGroup", 1);
   m_pPredictions = new ResponseVarGroup("Predictions");
   MEM_CHECK(m_pPredictions);

   //must compute parameter variance if interested in prediction stats.
   nrv = m_pPredictions->GetNumRespVars();
   if(nrv > 0)
   {
      //alert backups to the presence of predictions
      m_pMidBkup->SetResponseVarGroup(m_pPredictions);
      m_pLowBkup->SetResponseVarGroup(m_pPredictions);
      m_pHiBkup->SetResponseVarGroup(m_pPredictions);

      //allocate space for partial derivatives of the predictions
      NEW_PRINT("double *", nrv);
      m_pJacPred = new double*[nrv];
      for(i = 0; i < nrv; i++)
      {
         NEW_PRINT("double", m_NumParams);
         m_pJacPred[i] = new double[m_NumParams];
      }
      MEM_CHECK(m_pJacPred[i-1]);

      delete [] m_Pred;
      delete [] m_PredSD;
      delete [] m_PredCIupr;
      delete [] m_PredCIlwr;
      m_Pred   = new double[nrv];
      m_PredSD = new double[nrv];
      m_PredCIupr = new double[nrv];
      m_PredCIlwr = new double[nrv];
      m_StdErrFlag = true;
   }

   MPI_Comm_rank(MPI_COMM_WORLD, &id);

   m_Phi = m_pModel->Execute();
   m_StatsCount++;

   CalcJacobian(m_bOkToHoldParams, m_bOkToHoldObs, &negOne); //compute Jacobian, possibly in parallel

   if(id == 0)
   {
      n = m_NumObs - m_NumHeldObs;
      p = m_NumParams - m_NumHeldParams;
      AdjustJacobian();
      CalcResiduals();
      AdjustResiduals();
      CalcNormal();
      m_bInv = MatInv(m_pNormal, m_pInvNormal, p);
      if(m_bInv == false)
      {
         LogError(ERR_SING_MAT, "Could not invert normal matrix (J^T Q J). As such,");
         LogError(ERR_CONTINUE, "the following statistics will not be computed:");
         LogError(ERR_CONTINUE, "  (1)  parameter variance-covariance,");
         LogError(ERR_CONTINUE, "  (2)  parameter standard error, ");
         LogError(ERR_CONTINUE, "  (3)  parameter correlation, ");
         LogError(ERR_CONTINUE, "  (4)  linear confidence intervals, ");
         LogError(ERR_CONTINUE, "  (5)  the 'volume ratio' measure, ");
         LogError(ERR_CONTINUE, "  (6)  influential observations, ");
         LogError(ERR_CONTINUE, "  (7)  linearity measures, and ");
         LogError(ERR_CONTINUE, "  (8)  prediction statistics.");
         m_StdErrFlag = m_CorrCoefFlag = false;
         m_CIflag = false;                  
         m_CooksFlag = m_DfbetasFlag = false;
         m_BealeFlag = m_LinssenFlag = false;
      }
      CalcWeightedRy();
      CalcRawRy();
      m_Phi = AdjustObjFunc(m_Phi);
      m_Variance = m_Phi / (n - p);

      if(m_StdErrFlag == true)
      {
         for(i = 0; i < p; i++){
            for(j = 0; j < p; j++){
               m_pCovar[i][j] = m_pInvNormal[i][j] * m_Variance;
            }/* end for() */
         }/* end for() */
      }

      if(m_NormPlotFlag == true){CalcNormProbPlot();}   
      if(m_BestBoxCoxFlag == true){CalcBestBoxCox();}   

      if(m_RunsTestFlag == true){
         m_Runs.bSuccess = RunsTest(m_pResid, n, &(m_Runs.pos), 
                                    &(m_Runs.neg), &(m_Runs.runs),
                                    &(m_Runs.clwr), &(m_Runs.cupr)); }      
      if(m_AutorunFunctionFlag == true){
         AutorunFunctionTest(m_pResid, n, &(m_AR.r1), &(m_AR.var), &(m_AR.vpx), 
                             &(m_AR.med), &(m_AR.sur), &(m_AR.def), &(m_AR.n1),
                             &(m_AR.clwr), &(m_AR.cupr)); }

      if(m_MMRI_Flag == true) CalcMMRI(m_bInv);

      if(m_CIflag == true){ CalcCI();}

      if((m_BealeFlag == true) || (m_LinssenFlag == true)){CalcBealeAndLinssen();}

      if(m_CooksFlag == true){CalcCooksD();}

      if(m_DfbetasFlag == true){CalcDFBETAS();}

      if(m_SensFlag == true){CalcSensitivities();}

      if(m_pPredictions->GetNumRespVars() > 0){CalcPredictions(m_bInv, m_pCovar, m_NumParams);}
   }/* end if(master only) */
}/* end CalcStats() */

/******************************************************************************
CalcBestBoxCox()

Determine the optimal lambda value for a BoxCox transformation that conforms
the residuals to satisfy assumption of normality. This sets up and runs a 
separate Ostrich optimization on the internaa BoxCoxModel() objective function.
******************************************************************************/
void StatsClass::CalcBestBoxCox(void)
{
   int id;
   char cmd[DEF_STR_SZ];
   char * line;
   int max_line_size;
   m_BestBoxCoxVal = 0.00;

   MPI_Comm_rank(MPI_COMM_WORLD, &id);

   //only one proc. needs to perform calculation
   if(id != 0)
   {
      m_BestBoxCoxFlag = false;
      return; 
   }
   
   //Create BoxCoxIn.tpl
   #ifdef WIN32
   system("mkdir BoxCoxModel 2> NUL");
   FILE * pTpl = fopen(".\\BoxCoxModel\\BoxCoxIn.tpl", "w");
   #else
   system("mkdir BoxCoxModel 2> /dev/null");
   FILE * pTpl = fopen("./BoxCoxModel/BoxCoxIn.tpl", "w");
   #endif

   if(pTpl == NULL)
   {
      LogError(ERR_FILE_IO, "Unable to create BoxCoxIn.tpl! Best BoxCox lamba value not computed.");
      return;
   }/* end if() */

   fprintf(pTpl, "LAMBDA=lambda\nNUM_DATA_POINTS=%d\n", m_NumObs);
   ObservationGroup * pG;
   Observation * pO;
   int i;
   double x,y,w;
   pG = m_pModel->GetObsGroupPtr();
   if(pG != NULL)
   {
      for(i = 0; i < m_NumObs; i++)
      {
        pO=pG->GetObsPtr(i);
        if(pO != NULL)
        {
            w=GetObsWeight(pO);
            x=pO->GetMeasuredVal(false,false);
            y=pO->GetComputedVal(false,false);
            fprintf(pTpl, "%E\t%E\t%E\n", x, y, w);
        }
        else
        {
            printf("Observation #%d is NULL!\n", i);
        }
      }/* end for() */
   }
   else
   {
      printf("Observation Group is NULL!\n");
   }   
   fclose(pTpl);

   //create input file
   #ifdef WIN32
   FILE * pIn = fopen(".\\BoxCoxModel\\ostIn.txt", "w");
   #else
   FILE * pIn = fopen("./BoxCoxModel/ostIn.txt", "w");
   #endif

   if(pIn == NULL)
   {
      LogError(ERR_FILE_IO, "Unable to create ostIn.txt! Best BoxCox lamba value not computed.");
      return;
   }/* end if() */
fprintf(pIn, "ProgramType Powell\n\
ObjectiveFunction GCOP\n\
\n\
ModelSubdir    .\n\
NumDigitsOfPrecision 16\n\
\n\
BeginFilePairs\n\
BoxCoxIn.tpl ; BoxCoxIn.txt\n\
EndFilePairs\n\
\n\
CheckSensitivities no\n\
ModelExecutable BoxCox()\n\
\n\
BeginParams\n\
lambda 1 -3 +3 none none none\n\
EndParams\n\
\n\
BeginResponseVars\n\
F(x)    BoxCoxOut.txt ; ObjFunc   0   2   '='\n\
EndResponseVars\n\
\n\
BeginGCOP\n\
CostFunction F(x)\n\
PenaltyFunction APM\n\
EndGCOP\n\
\n\
BeginConstraints\n\
EndConstraints\n\
\n\
BeginPowellAlg\n\
ConvergenceVal 1E-10\n\
MaxIterations 200\n\
EndPowellAlg\n\
\n\
Begin1dSearch\n\
1dSearchConvergeVal 1.000000E-006\n\
1dSearchMethod Brent\n\
End1dSearch\n");
   fclose(pIn);
 
   //run Ostrich
   char OstExe[DEF_STR_SZ];
   strcpy(OstExe, GetOstExePath());
   MyStrRep(OstExe, "IsoFit", "Ostrich");
   MyStrRep(OstExe, "OstrichMPI", "Ostrich");
   MyStrRep(OstExe, "OstrichFMPI", "Ostrich");
   #ifdef WIN32
      sprintf(cmd, "cd BoxCoxModel & %s > NUL & cd ..", OstExe);
   #else
      sprintf(cmd, "cd BoxCoxModel; %s > /dev/null; cd ..", OstExe);
   #endif
   //printf("Launching Ostrich using the following command |%s|\n", cmd);
   system(cmd);

   //retrieve result   
   #ifdef WIN32
   max_line_size = GetMaxLineSizeInFile(".\\BoxCoxModel\\OstOutput0.txt");
   FILE * pOut = fopen(".\\BoxCoxModel\\OstOutput0.txt", "r");   
   #else
   max_line_size = GetMaxLineSizeInFile((char *)"./BoxCoxModel/OstOutput0.txt");
   FILE * pOut = fopen("./BoxCoxModel/OstOutput0.txt", "r");
   #endif
   line = new char[max_line_size];
   if(pOut == NULL)
   {
      LogError(ERR_FILE_IO, "Unable to open OstOutput0.txt! Best BoxCox lamba value not computed.");
      return;
   }/* end if() */
   const char * pTok = "lambda             : ";
   line[0] = NULLSTR;
   while(!feof(pOut))
   {
      fgets(line, max_line_size, pOut);
      if(strncmp(line, pTok, strlen(pTok)) == 0)
      {
         sscanf(line, "lambda             :  %lf\n", &m_BestBoxCoxVal);
         break;
      }/* end if() */
   }/* end while() */
   delete [] line;

   //return to working dir
   fclose(pOut);
}/* end CalcBestBoxCox() */

/******************************************************************************
CalcCI()

Calculate linear confidence interval of the estimated parameters, using 
students t-distribution. Also calculates the joint confidence ellipsoid - 
confidence interval volume ratio and the corresponding percentage point of the
equivalent confidence ellipsoid (see Draper & Smith, "Applied Regression 
Analysis, third edition", page 144-145 for details).
******************************************************************************/
void StatsClass::CalcCI(void)
{
   double ** coeff; //correlation coefficient matrix
   double stdErr; //standard error of parameter
   double tStat;  //t-statistic
   double alpha;  //1.00 - (m_CIpct/100)
   double p;      //1-alpha/2
   double est;    //estimated value
   double tmp;
   ParameterGroup * pGroup;
   ParameterABC * pParam;
   int v;
   int i, j; //parameter loop variable
   int obs, params;

   obs = m_NumObs - m_NumHeldObs;
   params = m_NumParams - m_NumHeldParams;

   pGroup = m_pModel->GetParamGroupPtr();

   alpha = 1.00 - (m_CIpct/100.00);
   p = 1.00 - (alpha / 2.00);
   tStat = StudentInvCDF((obs - params), p);

   j = 0;
   for(i = 0; i < m_NumParams; i++)
   {
      pParam = pGroup->GetParamPtr(i);
      if(m_bHoldParam[i] == false)
      {
         est = pParam->GetEstVal();

         stdErr = sqrt(m_pCovar[j][j]);
         m_pCIupr[j] = est + (tStat * stdErr);
         m_pCIlwr[j] = est - (tStat * stdErr);
         j++;
      }
   }/* end for() */

   NEW_PRINT("double *", params);
   coeff = new double *[params];
   MEM_CHECK(coeff);

   for(i = 0; i < params; i++)
   {
      NEW_PRINT("double", params);
      coeff[i] = new double [params];
   }
   MEM_CHECK(coeff[i-1]);

   for(i = 0; i < params; i++){
      for(j = 0; j < params; j++){
         coeff[i][j] = m_pCovar[i][j] / sqrt(m_pCovar[i][i]*m_pCovar[j][j]);    
      }/* end for() */
   }/* end for() */

   p = (double)params;
   v = (obs - params);
   //compute rhs of equation 5.5.6 in Draper & Smith (page 145)
   tmp = pow(exp(GammaLn((0.5 * p) + 1.00)), (2.00)/p);   
   tmp *= (4.00/(MY_PI*p));
   tmp *= FdistInvCDF(1, v, (1.00-alpha));
   /* using rhs result, compute the volume-equivalent ellipsoid percentage
   [i.e. 100*(1-theta), where theta is described in Draper and Smith] */
   m_EllipsePct = 100.00*FdistCDF((int)p, v, 0.00, tmp);

   for(i = 0; i < params; i++){delete [] coeff[i];}
   delete [] coeff;
}/* end ClacCI() */

/******************************************************************************
CalcBealeAndLinssen()

Calculate Beale and Linssen measures of non-linearity. The code is based on
the formulation of Linssen's measure given by Christensen and Cooley in 
"Evaluation of confidence intervals for a steady-state leaky aquifer model",
Advances in Water Resources, Vol. 22, No. 8, page 809, Equations 4-7. 
Beale's measure uses these sames equations but with the lineasrized obs. 
matrix in the denominator replaced with the computed observation matrix. 
Parameter sets are computed as recommended by Cooley and Naff in "Techniques
of Water-Resources Investigations of the USGS", Chapter B4:Regression Modeling
of Groundwater Flow (TWRI 3-B4), Page 174, Equation 5.6-14.
******************************************************************************/
void StatsClass::CalcBealeAndLinssen(void)
{
   double statB; //(Beale)
   double statL; //(Linssen)
   double numerSum;
   double denomSumB; //(Beale)
   double denomSumL; //(Linssen)
   ParameterGroup * pGroup;
   ModelBackup * pMod;
   ModelBackup * pTmp;
   int i, j, k;
   int vbi;
   int numSets; //number of parameter sets to use in computing the measure
   double sign;
   double tmp;
   double fstat;
   double * pParams;
   double * pObs;
   double * pParamSet; //parameter set (p x 1)
   double * pTrueObsSet; //model computed observation sets (n x 1)
   double * pAprxObsSet; //linearized approximation of obs sets (n x 1)
   //temporary matrices
   double * pDelta; //p x 1 vector
   double * pNumer; //n x 1 vector
   double ** pNumerT;//1 x n matrix
   double ** pTemp;//1 x n matrix
   double * pDenomB; //n x 1 vector (Beale)
   double ** pDenomBT;//1 x n matrix (Beale)
   double * pDenomL; //n x 1 vector (Linssen)
   double ** pDenomLT;//1 x n matrix (Linssen)
   int n, p;
   double * pResult;

   n = m_NumObs - m_NumHeldObs;
   p = m_NumParams - m_NumHeldParams;

   pGroup = m_pModel->GetParamGroupPtr();
   numSets = 2 * p;

   /*-------------------------------------
   Allocate space for the matrices
   -------------------------------------*/
   NEW_PRINT("ModelBackup", 1);
   pMod = new ModelBackup(m_pModel);

   NEW_PRINT("ModelBackup", 1);
   pTmp = new ModelBackup(m_pModel);

   NEW_PRINT("double", p);
   pParams = new double[p];

   NEW_PRINT("double", n);
   pObs = new double[n];

   NEW_PRINT("double", p);
   pParamSet = new double[p];

   NEW_PRINT("double", n);
   pTrueObsSet = new double [n];

   NEW_PRINT("double", n);
   pAprxObsSet = new double [n];

   NEW_PRINT("double", n);
   pNumer = new double[n];

   NEW_PRINT("double", n);
   pDenomB = new double[n];

   NEW_PRINT("double", n);
   pDenomL = new double[n];

   NEW_PRINT("double *", 1);
   pTemp = new double *[1];

   NEW_PRINT("double *", 1);
   pNumerT = new double *[1];

   NEW_PRINT("double *", 1);
   pDenomBT = new double *[1];

   NEW_PRINT("double *", 1);
   pDenomLT = new double *[1];

   NEW_PRINT("double", p);
   pDelta = new double[p];

   NEW_PRINT("double", 1);
   pResult = new double[1];
   MEM_CHECK(pResult);

   NEW_PRINT("double", n);
   pNumerT[0] = new double[n];

   NEW_PRINT("double", n);
   pDenomBT[0] = new double[n];

   NEW_PRINT("double", n);
   pDenomLT[0] = new double[n];

   NEW_PRINT("double", n);
   pTemp[0] = new double[n];
   MEM_CHECK(pTemp[0]);

   //95% confidence limit --> alpha = 5 and (1 - alpha) = 0.95
   fstat = FdistInvCDF(p, (n - p), 0.95);

   //compute thresholds, which are based on fstat
   m_NonLinThresh = 1.00/fstat;
   m_EffLinThresh = 0.09/fstat;

   /*
   Compute the statistic
   */
   pMod->Store();
   /* -------------------------------------------------------------------------
   Init. optimum parameter set and observations at this set, taking care to not 
   include insensitive observations and/or parameters.
   ------------------------------------------------------------------------- */
   j = 0;
   for(i = 0; i < m_NumParams; i++){
      if(m_bHoldParam[i] == false){
         pParams[j] = pGroup->GetParamPtr(i)->GetEstVal(); 
         j++; }
   }/* end for() */
   j = 0;
   for(i = 0; i < m_NumObs; i++){ 
      if(m_bHoldObs[i] == false){
         pObs[j] = pMod->GetObs(i, true, true);
         j++;}
   }/* end for() */

   //The statistic is computed over the sum of all sets   
   numerSum = 0.00;
   denomSumB = 0.00;
   denomSumL = 0.00;
   vbi = 0;
   for(j = 0; j < numSets; j++)
   {
      //as per eqn. 5.6-14, each set is +/- a deviation from the optimum
      if((j % 2) == 0) {sign = -1.00;}
      else {sign = +1.00;}

      //vbi indexes the parameter whose extreme values are to be computed
      if(((j % 2) == 0) && (j > 0)){ vbi++;}
  
      //compute parameter set according to eqn. 5.6-14 of TWRI 3-B4
      //printf("Parameters for Sample #%d\n", j);
      for(i = 0; i < p; i++)
      {
         tmp = sign * sqrt(p * fstat);
         tmp /= sqrt(m_pCovar[vbi][vbi]);
         tmp *= m_pCovar[vbi][i];

         pParamSet[i] = pParams[i] + tmp;
            
         pDelta[i] = pParamSet[i] - pParams[i];
      }/* end for() */

      k = 0;
      for(i = 0; i < m_NumParams; i++){
         if(m_bHoldParam[i] == false){
            pGroup->GetParamPtr(i)->SetEstVal(pParamSet[k]); 
            k++; }
      }/* end for() */
      
      //compute true observation set by executing model
      m_pModel->Execute();
      m_StatsCount++;
      pTmp->Store();
      
      //printf("Model computed obs. for Sample #%d\n", j);
      k = 0;
      for(i = 0; i < m_NumObs; i++) 
      {
         if(m_bHoldObs[i] == false){
            pTrueObsSet[k] = pTmp->GetObs(i, true, true);
            k++; }
      }/* end for() */

      //compute approximate observation set by linearization at optimal point
      VectMult(m_pJacob, pDelta, pAprxObsSet, n, p);  

      //printf("Linearized obs. for Sample #%d\n", j);
      for(i = 0; i < n; i++){
         pAprxObsSet[i] += pObs[i];
      }/* end for() */

      //compute numerator and denominators of eqn. 4 of Christensen and Cooley
      for(i = 0; i < n; i++)
      {
         pNumer[i] = pTrueObsSet[i] - pAprxObsSet[i];
         pNumerT[0][i] = pTrueObsSet[i] - pAprxObsSet[i];

         pDenomB[i] = pTrueObsSet[i] - pObs[i];
         pDenomBT[0][i] = pTrueObsSet[i] - pObs[i];

         pDenomL[i] = pAprxObsSet[i] - pObs[i];
         pDenomLT[0][i] = pAprxObsSet[i] - pObs[i];
      }/* end for() */

      /*----------------------------------------------------------
      accumulate sums in numer and denoms
      -----------------------------------------------------------*/
      VectMult(pNumerT, pNumer, pResult, 1, n);
      numerSum += pResult[0];

      VectMult(pDenomBT, pDenomB, pResult, 1, n);
      denomSumB += (pResult[0] * pResult[0]);

      VectMult(pDenomLT, pDenomL, pResult, 1, n);
      denomSumL += (pResult[0] * pResult[0]);
   }/* end for() */

   statB = (p * m_Variance)*(numerSum/denomSumB);
   statL = (p * m_Variance)*(numerSum/denomSumL);

   /*--------------------------------
   Free up matrices
   --------------------------------*/
   delete [] pParams;
   delete [] pObs;
   delete [] pParamSet;
   delete [] pTrueObsSet;
   delete [] pAprxObsSet;
   delete [] pNumer;
   delete [] pDenomB;
   delete [] pDenomL;
   delete [] pNumerT[0];
   delete [] pDenomBT[0];
   delete [] pDenomLT[0];
   delete [] pTemp[0];
   delete [] pNumerT;
   delete [] pDenomBT;
   delete [] pDenomLT;
   delete [] pTemp;
   delete [] pDelta;
   delete [] pResult;

   //restore original model configuration
   pMod->FullRestore();
   m_StatsCount++;

   delete pMod;
   delete pTmp;

   m_BealeStat = statB;
   m_LinssenStat = statL;
} /* end CalcBealeAndLinssen() */

/******************************************************************************
CalcHatAndChangeMatrices()

Compute the Hat matrix, as defined by Yager in "Detecting influential 
observations in non-linear regression modeling of groundwater flow", Water 
Resources Research, Vol. 34, Page 1624, Equation 6.

Also compute the 'change' (C) matrix, described by Yager (in WRR) and Belsley 
in Regression Diagnostics: Identifying Influential Data and Sources of 
Collinearity", J. Wiley & Sons, Page 13, Equation 2.3.
******************************************************************************/
void StatsClass::CalcHatAndChangeMatrices(void)
{
   double ** J; //Jacobian
   double ** JT; //transpose of Jacobian
   double ** pTmp; //intermediate matrix
   int i, j;
   int n, p;

   n = m_NumObs - m_NumHeldObs;
   p = m_NumParams - m_NumHeldParams;

   /*-------------------------------------------------------------
   Allocate and assign the square root matrix
   -------------------------------------------------------------*/
   NEW_PRINT("double *", n);
   pTmp = new double*[n];

   NEW_PRINT("double *", n);
   J  = new double *[n];

   NEW_PRINT("double *", p);
   JT = new double *[p];

   for(i = 0; i < n; i++)
   {
      NEW_PRINT("double", p);
      pTmp[i] = new double[p];
   }
   for(i = 0; i < n; i++)
   {
      NEW_PRINT("double", p);
      J[i] = new double[p];
   }
   for(i = 0; i < p; i++)
   {
      NEW_PRINT("double", n);
      JT[i] = new double[n];
   }
   MEM_CHECK(JT[i-1]);
   
   /*-------------------------------------------------------------
   Assign the J and JT matrices
   -------------------------------------------------------------*/
   for(i = 0; i < p; i++)
   {
      for(j = 0; j < n; j++)
      {
        J[j][i] = m_pJacob[j][i];
        JT[i][j] = J[j][i];
      }
   }/* end for() */

   /*-------------------------------------------------------------
   Assign the intermediate, Hat and Change matrices
   -------------------------------------------------------------*/
   MatMult(J, m_pInvNormal, pTmp, n, p, p);
   MatMult(pTmp, JT, m_pHat, n, p, n);
   MatMult(m_pInvNormal, JT, m_pChange, p, p, n);   
   
   /*-------------------------------------------------------------
   Free up matrices
   -------------------------------------------------------------*/
   for(i = 0; i < n; i++){delete [] J[i];}
   delete [] J;
   for(i = 0; i < p; i++){delete [] JT[i];}
   delete [] JT;
   for(i = 0; i < n; i++){delete [] pTmp[i];}
   delete [] pTmp;
} /* end CalcHatAndChangeMatrices() */

/******************************************************************************
ClacCooksD()

Calculate Cook's D, a measure of the influence of observations on model 
parameters. The code is based on the formulation given by Yager in "Detecting 
influential observations in non-linear regression modeling of groundwater flow", 
Water Resources Research, Vol. 34, Page 1624, Equation 5.
******************************************************************************/
void StatsClass::CalcCooksD(void)
{
   double ei, p, ss, Hii;
   double numer, denom;
   int i, n;

   n = m_NumObs - m_NumHeldObs;
   p = (double)(m_NumParams - m_NumHeldParams);

   //compute the Hat matrix
   CalcHatAndChangeMatrices();

   ss = m_Variance;
   m_CooksAvgLvg = 0.00;
   for(i = 0; i < n; i++)
   {      
      ei = m_pResid[i];
      Hii = m_pHat[i][i];
      numer = ei * ei * Hii;
      denom = p * ss * (1.00 - Hii) * (1.00 - Hii);
      m_pCooksD[i] = (numer / denom);
      m_CooksAvgLvg += Hii;
   }/* end for() */      
   m_CooksAvgLvg /= n;
   m_CooksInfluThresh = 4.00 / n;

   //now count the number of influential observations
   m_NumInfluLvg = 0;
   m_NumInfluCooks = 0;
   for(i = 0; i < n; i++)
   {
      if(fabs(m_pCooksD[i]) > m_CooksInfluThresh){ m_NumInfluCooks++;}
      if(fabs(m_pHat[i][i]) > m_CooksAvgLvg){ m_NumInfluLvg++;}
   }/* end for() */
} /* end CalcCooksD() */

/******************************************************************************
ClacDFBETAS()

Calculate DFBETAS, a measure of the influence of observations on model 
parameters. Computation is taken from the formula given by Yager in "Detecting 
influential observations in on-linear regression modeling of groundwater flow", 
Water Resources Research, Vol. 34, Page 1624, Equations 11-13.
******************************************************************************/
void StatsClass::CalcDFBETAS(void)
{
   double ssi, fi, ss, Cji, sumCjk, Cjk, n, p, Hii, val;
   int i, j, k;   
   bool influFlag; //if true, the obs is influential

   n = (double)(m_NumObs - m_NumHeldObs);
   p = (double)(m_NumParams - m_NumHeldParams);

   //compute the Hat matrix
   CalcHatAndChangeMatrices();

   ss = m_Variance;
   for(i = 0; i < (int)n; i++)
   {
      Hii = m_pHat[i][i];
      fi  = m_pResid[i];
      ssi = (1.00/(n - p - 1.00));
      ssi *= (((n - p)* ss) - ((fi*fi)/(1.00 - Hii)));

      for(j = 0; j < (int)p; j++)
      {
         Cji = m_pChange[j][i];

         sumCjk = 0.00;
         for(k = 0; k < (int)n; k++)
         {
            Cjk = m_pChange[j][k];
            sumCjk += (Cjk * Cjk);
         } /* end for() */

         val = (Cji / sqrt(sumCjk)) * (fi / (sqrt(ssi) * (1.00 - Hii)));
         m_pDFBETAS[i][j] = val;
      }/* end for() */
   }/* end for() */
   //Suggested threhold for influential DFBETAS is 2/sqrt(n)
   m_DfbetaInfluThresh = 2.00/sqrt(n);

   //count how many are influential
   m_NumInfluDfbeta = 0;
   for(i = 0; i < (int)n; i++)
   {
      influFlag = false;
      for(j = 0; j < (int)p; j++)
      {
         //set the flag if any one of the parameters are influenced....
         if(fabs(m_pDFBETAS[i][j]) > m_DfbetaInfluThresh){influFlag = true;}
      }/* end for() */
      if(influFlag == true){ m_NumInfluDfbeta++;}
   } /* end for() */
}/* end ClacDFBETAS() */

/******************************************************************************
CalcNormProbPlot()

Calculates data used in the normal probability plot. This data is a set of 
(x,y) coordinates. The x-coord is the expected value of the std. normal order 
statistics and is stored in m_pExpResid array. The y-coord is the ordered set
of weighted residuals and is stored in the m_pOrdResid array. Finally the 
correlation coefficient between the ordered weighted residuals and the order 
statisics is computated and stored in m_OrdCorrCoef.
******************************************************************************/
void StatsClass::CalcNormProbPlot(void)
{
   int i, n;
   double * pTau;
   double avgRes; //average weighted residual
   double numer, denom1, denom2, ei, m, pi;

   n = m_NumObs - m_NumHeldObs;

   /*----------------------------------------------
   Compute the ordered weighted residuals, and the 
   average weighted residual.
   ----------------------------------------------*/
   avgRes = 0.00;
   for(i = 0; i < n; i++)
   {      
      m_pOrdResid[i] = m_pResid[i];
      avgRes += m_pOrdResid[i];
   }/* end for() */
   avgRes /= n;

   SortInc(m_pOrdResid, n);

   /*----------------------------------------------
   Compute the expected values of std. norm. order 
   stats. using the Snedecor & Cochran approximation
   described by David W. Sabo (BCIT) in "Normal 
   Probability Plots", page #3. 
   ----------------------------------------------*/
   for(i = 0; i < n; i++)
   {
      pi = (double(i+1) - 0.375)/(n + 0.25);
      m_pExpResid[i] = StdNormInvCDF(pi);
   }/* end for() */   
   
#if 0
   /*----------------------------------------------
   Compute the 'tau' vector, as described by Mary C.
   Hill in "Methods and Guidelines for Effective 
   Model Calibration", page 23.
   ----------------------------------------------*/
   NEW_PRINT("double", n);
   pTau = new double[n];
   MEM_CHECK(pTau);

   for(i = 0; i < n; i++)
   {
      ui = (((double)(i+1) - 0.50)/n);
      pTau[i] = StdNormInvCDF(ui);
   }/* end for() */   
#endif
   pTau=m_pExpResid;

   /*----------------------------------------------
   Compute the numerator of Equation 5 of "Methods 
   and Guidelines for Effective Model Calibration",
   page 23.
   ----------------------------------------------*/
   numer = 0.00;
   for(i = 0; i < n; i++)
   {
      ei = m_pOrdResid[i];
      m = avgRes;      
      numer += ((ei - m)*pTau[i]);
   }/* end for() */   
   numer *= numer;

   /*----------------------------------------------
   Calc. the denominator of Equation 5 of "Methods 
   and Guidelines for Effective Model Calibration",
   page 23.
   ----------------------------------------------*/
   denom1 = 0.00;
   denom2 = 0.00;
   for(i = 0; i < n; i++)
   {
      ei = m_pOrdResid[i];
      m = avgRes;      
      denom1 += ((ei - m)*(ei - m));
      denom2 += (pTau[i]*pTau[i]);
   }/* end for() */   
   
   m_OrdCorrCoeff = (numer / (denom1 * denom2));

   //delete [] pTau;
}/* end CalcNormProbPlot() */

/******************************************************************************
CalcWeightedRy()

Calculate the Ry value for weighted observations vs. weighted measurements, 
this requires calculation of the variance due to the regression 
(sum of squared weighted differences between predicted observations and the 
average observation value).
******************************************************************************/
void StatsClass::CalcWeightedRy(void)
{
   int i, j;
   double avgObs; //average observation
   double avgEst; //average estimate
   double Yhat, Ytrue;
   double numer, denom, tmp1, tmp2;
   ObservationGroup * pObsGroup;
   Observation * pObs;
      
   pObsGroup = m_pModel->GetObsGroupPtr();

   /*----------------------------------------------
   Compute the weighted averages.
   ----------------------------------------------*/
   avgObs = 0.00;
   avgEst = 0.00;
   j = 0;
   for(i = 0; i < m_NumObs; i++)
   {
      if(m_bHoldObs[i] == false){
         pObs = pObsGroup->GetObsPtr(i);
         avgObs += pObs->GetMeasuredVal(true, true);
         avgEst += pObs->GetComputedVal(true, true);
         j++; }
   }/* end for() */

   avgObs /= (double)j;
   avgEst /= (double)j;

   //Calculation of Ry using equation 18, page 21 of Methods and Guidelines...
   //equalivalent to Rxy of Draper and Smith...
   numer = 0.00;
   tmp1 = 0.00;
   tmp2 = 0.00;
   j = 0;
   for(i = 0; i < m_NumObs; i++)
   {  
      if(m_bHoldObs[i] == false){    
         pObs = pObsGroup->GetObsPtr(i);
         Yhat = pObs->GetComputedVal(true, true);
         Ytrue = pObs->GetMeasuredVal(true, true);
         numer += ((Ytrue-avgObs)*(Yhat-avgEst));
         tmp1 += ((Ytrue-avgObs)*(Ytrue-avgObs));
         tmp2 += ((Yhat-avgEst)*(Yhat-avgEst));
         j++; }
   }/* end for() */
   denom = sqrt((tmp1*tmp2));
   m_WeightedRy = numer/denom;
}/* end CalcRy() */

/******************************************************************************
CalcRawRy()

Calculate the Raw Ry value for observations vs. measurements (i.e. no weighting) 
This requires calculation of the variance due to the regression 
(sum of squared differences between predicted observations and the average 
observation value).
******************************************************************************/
void StatsClass::CalcRawRy(void)
{
   int i, j;
   double avgObs; //average observation
   double avgEst; //average estimate
   double Yhat, Ytrue;
   double numer, denom, tmp1, tmp2;
   ObservationGroup * pObsGroup;
   Observation * pObs;
      
   pObsGroup = m_pModel->GetObsGroupPtr();

   /*----------------------------------------------
   Compute the averages.
   ----------------------------------------------*/
   avgObs = 0.00;
   avgEst = 0.00;
   j = 0;
   for(i = 0; i < m_NumObs; i++)
   {
      if(m_bHoldObs[i] == false){
         pObs = pObsGroup->GetObsPtr(i);
         avgObs += (pObs->GetMeasuredVal(false, false));
         avgEst += (pObs->GetComputedVal(false, false));
         j++; }
   }/* end for() */

   avgObs /= (double)j;
   avgEst /= (double)j;

   //Calculation of Ry using equation 18, page 21 of Methods and Guidelines...
   //equalivalent to Rxy of Draper and Smith...
   numer = 0.00;
   tmp1 = 0.00;
   tmp2 = 0.00;
   j = 0;
   for(i = 0; i < m_NumObs; i++)
   {  
      if(m_bHoldObs[i] == false){    
         pObs = pObsGroup->GetObsPtr(i);
         Yhat = (pObs->GetComputedVal(false, false));
         Ytrue =(pObs->GetMeasuredVal(false, false));
         numer += ((Ytrue-avgObs)*(Yhat-avgEst));
         tmp1 += ((Ytrue-avgObs)*(Ytrue-avgObs));
         tmp2 += ((Yhat-avgEst)*(Yhat-avgEst));
         j++; }
   }/* end for() */
   denom = sqrt((tmp1*tmp2));
   m_RawRy = numer/denom;
}/* end CalcRawRy() */

/******************************************************************************
CalcMMRI()

Compute alternative measures of model fit. These MMRI measures are designed
to help with model selection.

Takes 1 arg, a boolean that indicates whether or not the normal matrix is
invertible.
******************************************************************************/
void StatsClass::CalcMMRI(bool bInv)
{
   int n = m_NumObs - m_NumHeldObs;
   int k = m_NumParams - m_NumHeldParams + 1;
   double ss = m_Phi/n;

   m_MMRI.AIC  = (n * log(ss)) + (2.00 * k);
   if((n - k - 1) > 0)
   {
      m_bDOF = true;
      m_MMRI.AICc = (n * log(ss)) + (2.00 * k) + ((2.00 * k * (k + 1))/(n - k - 1));
      m_MMRI.AICu = (n * log((n*ss)/(n-k))) + (2.00 * k) + ((2.00 * k * (k + 1))/(n - k - 1));
   }
   else
   {
      m_bDOF = false;
   }
   m_MMRI.BIC  = (n * log(ss)) + (k * log((double)n));
   m_MMRI.HQ   = (n * log(ss)) + (2.00 * k * log(log((double)n)));
}/* enc CalcMMRI() */

/******************************************************************************
CalcSensitivities()

Compute parameter sensitivities. Sensitivity formulas are taken from:
"Methods and Guidelines for Effective Model Calibration" by Mary Hill, USGS
1998, page 14-16. 
******************************************************************************/
void StatsClass::CalcSensitivities(void)
{
   ParameterGroup * pGroup;
   ParameterABC * pParam;
   int i,j,k,n,p;
   int jj;
   double diff; //partial derivative from the Jacobian
   double bj;   //parameter value
   double wt;
   double ss;

   n = m_NumObs - m_NumHeldObs;
   p = m_NumParams - m_NumHeldParams;
   
   /*----------------------------------------------------
   Compute scaled sensitivities using full weight matrix 
   (Equation 9 in Hill, Page 15)
   ----------------------------------------------------*/
   pGroup = m_pModel->GetParamGroupPtr();
   for(i = 0; i < n; i++)
   {
      jj = 0;
      for(j = 0; j < m_NumParams; j++)
      {
         if(m_bHoldParam[j] == false){
            m_pScaledSens[i][jj] = 0.00;

            pParam = pGroup->GetParamPtr(j);
            bj = pParam->GetEstVal();

            for(k = 0; k < n; k++)
            {                  
               diff = m_pJacob[k][jj];
               if(i == k)
                 wt = 1.00;
               else
                 wt = 0.00;

               m_pScaledSens[i][jj] += (diff * bj * wt);
            }
            jj++;
         }/* end if() */
      }/* end for() */
   }/* end for() */
   
   /*----------------------------------------------------
   Compute composite scaled sensitivities.
   (Equation 10 in Hill, Page 15)
   ----------------------------------------------------*/
   for(j = 0; j < p; j++)
   {
      m_pCompScaledSens[j] = 0.00;

      for(i = 0; i < n; i++)
      {
         ss = m_pScaledSens[i][j] * m_pScaledSens[i][j];
         m_pCompScaledSens[j] += ss;         
      }/* end for() */
       
      m_pCompScaledSens[j] /= n;
      m_pCompScaledSens[j] = sqrt(m_pCompScaledSens[j]);
   }/* end for() */

   /*----------------------------------------------------
   Compute one-percent scaled sensitivities.
   (Equation 11 in Hill, Page 16)
   ----------------------------------------------------*/
   for(i = 0; i < n; i++)
   {
      jj = 0;
      for(j = 0; j < m_NumParams; j++)
      {
         if(m_bHoldParam[j] == false){
            m_pPctScaledSens[i][jj] = 0.00;

            pParam = pGroup->GetParamPtr(j);
            bj = pParam->GetEstVal();
            diff = m_pJacobUW[i][jj];

            m_pPctScaledSens[i][jj] += (diff * bj / 100.00);
            jj++;
         }/* end if() */
      }/* end for() */
   }/* end for() */
}/* end CalcSensitivities() */

/******************************************************************************
CalcPredictions()

Calculate linear confidence intervals of user-specified predicions, using 
students t-distribution (see Hill, "Methods and Guidelines for Effective Model
Calibration", page 29-31 for details.
******************************************************************************/
void StatsClass::CalcPredictions(bool bStats, double ** pV, int np)
{
   int i, j, k, nobs;
   double sum, sd, alpha, p, tStat, est;
   double p1, p2, V;
   int nrv = m_pPredictions->GetNumRespVars();
   
   if(nrv <= 0) return;

   nobs = m_NumObs;
   alpha = 1.00 - (m_CIpct/100.00);
   p = 1.00 - (alpha / 2.00);
   tStat = StudentInvCDF((nobs - np), p);

   for(i = 0; i < nrv; i++)
   {
      m_Pred[i] = est = m_pPredictions->GetRespVarPtr(i)->GetCurrentVal();

      if(bStats)
      {
         sum = 0.00;
         for(j = 0; j < np; j++)
         {
            for(k = 0; k < np; k++)
            {
               p1 = m_pJacPred[i][j];
               p2 = m_pJacPred[i][k];
               V  = pV[j][k];
               sum += p1*p2*V;
            }
         }
         m_PredSD[i] = sd = sqrt(sum);      
         m_PredCIupr[i] = est + (tStat * sd);
         m_PredCIlwr[i] = est - (tStat * sd);
      }/* end if(stats) */
   }/* for each prediction */

   if(!bStats)
   {
      delete [] m_PredSD; m_PredSD = NULL;
      delete [] m_PredCIupr; m_PredCIupr = NULL;
      delete [] m_PredCIlwr; m_PredCIlwr = NULL;
   }
}/* end ClacPredictions() */

/******************************************************************************
WriteMetrics()

Reports on the setup of the Math Class and also various run-time metrics.
******************************************************************************/
void StatsClass::WriteMetrics(FILE * pFile)
{
   fprintf(pFile, "\nFinite Difference Metrics\n");
   fprintf(pFile, "Difference Type    : ");

   if     (m_DiffType == FD_FORWARD){ fprintf(pFile, "Forward\n");}   
   else if(m_DiffType == FD_OUT_CEN){ fprintf(pFile, "Outside Central\n");}
   else if(m_DiffType == FD_PAR_CEN){ fprintf(pFile, "Parabolic Central\n");}
   else if(m_DiffType == FD_FIT_CEN){ fprintf(pFile, "Best-fit Central\n");}
   else                             { fprintf(pFile, "Unknown\n");}

   fprintf(pFile, "Increment Type    : ");

   if     (m_DiffIncType == FD_RANGE_REL){ fprintf(pFile, "Range-Relative\n");}   
   else if(m_DiffIncType == FD_VALUE_REL){ fprintf(pFile, "Value-Relative\n");}
   else if(m_DiffIncType == FD_ABSOLUTE) { fprintf(pFile, "Absolute\n");}
   else if(m_DiffIncType == FD_OPTIMAL)  { fprintf(pFile, "Optimal\n");}
   else                                  { fprintf(pFile, "Range-Relative\n");}
 
   fprintf(pFile, "Finite Difference Increments\n");
   for(int i = 0; i < m_NumParams; i++)
   {
      fprintf(pFile, "%-12s : ", GetParameterName(i));
      if(m_DiffIncType != FD_OPTIMAL){
         fprintf(pFile, "%lf\n", m_pDiffInc[i]);}
      else{
         fprintf(pFile, "optimal\n");}
   }
   fprintf(pFile, "Finite Difference Mimumum Increment : %E\n", m_MinInc);
   fprintf(pFile, "Jacobian Evals     : %d\n", m_DiffCount);
   fprintf(pFile, "Optimal Step Evals : %d\n", m_StepCount);
   fprintf(pFile, "Statistics Evals   : %d\n", m_StatsCount);
}/* end WriteMetrics() */

/******************************************************************************
WriteStats()

Write statisitics to the given file.
******************************************************************************/
void StatsClass::WriteStats(FILE * pFile)
{
   if(m_bNoStats) return;

   ObservationGroup * pObsGroup;
   Observation * pObs;
   ParameterGroup * pGroup;
   ParameterABC * pParam;
   int i, j, ii, jj, n, np;
   double coeff, low, hi;

   //temporary storage for obs. residuals
   double m, p, w, d; //measured, predicted, weight, difference

   n = m_NumObs - m_NumHeldObs;
   np = m_NumParams - m_NumHeldParams;

   pGroup = m_pModel->GetParamGroupPtr();
   pObsGroup = m_pModel->GetObsGroupPtr();

   fprintf(pFile, "\nStatistical Output\n");
   fprintf(pFile, "\nUntransformed WSSE : %E\n", ((WSSE *)(m_pModel->GetObjFuncPtr()))->CalcUntransformedObjFunc());
   if((m_NumHeldObs > 0) || (m_NumHeldParams > 0)){
      fprintf(pFile, "********************** NOTE **********************\n");
      if((m_bOkToHoldObs == true) && (m_bOkToHoldParams == true))
      {
         fprintf(pFile, "Insensitive observations (%d) and/or parameters (%d)\n", 
                        m_NumHeldObs, m_NumHeldParams);
         fprintf(pFile, "were detected and have not been included in the \n");
         fprintf(pFile, "following statistical calculations.\n");
      }
      else if(m_bOkToHoldObs == true) //hold obsrvations only
      {
         fprintf(pFile, "Insensitive observations (%d) were detected and have \n",
                        m_NumHeldObs);
         fprintf(pFile, "not been included in the following statistical \n");
         fprintf(pFile, "calculations.\n");
      }
      else //hold parameters only
      {
         fprintf(pFile, "Insensitive parameters (%d) were detected and have \n",
                        m_NumHeldParams);
         fprintf(pFile, "not been included in the following statistical \n");
         fprintf(pFile, "calculations.\n");
      }
      fprintf(pFile, "**************************************************\n");
      fprintf(pFile, "\nAdjusted Obj.Func. : %E\n", m_Phi);
      if(m_bOkToHoldParams == true)
      {
         fprintf(pFile, "\nParameter      Value            Sensitive?\n");
         for(i = 0; i < m_NumParams; i++)
         {
            fprintf(pFile, "%-14s ", pGroup->GetParamPtr(i)->GetName());
            pGroup->GetParamPtr(i)->Write(pFile, WRITE_SCI);
            if(m_bHoldParam[i] == false) 
               fprintf(pFile, "  YES\n");
            else
               fprintf(pFile, "  NO\n");
         }/* end for() */
      }/* end if() */
   }/* end if() */

   fprintf(pFile, "\nObservation Residuals\n");
   fprintf(pFile, "Observation    Measured       Simulated      Weight          Residual(Transformed and Weighted)");
   if(m_bOkToHoldObs == true){
      fprintf(pFile, "   Sensitive?");
   }
   fprintf(pFile, "\n");

   for(i = 0; i < m_NumObs; i++)
   {
      pObs = pObsGroup->GetObsPtr(i);
      m = pObs->GetMeasuredVal(false, false);
      p = pObs->GetComputedVal(false, false);
      w = GetObsWeight(pObs);
      d = pObs->CalcResidual(true, true);
      fprintf(pFile, "%-12s  %E  %E  %E  %+E   ", pObs->GetName(), m, p, w, d);
      if(m_bOkToHoldObs == true){
         if(m_bHoldObs[i] == false){ fprintf(pFile, "   YES");}
         else                      { fprintf(pFile, "   NO");}
      }
      fprintf(pFile, "\n");
   }/* end for() */   
   fprintf(pFile, "\nCorrelation between raw measured and simulated observations (no transformation or weighting)\n");
   fprintf(pFile, "Ry         : %6.3lf\n", m_RawRy);
   fprintf(pFile, "Ry-squared : %6.3lf\n", m_RawRy*m_RawRy);

   fprintf(pFile, "\nCorrelation between measured and simulated observations (with transformation and weighting)\n");
   fprintf(pFile, "Rw         : %6.3lf\n", m_WeightedRy);
   fprintf(pFile, "Rw-squared : %6.3lf\n", m_WeightedRy*m_WeightedRy);

   if(m_RunsTestFlag == true)
   {
      fprintf(pFile, "\nRuns Test on Residuals\n");
      fprintf(pFile, "NOTE: Residuals of zero are counted as positive.\n");
      if(m_Runs.bSuccess == true)
      {
         fprintf(pFile, "Positive Residuals : %d\n", m_Runs.pos);
         fprintf(pFile, "Negative Residuals : %d\n", m_Runs.neg);
         fprintf(pFile, "Number of Runs     : %d\n", m_Runs.runs);
         fprintf(pFile, "Lower-tail critical value (alpha=0.1) : %d\n", m_Runs.clwr);
         fprintf(pFile, "Upper-tail critical value (alpha=0.1) : %d\n", m_Runs.cupr);
         if((m_Runs.runs < m_Runs.clwr) || (m_Runs.runs > m_Runs.cupr)) fprintf(pFile, "Runs appear to be clustered (i.e. non-random)\n");
         else fprintf(pFile, "Runs appear to be randomly distributed\n");
      }
      else fprintf(pFile, "The Runs Test was unsuccessful\n");
   }

   if(m_AutorunFunctionFlag == true)
   {
      fprintf(pFile, "\nAutorun Function Test for Lag-1 Autocorrelation of Residuals\n");
      fprintf(pFile, "Lag-1 Autorun Function (r1)           : %lf\n", m_AR.r1);
      fprintf(pFile, "Variance of Lag-1 Autorun Function    : %lf\n", m_AR.var);
      fprintf(pFile, "Approximate Lag-1 Variance            : %lf\n", m_AR.vpx);
      fprintf(pFile, "Std. Dev. of Lag-1 Autorun Function   : %lf\n", sqrt(m_AR.var));
      fprintf(pFile, "Median Residual (m)                   : %lf\n", m_AR.med);
      fprintf(pFile, "Number of Surpluses (ei > m)          : %d\n", m_AR.sur);
      fprintf(pFile, "Number of Deficits  (ei <= m)         : %d\n", m_AR.def);
      fprintf(pFile, "Number of Lag-1 Surplus Pairs (n1)    : %d\n", m_AR.n1);         
      fprintf(pFile, "Lower-tail critical value (alpha=0.1) : %lf\n", m_AR.clwr);
      fprintf(pFile, "Upper-tail critical value (alpha=0.1) : %lf\n", m_AR.cupr);
      if((m_AR.r1 < m_AR.clwr) || (m_AR.r1 > m_AR.cupr)) fprintf(pFile, "Lag-1 residuals appear to be correlated (i.e. persistent)\n");
      else fprintf(pFile, "Lag-1 residuals do NOT appear to be correlated\n");
   }
   
   if(m_StdDevFlag == true)
   {
      fprintf(pFile, "\nError Variance and Standard Error of the Regression\n");
      fprintf(pFile, "S^2 : %E\n", m_Variance);
      fprintf(pFile, "S   : %E\n", sqrt(m_Variance));
   }/* end if() */

   if(m_MMRI_Flag == true)
   {
      fprintf(pFile, "\nMMRI (Alternative Measures of Model Fit)\n");
      fprintf(pFile, "Akaike Information Criterion           (AIC)  : %02lf\n", m_MMRI.AIC);
      if(m_bDOF == true)
      {
         fprintf(pFile, "Corrected Akaike Information Criterion (AICc) : %02lf\n", m_MMRI.AICc);
         fprintf(pFile, "Corrected Unbiased Akaike Criterion    (AICu) : %02lf\n", m_MMRI.AICu);
      }
      else
      {
         fprintf(pFile, "Corrected Akaike Information Criterion (AICc) : not computed\n");
         fprintf(pFile, "Corrected Unbiased Akaike Criterion    (AICu) : not computed\n");
      }
     
      fprintf(pFile, "Bayesian Information Criterion         (BIC)  : %02lf\n", m_MMRI.BIC);
      fprintf(pFile, "Hannan and Quinn's Criterion           (HQ)   : %02lf\n", m_MMRI.HQ);
   }

   if(m_StdErrFlag == true)
   {   
      fprintf(pFile, "\nParameter Variance-Covariance\n");      
      fprintf(pFile, "               ");
      for(i = 0; i < m_NumParams; i++)
      {
         pParam = pGroup->GetParamPtr(i);
         pParam->Write(pFile, WRITE_BNR);
         fprintf(pFile, " ");
      }
      fprintf(pFile, "\n");

      ii = 0;
      for(i = 0; i < m_NumParams; i++)
      {
         pParam = pGroup->GetParamPtr(i);
         pParam->Write(pFile, WRITE_BNR);

         jj = 0;
         for(j = 0; j < m_NumParams; j++)
         {
            if((m_bHoldParam[j] == false) && (m_bHoldParam[i] == false))
            {
               fprintf(pFile, "%+E  ", m_pCovar[ii][jj]);
            }
            else
            {
               fprintf(pFile, "not_computed    ");
            }

            if(m_bHoldParam[j] == false){ jj++;}
         }
         if(m_bHoldParam[i] == false){ ii++;}
         fprintf(pFile, "\n");
      }/* end for() */   

      fprintf(pFile, "\nParameter Standard Error\n");
      ii = 0;
      for(i = 0; i < m_NumParams; i++)
      {
         pParam = pGroup->GetParamPtr(i);
         pParam->Write(pFile, WRITE_BNR);
         if(m_bHoldParam[i] == false){
            fprintf(pFile, " : %E \n", sqrt(m_pCovar[ii][ii]));
            ii++;}
         else {
            fprintf(pFile, " : not_computed\n"); }
      }/* end for() */
   }/* end if() */

   if(m_CorrCoefFlag == true)
   {
      fprintf(pFile, "\nParameter Correlation\n");
      fprintf(pFile, "               ");
      pGroup->Write(pFile, WRITE_BNR);
      fprintf(pFile, "\n");

      ii = 0;
      for(i = 0; i < m_NumParams; i++)
      {
         pParam = pGroup->GetParamPtr(i);
         pParam->Write(pFile, WRITE_BNR);

         jj = 0;
         for(j = 0; j < m_NumParams; j++)
         {
            if((m_bHoldParam[i] == false) && (m_bHoldParam[j] == false)) {
               coeff = m_pCovar[ii][jj] / sqrt(m_pCovar[ii][ii]*m_pCovar[jj][jj]);    
               fprintf(pFile, "%+6.3lf         ", coeff); }
            else {
               fprintf(pFile, "n/a            "); }

            if(m_bHoldParam[j] == false){ jj++;}
         }
         if(m_bHoldParam[i] == false){ ii++;}
         fprintf(pFile, "\n");
      }/* end for() */   
   }/* end if() */

   if(m_CIflag == true)
   {
      fprintf(pFile, "\nLinear Confidence Intervals (%.2lf%%)\n", m_CIpct);
      fprintf(pFile, "Parameter      Lower Limit     Upper Limit\n");
      ii = 0;
      for(i = 0; i < m_NumParams; i++)
      {
         pParam = pGroup->GetParamPtr(i);
         pParam->Write(pFile, WRITE_TX_BNR);
         if(m_bHoldParam[i] == false){
            low = pParam->ConvertOutVal(m_pCIlwr[ii]);
            hi = pParam->ConvertOutVal(m_pCIupr[ii]);
            fprintf(pFile, "%+E  %+E\n", low, hi);
            ii++;}
         else{
            fprintf(pFile, "not_computed    not_computed\n");}
      }/* end for() */
   }/* end if() */

   if((m_BealeFlag == true) || (m_LinssenFlag == true))
   {
      fprintf(pFile, "\nNon-Linearity Measures\n");
   }/* end if() */

   if(m_BealeFlag == true)
   {
      fprintf(pFile, "Beale (N)     : %E\n", m_BealeStat);
      if(m_BealeStat > m_NonLinThresh)
      {
         fprintf(pFile, "Assessment    : Non-Linear\n");
      }
      else if(m_BealeStat < m_EffLinThresh)
      {
         fprintf(pFile, "Assessment    : Linear\n");
      }
      else
      {
         fprintf(pFile, "Assessment    : Uncertain\n");
      }
   }/* end if() */
   if(m_LinssenFlag == true)
   {
      if(m_BealeFlag == true){fprintf(pFile, "\n");}

      fprintf(pFile, "Linssen (M^2) : %E\n", m_LinssenStat);
      if(m_LinssenStat > m_NonLinThresh)
      {
         fprintf(pFile, "Assessment    : Non-Linear\n");
      }
      else if(m_LinssenStat < m_EffLinThresh)
      {
         fprintf(pFile, "Assessment    : Linear\n");
      }
      else
      {
         fprintf(pFile, "Assessment    : Uncertain\n");
      }
   }/* end if() */
   if((m_BealeFlag == true) || (m_LinssenFlag == true))
   {
      fprintf(pFile, "\nThresholds for N and/or M^2\n");
      fprintf(pFile, "Non-linear : > %E\n", m_NonLinThresh);
      fprintf(pFile, "Linear     : < %E\n",m_EffLinThresh);
   }/* end if() */

   if(m_NormPlotFlag == true)
   {
      fprintf(pFile, "\nNormalized Residuals\n");
      fprintf(pFile, "r_expected      r_ordered\n");
      
      for(i = 0; i < n; i++)
      {      
         fprintf(pFile, "%+E  %+E\n", m_pExpResid[i], m_pOrdResid[i]);
      }/* end for() */      

      double rcrit=GetCritValNormPPCC(n);

      fprintf(pFile, "\nNormal probability correlation coefficient\n");
      fprintf(pFile, "R2N                  : %6.4lf\n", m_OrdCorrCoeff);
      fprintf(pFile, "RN                   : %6.4lf\n", sqrt(m_OrdCorrCoeff));
      fprintf(pFile, "RN Critical Value    : %6.4lf\n", rcrit);
      fprintf(pFile, "Normality Assessment : ");
      if(rcrit > sqrt(m_OrdCorrCoeff))
         fprintf(pFile, "Residuals do NOT appear to be normally distributed\n");
      else
         fprintf(pFile, "Residuals appear to be normally distributed\n");

      double mean = CalcMean(m_pOrdResid, n);
      double median = CalcMedian(m_pOrdResid, n);
      double sd = CalcStdDev(m_pOrdResid, n, CENTRAL_TEND_MEAN);
      double skewness = CalcSkewness(m_pOrdResid, n);
      double kurtosis = CalcKurtosis(m_pOrdResid, n);
      fprintf(pFile, "\nSample Statistics for Residuals\n");
      fprintf(pFile, "Minumum       : %E\n", m_pOrdResid[0]);
      fprintf(pFile, "Maximum       : %E\n", m_pOrdResid[n-1]);
      fprintf(pFile, "Mean          : %E\n", mean);
      fprintf(pFile, "Median        : %E\n", median);
      fprintf(pFile, "Std Deviation : %E\n", sd);
      fprintf(pFile, "Skewness      : %6.3lf\n", skewness);
      fprintf(pFile, "Kurtosis      : %6.3lf\n", kurtosis);
      fprintf(pFile, "(Skewness and Kurtosis should be close to 0 if residuals are normally distributed)\n");

   }/* end if() */

   if(m_BestBoxCoxFlag == true)
   {
      fprintf(pFile, "\nEstimated Optimal Box-Cox Transformation\n");
     fprintf(pFile, "Lambda : %lf\n", m_BestBoxCoxVal);
   }/* end if() */

   if((m_CooksFlag == true) || (m_DfbetasFlag == true))
   {
      fprintf(pFile, "\nMeasures of Observation Influence\n");
   }/* end if() */

   if(m_CooksFlag == true)
   {
      fprintf(pFile, "\nCook's D\n");
      fprintf(pFile, "Observation    Leverage   infl.?  Di         infl.?\n");
      ii = 0;
      for(i = 0; i < m_NumObs; i++)   
      {
         pObs = pObsGroup->GetObsPtr(i);
         if(m_bHoldObs[i] == false) {
            fprintf(pFile, "%-12s  %.2E  ", pObs->GetName(), m_pHat[ii][ii]);
            if(fabs(m_pHat[ii][ii]) > m_CooksAvgLvg){ fprintf(pFile, "yes     ");}
            else{ fprintf(pFile, "no      ");}
            fprintf(pFile, "%.2E  ", m_pCooksD[ii]);
            if(fabs(m_pCooksD[ii]) > m_CooksInfluThresh){ fprintf(pFile, "yes\n");}
            else{ fprintf(pFile, "no\n");} 
            ii++;}
      }/* end for() */
      
      fprintf(pFile, "\nNumber of  influential Leverage : %d\n", m_NumInfluLvg);
      fprintf(pFile, "Number of influential Di        : %d\n", m_NumInfluCooks);

      fprintf(pFile, "\nThresholds for Cook's D\n");
      fprintf(pFile, "Di       > %.2E\n", m_CooksInfluThresh);
      fprintf(pFile, "Leverage > %.2E\n", m_CooksAvgLvg);
   }/* end if() */

   if(m_DfbetasFlag == true)
   {
      fprintf(pFile, "\nDFBETAS\n");
      fprintf(pFile, "Observation    ");
      for(i = 0; i < m_NumParams; i++)
      {
         if(m_bHoldParam[i] == false){
            pParam = pGroup->GetParamPtr(i);
            pParam->Write(pFile, WRITE_BNR);
            fprintf(pFile, "infl.?  ");}
      }/* end for() */
      fprintf(pFile, "\n");

      ii = 0;
      for(i = 0; i < m_NumObs; i++)
      {
         pObs = pObsGroup->GetObsPtr(i);

         if(m_bHoldObs[i] == false) {
            fprintf(pFile, "%-12s  ", pObs->GetName());

            jj = 0;
            for(j = 0; j < m_NumParams; j++)
            {
               if(m_bHoldParam[j] == false){
                  fprintf(pFile, "%+.2E     ", m_pDFBETAS[ii][jj]);
                  if(fabs(m_pDFBETAS[ii][jj]) > m_DfbetaInfluThresh){fprintf(pFile, "yes     ");}
                  else{fprintf(pFile, "no      ");}
                  jj++;}
            }/* end for() */
            fprintf(pFile, "\n");
            ii++;}
      }/* end for() */      

      fprintf(pFile, "\nNumber of influential DFBETAS : %d\n", m_NumInfluDfbeta);

      fprintf(pFile, "\nThreshold for DFBETAS\n");
      fprintf(pFile, "|DFBETASij| > %.2E\n", m_DfbetaInfluThresh);      
   }/* end if() */

   if(m_SensFlag == true)
   {
      fprintf(pFile, "\nParameter Sensitivities\n");
      fprintf(pFile, "\nDimensionless Scaled Sensitivities\n");
      fprintf(pFile, "Observation    ");
      pGroup->Write(pFile, WRITE_BNR);
      fprintf(pFile, "\n");
      ii = 0;
      for(i = 0; i < m_NumObs; i++)
      {
         if(m_bHoldObs[i] == false){
            pObs = pObsGroup->GetObsPtr(i);
            fprintf(pFile, "%-12s  ", pObs->GetName());
            jj = 0;
            for(j = 0; j < m_NumParams; j++)
            {
               if(m_bHoldParam[j] == false){
                  fprintf(pFile, "%+.5E  ", m_pScaledSens[ii][jj]);
                  jj++;}
               else{
                  fprintf(pFile, "not_computed   ");}
            }/* end for() */
            fprintf(pFile, "\n");
            ii++;}
      }/* end for() */

      fprintf(pFile, "\n1-Percent Scaled Sensitivities\n");
      fprintf(pFile, "Observation    ");
      pGroup->Write(pFile, WRITE_BNR);
      fprintf(pFile, "\n");
      ii = 0;
      for(i = 0; i < m_NumObs; i++)
      {
         if(m_bHoldObs[i] == false){         
            pObs = pObsGroup->GetObsPtr(i);
            fprintf(pFile, "%-12s  ", pObs->GetName());
            jj = 0;
            for(j = 0; j < m_NumParams; j++)
            {
               if(m_bHoldParam[j] == false){
                  fprintf(pFile, "%+.5E  ", m_pPctScaledSens[ii][jj]);
                  jj++;}
               else{
                  fprintf(pFile, "not_computed   ");}
            }/* end for() */
            fprintf(pFile, "\n");
            ii++;}
      }/* end for() */

      fprintf(pFile, "\nComposite Scaled Sensitivities\n");
      ii = 0;
      for(i = 0; i < m_NumParams; i++)
      {
         pParam = pGroup->GetParamPtr(i);
         pParam->Write(pFile, WRITE_BNR);

         if(m_bHoldParam[i] == false){
            fprintf(pFile, " : %E\n", m_pCompScaledSens[ii]);
            ii++;}
         else{
            fprintf(pFile, " : not_computed\n");}
      }/* end for() */
   } /* end if() */

   if(m_MatricesFlag == true)
   {
      fprintf(pFile, "\nMatrices\n");
      fprintf(pFile, "\nJacobian Matrix (note: includes Transformation, if applicable)\n");
      fprintf(pFile, "Observation    ");
      pGroup->Write(pFile, WRITE_BNR);
      fprintf(pFile, "\n");
      ii = 0;
      for(i = 0; i < m_NumObs; i++)
      {
         pObs = pObsGroup->GetObsPtr(i);
         fprintf(pFile, "%-12s  ", pObs->GetName());
         jj = 0;
         for(j = 0; j < m_NumParams; j++)
         {
            if((m_bHoldObs[i] == false) && (m_bHoldParam[j] == false)){
               fprintf(pFile, "%+E ", m_pJacobUW[ii][jj]);}
            else{
               fprintf(pFile, "%+E ", 0.00);}

            if(m_bHoldParam[j] == false){ jj++;}
         }/* end for() */
         fprintf(pFile, "\n");
         if(m_bHoldObs[i] == false){ ii++;}
      }/* end for() */

      fprintf(pFile, "\nNormal Matrix\n");
      for(i = 0; i < np; i++)
      {
         for(j = 0; j < np; j++)
         {
            fprintf(pFile, "%+E  ", m_pNormal[i][j]);
         }/* end for() */
         fprintf(pFile, "\n");
      }/* end for() */
      if(m_bInv == true)
      {
         fprintf(pFile, "\nInverse Normal Matrix\n");
         for(i = 0; i < np; i++)
         {
            for(j = 0; j < np; j++)
            {
               fprintf(pFile, "%+E  ", m_pInvNormal[i][j]);
            }/* end for() */
            fprintf(pFile, "\n");
         }/* end for() */
      }/* end if() */
   }/* end if() */

   if(m_pPredictions != NULL)
   {
      int nrv = m_pPredictions->GetNumRespVars();
      UnchangeableString name;
      if(nrv > 0)
      {
         fprintf(pFile, "\nLinear Confidence Intervals on Predictions (%.2lf%%)\n", m_CIpct);
         fprintf(pFile, "Prediction       Expected Value  Std. Deviation  Lower Limit     Upper Limit\n");
         for(i = 0; i < nrv; i++)
         {
            name = m_pPredictions->GetRespVarPtr(i)->GetName();
            if(m_PredSD != NULL)
               fprintf(pFile, "%-15s  %E  %E  %E  %E\n", name, m_Pred[i], m_PredSD[i], m_PredCIlwr[i], m_PredCIupr[i]);
            else
               fprintf(pFile, "%-15s  %E  not_computed    not_computed    not_computed\n", name, m_Pred[i]);
         }/* end for() */
      }/* end if() */
   }/* end if() */
} /* end WriteStats() */

/******************************************************************************
WriteResiduals()

Write residuals at a given step or iteration to a file. This will be the current 
best set of residuals discovered by a given processor.
******************************************************************************/
void StatsClass::WriteResiduals(int step, char * prefix)
{
   ParameterGroup * pParamGroup;
   ObservationGroup * pObsGroup;
   Observation * pObs;
   int i, n, rank, t1, t2;
   double fcurbest, fprevbest, fcur;
   double m, p, w, d; //measured, predicted, weight, difference
   char fname[1000];
   char pname[1000];
   char cmd[1000];
   FILE * pFile;

   n = m_NumObs - m_NumHeldObs;
   pObsGroup = m_pModel->GetObsGroupPtr();
   pParamGroup = m_pModel->GetParamGroupPtr();

   if(m_bNoStats) return;

   if(m_bWriteIterationResiduals == false) return;

   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   sprintf(fname, "OstResiduals%s_P%03d_S%03d.txt", prefix, rank, step);
   pFile = fopen(fname, "r");

   fcurbest = HUGE_VAL;
   fprevbest = HUGE_VAL;
   fcur = m_pModel->GetObjFuncVal();

   //extract current best if file already exists
   if (pFile != NULL)
   {
      fscanf(pFile, "Iteration/Step : %d\n", &t1);
      fscanf(pFile, "Processor/Rank : %d\n", &t2);
      fscanf(pFile, "Min WSSE       : %lf\n", &fcurbest);
      fclose(pFile);
   }
   //otherwise extract best from previous step, if it exists
   else
   {
      sprintf(pname, "OstResiduals%s_P%03d_S%03d.txt", prefix, rank, step-1);
      pFile = fopen(pname, "r");

      if(pFile != NULL)
      {
         fscanf(pFile, "Iteration/Step : %d\n", &t1);
         fscanf(pFile, "Processor/Rank : %d\n", &t2);
         fscanf(pFile, "Min WSSE       : %lf\n", &fprevbest);
         fclose(pFile);
      }
   }

   // currently recorded value is already best
   if (fcurbest < fcur)
   {
      return;
   }
   // previously recorded value is best
   else if (fprevbest < fcur) 
   {
      #ifdef WIN32
         sprintf(cmd, "copy %s %s", pname, fname);
      #else
         sprintf(cmd, "cp %s %s", pname, fname);
      #endif
      system(cmd);

      return;
   }
   else
   {
      pFile = fopen(fname, "w");
      if(pFile == NULL) return;
   
      fprintf(pFile, "Iteration/Step : %d\n", step);
      fprintf(pFile, "Processor/Rank : %d\n", rank);
      fprintf(pFile, "Min WSSE       : %E\n", fcur);
      fprintf(pFile, "\nParameter Values\n");
      pParamGroup->Write(pFile, WRITE_OPT);
      fprintf(pFile, "\n\nObservation Residuals\n");
      fprintf(pFile, "\nObservation   Measured      Simulated     Weight        Residual        ");
      if(m_bOkToHoldObs == true)
      {
         fprintf(pFile, "   Sensitive?");
      }
      fprintf(pFile, "\n");

      for(i = 0; i < m_NumObs; i++)
      {
         pObs = pObsGroup->GetObsPtr(i);
         m = pObs->GetMeasuredVal(false, false);
         p = pObs->GetComputedVal(false, false);
         w = GetObsWeight(pObs);
         d = pObs->CalcResidual(true, true);
         fprintf(pFile, "%-12s  %E  %E  %E  %+E   ", pObs->GetName(), m, p, w, d);
         if(m_bOkToHoldObs == true)
         {
            if(m_bHoldObs[i] == false){ fprintf(pFile, "   YES");}
            else                      { fprintf(pFile, "   NO");}
         }
         fprintf(pFile, "\n");
      }/* end for() */   
      fclose(pFile);
   }
}/* end WriteResiduals() */

/******************************************************************************
STATS_Program()

Compute the statisitics of the parameter set defined in the input file.
******************************************************************************/
void STATS_Program(int argc, StringType argv[])
{
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model;

   NEW_PRINT("StatsClass", 1);
   StatsClass * stats = new StatsClass(model);
   MEM_CHECK(stats);
   RegisterStatsPtr(stats);

   FILE * pFile;
   int id;
   char outName[DEF_STR_SZ];

   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   sprintf(outName, "OstOutput%d.txt", id);

   //write setup to file
   pFile = fopen(outName, "w");
   fprintf(pFile, "Ostrich Setup\n");
   fprintf(pFile, "Model: %s\n", model->GetModelStr());
   fprintf(pFile, "Algorithm: Regression Statistics\n");
   fprintf(pFile, "Objective Function: %s\n", model->GetObjFuncStr());
   fprintf(pFile, "Number of Parameters: %d\n", model->GetParamGroupPtr()->GetNumParams());
   fprintf(pFile, "Number of Observations: ");
   if(model->GetObsGroupPtr() == NULL){fprintf(pFile, "0\n");}
   else {fprintf(pFile, "%d\n", model->GetObsGroupPtr()->GetNumObs());}
   fclose(pFile);

   //write setup to stdout
   fprintf(stdout, "Ostrich Setup\n");
   fprintf(stdout, "Model: %s\n", model->GetModelStr());
   fprintf(stdout, "Algorithm: Regression Statistics\n");
   fprintf(stdout, "Objective Function: %s\n", model->GetObjFuncStr());
   fprintf(stdout, "Number of Parameters: %d\n", model->GetParamGroupPtr()->GetNumParams());
   fprintf(stdout, "Number of Observations: ");
   if(model->GetObsGroupPtr() == NULL){fprintf(stdout, "0\n");}
   else {fprintf(stdout, "%d\n", model->GetObsGroupPtr()->GetNumObs());}
   
   stats->CalcStats();
   stats->WriteStats(stdout);
   
   pFile = fopen(outName, "a");
   stats->WriteStats(pFile);
   fclose(pFile);

   delete stats;
   delete model;
}/* end STATS_Program() */

/******************************************************************************
Jacobian_Program()

Compute the jacobian of the parameter set defined in the input file.
******************************************************************************/
void Jacobian_Program(int argc, StringType argv[])
{
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model;

   NEW_PRINT("StatsClass", 1);
   StatsClass * stats = new StatsClass(model);
   MEM_CHECK(stats);
   RegisterStatsPtr(stats);

   FILE * pFile;
   int id, j;
   double negOne = -1.00;
   char outName[DEF_STR_SZ];
   char tmp[DEF_STR_SZ];
   const char * inFile = GetOstFileName();
   char * line, * pTok;

   //allocate space for the parameter list
   int num = model->GetParamGroupPtr()->GetNumParams();
   double * pVals = new double[num];

   /* read in user-specified parameter set */
   pFile = fopen(inFile, "r");
   FindToken(pFile, "BeginInitParams", inFile);
   line = GetNxtDataLine(pFile, inFile);

   pTok = line;
   //extract values, one-by-one, making any necessary conversions
   for(int k = 0; k < num; k++)
   {
      j = ExtractString(pTok, tmp);
      j = ValidateExtraction(j, k, num, "Jacobian_Program()");
      pTok += j;            
      pVals[k] = model->GetParamGroupPtr()->GetParamPtr(k)->ConvertInVal(atof(tmp));
   }/* end for() */
   model->GetParamGroupPtr()->WriteParams(pVals);
   delete [] pVals;

   FindToken(pFile, "EndInitParams", inFile);
   fclose(pFile);


   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   sprintf(outName, "OstOutput%d.txt", id);

   //write setup to file
   pFile = fopen(outName, "w");
   fprintf(pFile, "Ostrich Setup\n");
   fprintf(pFile, "Model: %s\n", model->GetModelStr());
   fprintf(pFile, "Algorithm: Jacobian Calculation\n");
   fprintf(pFile, "Objective Function: %s\n", model->GetObjFuncStr());
   fprintf(pFile, "Number of Parameters: %d\n", model->GetParamGroupPtr()->GetNumParams());
   fprintf(pFile, "Number of Observations: ");
   if(model->GetObsGroupPtr() == NULL){fprintf(pFile, "0\n");}
   else {fprintf(pFile, "%d\n", model->GetObsGroupPtr()->GetNumObs());}
   fprintf(pFile, "Jacobian Matrix written to OstJacobian.txt\n");
   fclose(pFile);

   //write setup to stdout
   fprintf(stdout, "Ostrich Setup\n");
   fprintf(stdout, "Model: %s\n", model->GetModelStr());
   fprintf(stdout, "Algorithm: Jacobian Calculation\n");
   fprintf(stdout, "Objective Function: %s\n", model->GetObjFuncStr());
   fprintf(stdout, "Number of Parameters: %d\n", model->GetParamGroupPtr()->GetNumParams());
   fprintf(stdout, "Number of Observations: ");
   if(model->GetObsGroupPtr() == NULL){fprintf(stdout, "0\n");}
   else {fprintf(stdout, "%d\n", model->GetObsGroupPtr()->GetNumObs());}
   fprintf(stdout, "Jacobian Matrix written to OstJacobian.txt\n");

   model->Execute();
   Unchangeable2DArray pJ = stats->CalcJacobian(false, false, &negOne); //compute Jacobian, possibly in parallel

   if(id == 0)
   {
     FILE * pOut = fopen("OstJacobian.txt", "w");
     for(int i = 0; i < model->GetObsGroupPtr()->GetNumObs(); i++)
     {
       for(int j = 0; j < model->GetParamGroupPtr()->GetNumParams(); j++)
       {
         fprintf(pOut, "%.14E ", pJ[i][j]);
       }/* end for() */
       fprintf(pOut, "\n");
     }/* end for() */
     fclose(pOut);
   }/* end if() */
   
   pFile = fopen(outName, "a");
   stats->WriteMetrics(pFile);
   fclose(pFile);
   stats->WriteMetrics(stdout);
   
   delete stats;
   delete model;
}/* end Jacobian_Program() */

/******************************************************************************
EVAL_Program()

Evaluate objective function using a list of predefined parameter values.
******************************************************************************/
void EVAL_Program(int argc, StringType argv[])
{
   int samples_per_iter, num_left, count;
   double viol;
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model;

   FILE * pFile;
   int i, j, k, id, size, num, bi, np;
   double ** pList;
   double val, best;
   char tmp[DEF_STR_SZ];
   const char * inFile = GetOstFileName();
   char * line, * pTok;

   /* initialize parameter sets to specied values */
   pFile = fopen(inFile, "r");
   FindToken(pFile, "BeginInitParams", inFile);
   FindToken(pFile, "EndInitParams", inFile);
   rewind(pFile);

   //allocate space for the parameter list
   num = model->GetParamGroupPtr()->GetNumParams();

   //count the number of entries
   FindToken(pFile, "BeginInitParams", inFile);
   line = GetNxtDataLine(pFile, inFile);
   size = 0;
   while(strstr(line, "EndInitParams") == NULL)
   {
      size++;
      line = GetNxtDataLine(pFile, inFile);
   }/* end while() */

   //allocate space for entries
   if(size > 0)
   {
      NEW_PRINT("double *", size);
      pList = new double * [size];
      MEM_CHECK(pList);
      for(i = 0; i < size; i++)
      { 
         NEW_PRINT("double", num+1);
         pList[i] = new double[num+1];
         MEM_CHECK(pList[i]);
      }
   }/* end if() */

   //read in entries
   rewind(pFile);
   FindToken(pFile, "BeginInitParams", inFile);
   line = GetNxtDataLine(pFile, inFile);
   i = 0;
   while(strstr(line, "EndInitParams") == NULL)
   {
      pTok = line;
      //extract values, one-by-one, making any necessary conversions
      for(k = 0; k < num; k++)
      {
         j = ExtractString(pTok, tmp);
         j = ValidateExtraction(j, k, num, "EVAL_Program()");
         pTok += j;            
         pList[i][k] = model->GetParamGroupPtr()->GetParamPtr(k)->ConvertInVal(atof(tmp));
      }/* end for() */                  
      i++;
      line = GetNxtDataLine(pFile, inFile);
   }/* end while() */

   /*
   ----------------------------------------------------------------------
   Read in flag to use penalty function for infeasible parameter settings
   ----------------------------------------------------------------------
   */   
   bool bUsePenalty = true; //default is to use the penalty
   char tmp1[DEF_STR_SZ];
   char tmp2[DEF_STR_SZ];
   rewind(pFile);
   double penalty = -1.00;
   if(CheckToken(pFile, "PenalizeInfeasibleParameters", inFile) == true)
   {   
      line = GetCurDataLine();
      sscanf(line, "%s %s %lf", tmp1, tmp2, &penalty);
      MyStrLwr(tmp2);
      if(strncmp(tmp2, "no", 2) == 0) {bUsePenalty = false;}
   }/* end if() */
   if(penalty < 0.00) penalty = 1.00;

   fclose(pFile);
   
   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   MPI_Comm_size(MPI_COMM_WORLD, &np);

   if(id == 0)
   {
      WriteSetup(model, "Model Evaluations");
      //write banner
      WriteBanner(model, "iter   best value     ", "Percent Complete");
   }/* end if() */

   //insert warm start solution, if desired
   int istart = 0;
   if (model->CheckWarmStart() == true)
   {
      istart = ResumeEvaluations(model, id, np, pList[0]);
   }

   if(np == 1) //serial execution
   {      
      //perform model evaluations
      if(istart == 0)
      {
         bi = 0;
         best = NEARLY_HUGE;
      }
      else
      {
         bi = 0;
         best = pList[0][num];
      }
      
      count = 0;
      num_left = size - istart;
      for(i = istart; i < size; i++)
      {
         if(count == 0)
         {
            if(num_left > 10) samples_per_iter = 10;
            else samples_per_iter = num_left;
            WriteInnerEval(WRITE_USR, samples_per_iter, '.');
         }

         viol = model->GetParamGroupPtr()->WriteParams(pList[i]);
         if(!bUsePenalty) viol = 0.00;
         else viol *= penalty;
         val = ((Model *)model)->Execute(viol);
         num_left--;
         WriteInnerEval(++count, 0, '.');

         if(val < best)
         { 
            best = val;
            bi = i;
         }

         if(count == samples_per_iter)
         {
            count = 0;
            WriteInnerEval(WRITE_ENDED, 0, '.');

            model->GetParamGroupPtr()->WriteParams(pList[bi]);
            WriteRecord(model, i+1, best, 100.00*(1.00 - (double)num_left/(double)size));
         }
      }
   }/* end if() */
   else /* parallel execution */
   {
      i = istart;
      num_left = size - istart;
      while(num_left > 0)
      {
         bi = EvalInitParamsParallel(np, id, pList, size, model, &num_left);
         if(id == 0)
         {
            model->GetParamGroupPtr()->WriteParams(pList[bi]);
            best = pList[bi][num];
            WriteRecord(model, i+1, best, 100.00*(1.00 - (double)num_left/(double)size));
         }
         i++;
      }
   }/* end else() */
 
   for(i = 0; i < size; i++)
   {
      delete [] pList[i];
   }
   delete [] pList;

   delete model;
}/* end EVAL_Program() */

/******************************************************************************
EvalInitParamsParallel()

Compute objective function of entire set of samples in parallel. Each processor 
evaluates a predetermined number of samples, based on their processor id.

Returns index of the best (i.e. lowest objective function) parameter set.
******************************************************************************/
int EvalInitParamsParallel(int np, int id, double ** pList, int size, 
                           ModelABC * pModel, int * num_left)
{  
   double * F = new double[np];
   double Fx;
   int i ,j, samples_per_iter;
   int num = pModel->GetParamGroupPtr()->GetNumParams();
   ParameterGroup * pGroup;
  
   if(*num_left > np) samples_per_iter = np;
   else samples_per_iter = *num_left;

   //perform parallel evaluations
   j = size - *num_left;
   pGroup = pModel->GetParamGroupPtr();
   for(i = 0; i < samples_per_iter; i++) 
   { 
      if((i % np) == id)
      {
         pGroup->WriteParams(pList[i+j]);
         F[i] = pList[i+j][num] = pModel->Execute();
        }/* end if() */
   }/* end for() */

   *num_left = *num_left - samples_per_iter;

   //gather results
   for(i = 0; i < samples_per_iter; i++)
   {
      //receive someones F(x)
      Fx = F[i];
      MPI_Bcast((void *)(&Fx), 1, MPI_DOUBLE, i, MPI_COMM_WORLD);
      pList[i+j][num] = Fx;
   }/* end for() */

   //determine the 'best'
   int istart = j;
   int bi = 0;
   double best = pList[0][num];
   j = size - *num_left;
   for(i = istart; i < j; i++)
   {
      if(pList[i][num] < best) 
      {
         bi = i;
         best = pList[i][num];
      }
   }/* end for() */

   delete [] F;
   return bi;
}/* end EvalInitParamsParallel() */

/******************************************************************************
ResumeEvaluations()

Read the solutions from a previous run. Returns index of next solution.
******************************************************************************/
int ResumeEvaluations(ModelABC * pModel, int id, int nprocs, double * pbest)
{
   int retval = 0;
   int np = pModel->GetParamGroupPtr()->GetNumParams();
   int newcount = SimpleWarmStart(np, pbest);
   ((Model *)pModel)->SetCounter(newcount);
   if(nprocs == 1)
   {
      retval = newcount;
   }
   else
   {
      MPI_Allreduce(&newcount, &retval, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   }
   return retval;
}/* end ResumeEvaluations() */
