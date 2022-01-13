/******************************************************************************
File      : ObjectiveFunction.cpp
Author    : L. Shawn Matott and Vijaykumar Raghavan
Copyright : 2003, L. Shawn Matott and Vijaykumar Raghavan

Computes the objective function, which can either be weighted sum of squared 
errors (WSSE) or sum of the absolute weighted error (SAWE).

Version History
06-12-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
03-05-04    lsm   added GetObjFuncStr()
08-17-04    lsm   added reporting of memory allocations
12-18-04    lsm   added more descriptive error messages to USER obj. func.
01-10-05    lsm   removed some unused member variables
******************************************************************************/
#include <string.h>
#include <math.h>

#include "ObjectiveFunction.h"
#include "ObservationGroup.h"
#include "Observation.h"

#include "Exception.h"
#include "Utility.h"

/* ---------------------------------------------------------------------------
Box-Cox variables are set by WSSE constructor and shared with other modules
so that they can transform residuals as needed.
--------------------------------------------------------------------------- */
double gBoxCoxParam = 1.00;
bool   gBoxCoxFlag = false;

/******************************************************************************
WSSE::CTOR

Sets the observation group pointer.
******************************************************************************/
WSSE::WSSE(ObservationGroup * pObsGroup, bool boxCoxFlag, double boxCoxVal)
{
   gBoxCoxFlag = boxCoxFlag;
   gBoxCoxParam = boxCoxVal;
   m_pObsGroup = pObsGroup;
   strcpy(m_ObjFuncStr, "WSSE");
   IncCtorCount();
}/* end WSSE CTOR */

/******************************************************************************
WSSE::Destroy()
******************************************************************************/
void WSSE::Destroy(void)
{
   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
WSSE::GetCostFuncStr()
******************************************************************************/
UnchangeableString WSSE::GetCostFuncStr(int which)
{
   int nGroups = m_pObsGroup->GetNumGroups();
   if(nGroups <= 1) return m_ObjFuncStr;
   return(m_pObsGroup->GetGroup(which));
}/* end GetCostFuncStr() */

/******************************************************************************
WSSE::CalcMultiObjFunc()

Compute WSSE of each observation group.
******************************************************************************/
int WSSE::CalcMultiObjFunc(double * pF, int nObj)
{
   int nGroups = m_pObsGroup->GetNumGroups();
   if((pF == NULL) && (nObj == -1)) return nGroups;
   
   double error;
   double sum;
   int numObs;
   int i;
   Observation * pObs;
   numObs = m_pObsGroup->GetNumObs();

   for(int whichObj = 0; whichObj < nGroups; whichObj++)
   {
      sum = 0.00;
      UnchangeableString groupStr = m_pObsGroup->GetGroup(whichObj);

      for(i = 0; i < numObs; i++)
      {
         pObs = m_pObsGroup->GetObsPtr(i);
         if(strcmp(pObs->GetGroup(), groupStr) == 0)
         {
            error = pObs->CalcResidual(true, true);
            sum += (error * error);
         }
      } /* end for() */
      pF[whichObj] = sum;
   } /* end for() */
   return nGroups;
}/* end CalcMultiObjFunc() */

/******************************************************************************
WSSE::CalcObjFunc()

Computes the objective function and returns the result.
******************************************************************************/
double WSSE::CalcObjFunc(void)
{
   double error;
   double sum;
   int numObs;
   int i;
   Observation * pObs;

   sum = 0.00;
   numObs = m_pObsGroup->GetNumObs();

   for(i = 0; i < numObs; i++)
   {
      pObs = m_pObsGroup->GetObsPtr(i);
      error = pObs->CalcResidual(true, true);
      sum += (error * error);
   } /* end for() */
   return sum;
} /* end WSSE::CalcObjFunc() */

/******************************************************************************
WSSE::CalcObjFunc()

Computes the WSSE objective function without any transformation applied.
******************************************************************************/
double WSSE::CalcUntransformedObjFunc(void)
{
   double measured, computed, weight, error;
   double sum;
   int numObs;
   int i;
   Observation * pObs;

   sum = 0.00;
   numObs = m_pObsGroup->GetNumObs();

   for(i = 0; i < numObs; i++)
   {
      pObs = m_pObsGroup->GetObsPtr(i);
      measured = pObs->GetMeasuredVal(false, false);
      computed = pObs->GetComputedVal(false, false);
      weight = GetObsWeight(pObs);      
      error = weight*(measured - computed);
      sum += (error * error);
   } /* end for() */
   return sum;
} /* end WSSE::CalcUntransformedObjFunc() */

/******************************************************************************
BoxCox()

Perform Box-Cox transformation on the given input value.
******************************************************************************/
double BoxCox(double y)
{  
  if (gBoxCoxFlag == false) return y;

  //y must be positive, if not: don't perform transformation and log error
  if(y < 0.00)
  {
    LogError(ERR_BAD_ARGS, "Couldn't perform Box-Cox transformation, data is non-positive!");
    return y;
  }

  double lambda = gBoxCoxParam;
  double h = 0.00;
  if(lambda != 0)
  {
    h = (pow(y, lambda) - 1.00)/lambda;
  }
  else //natural log transformation
  {
    h = log(y);
  }
  return h;
} /* end BoxCox() */

/******************************************************************************
UnWeightJacobian()

Remove the weight term from a weighted jacobian entry.
******************************************************************************/
double UnWeightJacobian(double J, double w)
{  
  if (gBoxCoxFlag == false) return (J/w);

  if(gBoxCoxParam != 0)
  {
    return (J/(pow(w, gBoxCoxParam)));
  }
  else //natural log transformation, weights already removed from differencing
  {
    return J;
  }
} /* end UnWeightJacobian() */

/******************************************************************************
WSSE::WriteSetupToFile()

Output summary of setup.
******************************************************************************/
void WSSE::WriteSetupToFile(FILE * pFile)
{ 
   double b = gBoxCoxParam;
   if(gBoxCoxFlag == true)
   {
      fprintf(pFile, "WSSE calibration using a Transform-Both-Sides approach.\n");
      fprintf(pFile, "Box-Cox Parameter (b) : %lf\n", b);
      fprintf(pFile, "Box-Cox Formula : h(y,b) = ");
      if(b != 0.00)  
         fprintf(pFile, "(y^b - 1) / b\n");
      else
         fprintf(pFile, "log(y)\n");
   }
}/* end WriteSetupToFile() */

/******************************************************************************
SAWE::Destroy()
******************************************************************************/
void SAWE::Destroy(void)
{
   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
WSSE::CTOR

Sets the observation group pointer.
******************************************************************************/
SAWE::SAWE(ObservationGroup * pObsGroup)
{
  m_pObsGroup = pObsGroup;
  strcpy(m_ObjFuncStr, "SAWE");
  IncCtorCount();
}/* end SAWE CTOR */

/******************************************************************************
SAWE::CalcObjFunc()

Computes the objective function and returns the result.
******************************************************************************/
double SAWE::CalcObjFunc(void)
{
   int i;
   int numObs;
   double error;   
   double sumOfErrors;
   Observation * pObs;

   numObs = m_pObsGroup->GetNumObs();
   sumOfErrors  = 0.00;

   for (i = 0; i < numObs; i++)
   {
      pObs = m_pObsGroup->GetObsPtr(i);
      error = pObs->CalcResidual(true, true);
      if (error < 0) { error *= -1.0; }
      sumOfErrors += error;
   }/* end for() */

   return sumOfErrors;
} /* end SAWE::CalcObjFunc() */

/******************************************************************************
UserObjFunc::CTOR

Set the name of the output file of user-defined obj. function program, where
the objective function value is stored.
******************************************************************************/
UserObjFunc::UserObjFunc(IroncladString pFileName)
{
   int len;
   m_pObsGroup = NULL;
   strcpy(m_ObjFuncStr, "USER");

   len = (int)strlen(pFileName) + 1;
   NEW_PRINT("char", len);
   m_FileName = new char[len];
   MEM_CHECK(m_FileName);

   strcpy(m_FileName, pFileName);

   m_FileStr = NULL;

   IncCtorCount();
}/* end UserObjFunc CTOR */

/******************************************************************************
UserObjFunc::Destroy()
******************************************************************************/
void UserObjFunc::Destroy(void)
{
   delete [] m_FileName;
   delete [] m_FileStr;
   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
UserObjFunc::CalcObjFunc()

Computes the objective function and returns the result.
******************************************************************************/
double UserObjFunc::CalcObjFunc(void)
{
   UnchangeableString curPos; //current position within the file string
   char tmp[DEF_STR_SZ], valStr[DEF_STR_SZ];
   double val;
   int i;

   //read the output of the user-defined executable
   FileToString();
   
   //verify required output
   curPos = strstr(m_FileStr, "OST_ObjFuncVal");
   if(curPos == NULL)
   {      
      LogError(ERR_FILE_IO, "Couldn't locate OST_ObjFuncVal tag-string in model output");
      ExitProgram(1);
   }/* end if() */

   //extract the objective function value, use the last occurence
   valStr[0] = (char)NULL;
   curPos = m_FileStr;
   while(curPos != NULL)
   {
      curPos = strstr(curPos, "OST_ObjFuncVal");
      if(curPos != NULL) 
      { 
         sscanf(curPos, "%s %s", tmp, valStr);
         if(valStr[0] == (char)NULL)
         {
            LogError(ERR_FILE_IO, "Couldn't locate objective function value for model output");
            ExitProgram(1);
         }
         val = atof(valStr);
         curPos++;
      }
   }/* end while() */

   //check for model errors
   curPos = strstr(m_FileStr, "OST_ModelErrCode");
   if(curPos != NULL) 
   {      
      if(strstr(curPos, "no_errors") == NULL)
      {   
         for(i = 0; curPos[i] != '\n'; i++){ tmp[i] = curPos[i];}
         tmp[i] = (char)NULL;
         LogError(ERR_MODL_EXE, tmp);
      }/* end if() */
   }/* end if() */

   return val;     
} /* end UserObjFunc::CalcObjFunc() */

/******************************************************************************
FileToString()

Reads a file into a string.
******************************************************************************/
void UserObjFunc::FileToString(void)
{
   int fileSize;
   int i;
   FILE * pFile;

   pFile = fopen(m_FileName, "r");
   if(pFile == NULL)
   {
      FileOpenFailure("UserObjFunc::FileToString", m_FileName);
   }/* end if() */

   /*
   count number of chars in file, 
   so that fileStr can be sized
   */
   fileSize = 0;
   while(feof(pFile) == 0) 
   {
      fileSize++;
      fgetc(pFile);
   }/* end while() */   
   fileSize--;

   //size fileStr
   if(m_FileStr != NULL)
   {
      delete [] m_FileStr;
   }/* end if() */
   NEW_PRINT("char", fileSize+1);
   m_FileStr = new char[fileSize+1];
   MEM_CHECK(m_FileStr);

   //fill fileStr
   rewind(pFile);
   for(i = 0; i < fileSize; i++)
   {
      m_FileStr[i] = (char)(fgetc(pFile));
   }/* end for() */
   m_FileStr[i] = 0;

   fclose(pFile);
} /* end FileToString() */

