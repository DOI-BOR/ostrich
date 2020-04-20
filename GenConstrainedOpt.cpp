/******************************************************************************
File      : GenConstrainedOpt.cpp
Author    : L. Shawn Matott
Copyright : 2005, L. Shawn Matott

Defines a general constrained optimization extension to the ObjectiveFunction class.

This class supports a veriety of cost and constraint formulations, allowing users to
define fairly generic objective functions without having to write a separate driver 
program.

This class instantiates a set of constraint classes which can be 
combined with the system cost using a user-selected penalty method (additive penalty, 
multiplicative penalty, etc.). Cost and constraints are made up of response variables
which are functions of model output and/or model parameters.
   
Version History
01-10-05    lsm   created
******************************************************************************/
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "GenConstrainedOpt.h"
#include "ObservationGroup.h"
#include "ResponseVarGroup.h"
#include "RespVarABC.h"
#include "ParameterGroup.h"
#include "ObjectiveFunction.h"
#include "ConstraintABC.h"

#include "Exception.h"
#include "Utility.h"

/*
A mapping between penalty methods and human readable strings.
*/
IroncladString PenMethMap[NUM_PEN_METHS] = 
{
   "Additive Penalty Method (APM)",
   "Multiplicative Penalty Method (MPM)",
   "Exponential Penalty Method (EPM)"
};
//provides access to the mapping
IroncladString GetPenMethStr(LmtPenType i){ return PenMethMap[i];}

/******************************************************************************
GetResponseVarGroup()
******************************************************************************/
void * GCOP::GetResponseVarGroup(void) 
{ 
	return (void *)m_pRespGroup; 
}/* end GetResponseVarGroup() */

/******************************************************************************
GCOP::CTOR

Sets up the general constrained optimizer.
******************************************************************************/
GCOP::GCOP(ParameterGroup * pParamGroup)
{
   m_pObsGroup = NULL;
   m_pParamGroup = pParamGroup;
   strcpy(m_ObjFuncStr, "GCOP");
   m_PenType = PEN_TYPE_MPM;
   m_pConstraints = NULL;
   m_pRespGroup = NULL;      

   InitFromFile();

   IncCtorCount();
}/* end GCOP CTOR */

/******************************************************************************
GCOP::WriteSetupToFile()

Output summary of setup.
******************************************************************************/
void GCOP::WriteSetupToFile(FILE * pFile)
{
   int count;
   ConstraintABC * pCur;

   //count constraints
   count = 0;
   pCur = m_pConstraints;
   while(pCur != NULL){count++; pCur = pCur->GetNext();}
   
   fprintf(pFile, "Number of Resp. Vars        : %d\n", m_pRespGroup->GetNumRespVars());
   fprintf(pFile, "Number of Tied Resp. Vars   : %d\n", m_pRespGroup->GetNumTiedRespVars());
   fprintf(pFile, "Number of Constraints       : %d\n", count);
   fprintf(pFile, "Penalty Method              : %s\n", GetPenMethStr(m_PenType));
}/* end WriteSetupToFile() */

/******************************************************************************
GCOP::InitFromFile()

Initialize the GCOP classes by parsing the information in the input file.
******************************************************************************/
void GCOP::InitFromFile(void)
{
   const char * start_tag = "BeginGCOP";
   const char * end_tag   = "EndGCOP";
   int    start_len = (int)strlen(start_tag);
   int    end_len   = (int)strlen(end_tag);
   int whichObj;

   FILE * pFile;
   char * lineStr;
   char tmp1[DEF_STR_SZ], tmp2[DEF_STR_SZ], costStr[DEF_STR_SZ];
   IroncladString fileName = GetInFileName();

   //read in response variables
   InitResponseVars();

   pFile = fopen(fileName, "r");
   costStr[0] = 0;
   
   if(pFile == NULL)
   {
      FileOpenFailure("GCOP::InitFromFile", fileName);
   }/* end if() */

   //make sure correct tokens are present
   FindToken(pFile, start_tag, fileName);
   FindToken(pFile, end_tag, fileName);
   rewind(pFile);

   //count the number of cost functions
   FindToken(pFile, start_tag, fileName);
   lineStr = GetNxtDataLine(pFile, fileName);
   m_NumMultiObjCostFuncs = 0;
   while(strncmp(lineStr, end_tag, end_len) != 0)
   {      
      if(strncmp(lineStr, "CostFunction", 12) == 0)
      {
         m_NumMultiObjCostFuncs++;
      }/* end if() --> Cost Function */
      lineStr = GetNxtDataLine(pFile, fileName);
   }/* end while() */
   rewind(pFile);

   if(m_NumMultiObjCostFuncs == 0)
   {
      LogError(ERR_FILE_IO, "No Cost Function was defined");
      ExitProgram(1);
   }/* end if() */

   m_pMultiObjCostFunc = NULL;
   if(m_NumMultiObjCostFuncs > 0)
   {
      m_pMultiObjCostFunc = new RespVarABC *[m_NumMultiObjCostFuncs];
      for(whichObj = 0; whichObj < m_NumMultiObjCostFuncs; whichObj++)
      {
         m_pMultiObjCostFunc[whichObj] = NULL;
      }/* end for() */
   }/* end if() */

   FindToken(pFile, start_tag, fileName);
   lineStr = GetNxtDataLine(pFile, fileName);
   whichObj = 0;
   while(strncmp(lineStr, end_tag, end_len) != 0)
   {      
      //read in configuration parameters
      //read in penalty function type
      if(strncmp(lineStr, "PenaltyFunction", 15) == 0)
      {
         sscanf(lineStr, "%s %s", tmp1, tmp2);
         MyStrLwr(tmp2);
         if(strcmp(tmp2, "apm") == 0)
            { m_PenType = PEN_TYPE_APM;}
         else if(strcmp(tmp2, "mpm") == 0)
            { m_PenType = PEN_TYPE_MPM;}
         else if(strcmp(tmp2, "epm") == 0)
            { m_PenType = PEN_TYPE_EPM;}
         else
         {
            sprintf(tmp1, "GCOP::InitFromFile() invalid Penalty Function: |%s|", tmp2);
            LogError(ERR_FILE_IO, tmp1);
         }         
      }/* end if() --> Penalty Function Type */

      else if(strncmp(lineStr, "CostFunction", 12) == 0)
      {
         sscanf(lineStr, "%s %s", tmp1, tmp2);
         strcpy(costStr, tmp2);

         m_pMultiObjCostFunc[whichObj] = m_pRespGroup->GetRespVarPtr(costStr);
         if(whichObj == 0)
         {
            m_pCostFunc = m_pMultiObjCostFunc[whichObj];
         }         
         if(m_pMultiObjCostFunc[whichObj] == NULL)
         {
            sprintf(tmp1, "GCOP::InitFromFile(): CostFunction |%s| is not a response variable", costStr);
            LogError(ERR_FILE_IO, tmp1);
            ExitProgram(1);
         }
         whichObj++;
      }/* end if() --> Cost Function */
      else
      {
         sprintf(tmp1, "GCOP::InitFromFile(): unknown token |%s|", lineStr);
         LogError(ERR_FILE_IO, tmp1);
      }

      lineStr = GetNxtDataLine(pFile, fileName);
   }/* end while() */

   fclose(pFile);

   //read in constraints
   InitConstraints();

   /* Display constraint information
   ConstraintABC * pCur;
   if(m_pConstraints != NULL)
   {
      pCur = m_pConstraints;
      pCur->Write(stdout, WRITE_BNR);
      fprintf(stdout, "\n");
      while(pCur != NULL)
      {
         pCur->Write(stdout, WRITE_SCI);
         fprintf(stdout, "\n");
         pCur = pCur->GetNext();
      }
   } */
}/* end InitFromFile()*/

/******************************************************************************
GCOP::InitResponseVars()

Initialize all response variables, which are the basis for the constraints.
******************************************************************************/
void GCOP::InitResponseVars(void)
{
   NEW_PRINT("ResponseVarGroup", 1);
   m_pRespGroup = new ResponseVarGroup();
   MEM_CHECK(m_pRespGroup);
}/* end InitResponseVars()*/

/******************************************************************************
GCOP::InitConstraints()

Initialize all constraints by parsing the information in the "Constraints" 
section of the input file.
******************************************************************************/
void GCOP::InitConstraints(void)
{
   FILE * pFile;
   StringType * names;
   RespVarABC * pLoc1;
   char * pTok, * pOld;
   int i, j, np, len;
   double conv, lwr, upr;
   char * lineStr, nameStr[DEF_STR_SZ], typeStr[DEF_STR_SZ];
   char tmp1[DEF_STR_SZ];
   IroncladString fileName = GetInFileName();

   pFile = fopen(fileName, "r");

   if(pFile == NULL)
   {
      FileOpenFailure("GCOP::InitConstraints", fileName);
   }/* end if() */

   if(CheckToken(pFile, "BeginConstraints", fileName) == true)
   {
      FindToken(pFile, "EndConstraints", fileName);
      rewind(pFile);

      FindToken(pFile, "BeginConstraints", fileName);
      lineStr = GetNxtDataLine(pFile, fileName);

      while(strcmp(lineStr, "EndConstraints") != 0)
      {
         sscanf(lineStr, "%s %s", nameStr, typeStr);
         MyStrLwr(typeStr);

         /*------------------------------------------------
         Capacity constraint
         ------------------------------------------------*/
         if(strcmp(typeStr, "capacity") == 0)
         {
            pTok = lineStr;
            //extract name of constraint (no spaces allowed)
            j = ExtractString(pTok, nameStr);
            j = ValidateExtraction(j, 1, 1, "GCOP::InitConstraints()");
            pTok += j;
            //extract type  
            j = ExtractString(pTok, typeStr);
            j = ValidateExtraction(j, 1, 1, "GCOP::InitConstraints()");
            pTok += j;
            //extract conversion factor
            j = ExtractString(pTok, tmp1);
            j = ValidateExtraction(j, 1, 1, "GCOP::InitConstraints()");
            conv = atof(tmp1);
            pTok += j;
            //extract lower bound
            j = ExtractString(pTok, tmp1);
            j = ValidateExtraction(j, 1, 1, "GCOP::InitConstraints()");
            lwr = atof(tmp1);
            pTok += j;
            //extract upper bound
            j = ExtractString(pTok, tmp1);
            j = ValidateExtraction(j, 1, 1, "GCOP::InitConstraints()");
            upr = atof(tmp1);
            pTok += j;
            pOld = pTok;
            //count parameter names
            np = 0;
            while(*pTok != (char)NULL){ if(*pTok == ','){ np++;} pTok++;}
            if(*(pTok-1) != ','){ np++;} //in case the last name isn't terminated by a comma
            NEW_PRINT("StringType", np);
            names = new StringType[np];
            MEM_CHECK(names);
            for(i = 0; i < np; i++){names[i] = NULL;}

            //read parameter name list
            pTok = pOld;
            for(i = 0; i < np; i++)
            {               
               j = ExtractColString(pTok, tmp1, ',');
               j = ValidateExtraction(j, i, np, "GCOP::InitConstraints()");
               MyTrim(tmp1);
               pTok += j;
               len = (int)strlen(tmp1)+1;
               NEW_PRINT("char", len);
               names[i] = new char[len];
               MEM_CHECK(names[i]);
               strcpy(names[i], tmp1);
            }
            //create constraint
            NEW_PRINT("CapacityConstraint", 1);
            CapacityConstraint * pNewCC = new CapacityConstraint(nameStr, names, np, m_pParamGroup, lwr, upr, conv);
            MEM_CHECK(pNewCC);
            //insert into linked list
            if(m_pConstraints == NULL){ m_pConstraints = pNewCC;}
            else{ m_pConstraints->AddConstraint(pNewCC);}

            //free up names list
            for(i = 0; i < np; i++){ delete [] names[i];}
            delete [] names;
         }/* end if() ---> Capacity Constraint */

         /*------------------------------------------------
         General constraint
         ------------------------------------------------*/
         else if (strcmp(typeStr, "general") == 0)
         {
            sscanf(lineStr, "%s %s %lf %lf %lf %s", nameStr, typeStr, &conv, &lwr, &upr, tmp1);
            pLoc1 = m_pRespGroup->GetRespVarPtr(tmp1);
            if(pLoc1 == NULL)
            {
               sprintf(typeStr, "GCOP::InitConstraints() unknown response variable |%s|", tmp1);
               LogError(ERR_FILE_IO, typeStr);
               ExitProgram(1);
            }/* end if() */
            //create constraint
            NEW_PRINT("GeneralConstraint", 1);
            GeneralConstraint * pNewGN = new GeneralConstraint(nameStr, pLoc1, lwr, upr, conv);
            MEM_CHECK(pNewGN);
            //insert into linked list
            if(m_pConstraints == NULL){ m_pConstraints = pNewGN;}
            else{ m_pConstraints->AddConstraint(pNewGN);}
         }/* end else if() ---> General Constraint */
         else
         {
            sprintf(nameStr, "GCOP::InitConstrints() unknown type |%s|", typeStr);
            LogError(ERR_FILE_IO, nameStr);
         }         
         lineStr = GetNxtDataLine(pFile, fileName);
      }/* end while() --> read in constraints */
   } /* end if() */
   else
   {
      LogError(ERR_FILE_IO, "No constraints specified.");
   }/* end else() */

   fclose(pFile);
}/* end InitConstraints()*/

/******************************************************************************
GCOP::Destroy()
******************************************************************************/
void GCOP::Destroy(void)
{
   delete m_pConstraints;
   delete m_pRespGroup;
   if(m_pMultiObjCostFunc!= NULL)
      delete [] m_pMultiObjCostFunc;
   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
GCOP::CalcObjFunc()

Computes the objective function and returns the result.
******************************************************************************/
double GCOP::CalcObjFunc(void)
{
   static bool firstTime = true;
   FILE * pFile;
   double cost = 0.00;
   double penalty = 0.00;
  
   ConstraintABC * pCur;
   
   m_pRespGroup->ExtractVals();
   
   if(firstTime == true)
   {
      pFile = fopen("OstGcopOut.txt", "w");
      fprintf(pFile, "True Cost \tPenalty \tAdjusted Cost\n");
      firstTime = false;
   }
   else
   { 
      pFile = fopen("OstGcopOut.txt", "a+");
   }

   //compute and output the un-penalized cost
   cost = m_pCostFunc->GetCurrentVal();
   fprintf(pFile, "%E\t", cost);

   /* compute constraint penalties */
   if(m_pConstraints != NULL)
   {
      pCur = m_pConstraints;
      while(pCur != NULL)
      {
         penalty += pCur->CalcPenalty();
         pCur = pCur->GetNext();
      }
   }
  
   fprintf(pFile, "%E\t", penalty);

   /* assess penalty using APM, MPM or EPM */
   if(penalty != 0.00)
   {
      switch(m_PenType)
      {
         case(PEN_TYPE_APM) :
         {
            cost += penalty;
            break;
         }
         case(PEN_TYPE_MPM) :
         {
            cost = MyMax(cost, penalty) * (1.00 + penalty);
            break;
         }         
         case(PEN_TYPE_EPM) :
         {
            if(penalty >= NEARLY_HUGE_LN_EXP){ cost = NEARLY_HUGE; break;}            
            cost = MyMax(cost, penalty) * exp(penalty);
            break;
         }         
      }/* end switch() */
   }/* end if() */

   fprintf(pFile, "%E\n", cost);
   fclose(pFile);
   return cost;
} /* end GCOP::CalcObjFunc() */

/******************************************************************************
GCOP::CalcMultiObjFunc()

Computes the multi-objective function vector and returns the result.

If pF == NULL and nObj == -1, just returns the number of objectives
******************************************************************************/
int GCOP::CalcMultiObjFunc(double * pF, int nObj)
{
   if((pF == NULL) && (nObj == -1)) return m_NumMultiObjCostFuncs;

   static bool firstTime = true;
   FILE * pFile;
   double cost = 0.00;
   double penalty = 0.00;
   int whichObj;
   ConstraintABC * pCur;
   char outfileName[DEF_STR_SZ];
   
   m_pRespGroup->ExtractVals();

   for(whichObj = 0; whichObj < nObj; whichObj++)
   {    
      sprintf(outfileName, "OstGcopOut_%s.txt", m_pMultiObjCostFunc[whichObj]->GetName());
      if(firstTime == true)
      {
         pFile = fopen(outfileName, "w");
         fprintf(pFile, "True Cost \tPenalty \tAdjusted Cost\n");         
      }
      else
      { 
         pFile = fopen(outfileName, "a+");
      }

      //compute and output the un-penalized cost
      cost = m_pMultiObjCostFunc[whichObj]->GetCurrentVal();
      fprintf(pFile, "%E\t", cost);

      /* compute constraint penalties */
      penalty = 0.00;
      if(m_pConstraints != NULL)
      {
         pCur = m_pConstraints;         
         while(pCur != NULL)
         {
            penalty += pCur->CalcPenalty();
            pCur = pCur->GetNext();
         }
      }/* end if() */
  
      fprintf(pFile, "%E\t", penalty);

      /* assess penalty using APM, MPM or EPM */
      if(penalty != 0.00)
      {
         switch(m_PenType)
         {
            case(PEN_TYPE_APM) :
            {
               cost += penalty;
               break;
            }
            case(PEN_TYPE_MPM) :
            {
               cost = MyMax(cost, penalty) * (1.00 + penalty);
               break;
            }         
            case(PEN_TYPE_EPM) :
            {
               if(penalty >= NEARLY_HUGE_LN_EXP){ cost = NEARLY_HUGE; break;}            
               cost = MyMax(cost, penalty) * exp(penalty);
               break;
            }         
         }/* end switch() */
      }/* end if() */

      fprintf(pFile, "%E\n", cost);
      fclose(pFile);
      pF[whichObj] = cost;
   }/* end for(whichObj) */
   firstTime = false;

   return m_NumMultiObjCostFuncs;
} /* end GCOP::CalcMultiObjFunc() */

/******************************************************************************
GCOP::WriteConstraints()

Display constraint information.
******************************************************************************/
void GCOP::WriteConstraints(FILE * pFile, int type)
{
   ConstraintABC * pCur;

   pCur = m_pConstraints;
   
   while(pCur != NULL)
   {
      if(type == WRITE_BNR)
      { 
         pCur->Write(pFile, type); 
         fprintf(pFile, "\n");
         return;
      }

      pCur->Write(pFile, type);
      fprintf(pFile, "\n");
      pCur = pCur->GetNext();            
   }/* end while() */
}/* end WriteConstraints() */

/******************************************************************************
GCOP::GetConstraintPtr()

Retrieve constraint associated with pName
******************************************************************************/
ConstraintABC * GCOP::GetConstraintPtr(IroncladString pName)
{
   ConstraintABC * pCur;

   for(pCur = m_pConstraints; pCur != NULL; pCur = pCur->GetNext())
   {
      if(strcmp(pCur->GetName(), pName) == 0) return pCur;
   }	
   return NULL;
}/* end GetConstraintPtr() */

