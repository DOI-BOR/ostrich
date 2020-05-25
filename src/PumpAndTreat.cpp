/****************************************************************************
File      : PumpAndTreat.cpp
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
   Capacity constraints which limit total of a set of parameters
   
Version History
05-07-04    lsm   created
01-10-05    lsm   Generalized the PatoConstraintABC class and modified to interface 
                  with abstract response variables (RespVarABC)
03-21-05    lsm   Added support for Mayer cost formulation.
******************************************************************************/
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "PumpAndTreat.h"
#include "ResponseVarGroup.h"
#include "RespVarABC.h"
#include "ConstraintABC.h"
#include "ObservationGroup.h"
#include "GenConstrainedOpt.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"

#include "Exception.h"
#include "Utility.h"

/*
A mapping between cost functions and human readable strings.
*/
IroncladString CostFuncMap[NUM_COST_FUNCS] = 
{
   "Minimize Total Q",
   "Minimize Operational Costs",
   "Minimize Capital and Operational Costs",
   "Minimize Capital, Operational and Treatment Costs"
};

/* Maintenance adjustment factors for 30 years from RS Means reference */
const double gAdjFact[30] = 
{
   0.01, 0.02, 0.01, 0.02, 0.05, 0.01, 0.02, 0.01, 0.02, 0.1, 
   0.01, 0.02, 0.01, 0.02, 0.05, 0.01, 0.02, 0.01, 0.02, 0.1, 
   0.01, 0.02, 0.01, 0.02, 0.05, 0.01, 0.02, 0.01, 0.02, 0.1
};

/*--------------------------------------------
Extraction well cost table RS Means reference 
0 : Qmin
1 : Qmax
2 : Lmin
3 : Lmax
4 : Cost
---------------------------------------------*/
const double gPumpCosts[77][5] = 
{
   {0.00, 7.00,   0.00, 140.00, 1828.00}, {0.00, 7.00, 140.00, 240.00, 1981.00},
   {0.00, 7.00, 240.00, 340.00, 2302.00}, {0.00, 7.00, 340.00, 520.00, 2604.00},
   {0.00, 7.00, 520.00, 800.00, 3409.00}, 
   
   {7.00, 14.00,   0.00,  80.00, 1584.00}, {7.00, 14.00,  80.00, 140.00, 1685.00},
   {7.00, 14.00, 140.00, 220.00, 1955.00}, {7.00, 14.00, 220.00, 280.00, 2179.00},
   {7.00, 14.00, 280.00, 460.00, 2613.00}, {7.00, 14.00, 460.00, 600.00, 3000.00},
   {7.00, 14.00, 600.00, 800.00, 4212.00},

   {14.00, 20.00,   0.00,   80.00, 1594.00}, {14.00, 20.00,  80.00, 160.00, 1824.00},
   {14.00, 20.00, 160.00,  240.00, 2013.00}, {14.00, 20.00, 240.00, 300.00, 2278.00},
   {14.00, 20.00, 300.00,  400.00, 2675.00}, {14.00, 20.00, 400.00, 600.00, 3566.00},
   {14.00, 20.00, 600.00, 1000.00, 4949.00}, 

   {20.00, 32.00,   0.00,  60.00, 1632.00}, {20.00, 32.00,  60.00, 120.00, 1803.00},
   {20.00, 32.00, 120.00, 160.00, 1940.00}, {20.00, 32.00, 160.00, 200.00, 2157.00},
   {20.00, 32.00, 200.00, 280.00, 2433.00}, {20.00, 32.00, 280.00, 340.00, 3245.00},
   {20.00, 32.00, 340.00, 600.00, 4346.00}, 

   {32.00, 55.00,   0.00,  20.00, 2031.00}, {32.00, 55.00,  20.00, 100.00, 2168.00},
   {32.00, 55.00, 100.00, 160.00, 2305.00}, {32.00, 55.00, 160.00, 220.00, 3148.00},
   {32.00, 55.00, 220.00, 340.00, 4451.00}, {32.00, 55.00, 340.00, 600.00, 5909.00},
   {32.00, 55.00, 600.00, 800.00, 8246.00}, 

   {55.00, 95.00,   0.00,  40.00, 2170.00}, {55.00, 95.00,  40.00, 100.00, 3042.00},
   {55.00, 95.00, 100.00, 220.00, 4113.00}, {55.00, 95.00, 220.00, 300.00, 5072.00},
   {55.00, 95.00, 300.00, 400.00, 6794.00}, 

   {95.00, 200.00,   0.00,   50.00,  2281.00}, {95.00, 200.00,   50.00,  100.00,  3779.00},
   {95.00, 200.00, 100.00,  150.00,  4481.00}, {95.00, 200.00,  150.00,  200.00,  6001.00},
   {95.00, 200.00, 200.00,  300.00,  7400.00}, {95.00, 200.00,  300.00,  400.00,  9360.00}, 
   {95.00, 200.00, 400.00,  500.00, 10505.00}, {95.00, 200.00,  500.00,  600.00, 13800.00},
   {95.00, 200.00, 600.00,  725.00, 16758.00}, {95.00, 200.00,  725.00,  950.00, 20507.00},
   {95.00, 200.00, 950.00, 1100.00, 35430.00}, {95.00, 200.00, 1100.00, 1400.00, 42382.00}, 

   {200.00, 410.00,   0.00,   75.00,  4067.00}, {200.00, 410.00,   75.00,  150.00,  6657.00},
   {200.00, 410.00, 150.00,  175.00,  8212.00}, {200.00, 410.00,  175.00,  225.00,  9076.00},
   {200.00, 410.00, 225.00,  300.00, 11989.00}, {200.00, 410.00,  300.00,  400.00, 13448.00}, 
   {200.00, 410.00, 400.00,  500.00, 15447.00}, {200.00, 410.00,  500.00,  600.00, 20787.00},
   {200.00, 410.00, 600.00,  750.00, 26495.00}, 

   {410.00, 680.00,   0.00,   50.00,   6527.00}, {410.00, 680.00,   50.00,  125.00,  8806.00},
   {410.00, 680.00, 125.00,  200.00,  12010.00}, {410.00, 680.00,  200.00,  275.00, 12972.00},
   {410.00, 680.00, 275.00,  350.00,  16018.00}, {410.00, 680.00,  350.00,  400.00, 18306.00}, 
   {410.00, 680.00, 400.00,  500.00,  23349.00}, {410.00, 680.00,  500.00,  700.00, 30504.00},
   {410.00, 680.00, 700.00,  900.00,  40264.00}, 

   {680.00, 1400.00,   0.00,  100.00,  13381.00}, {680.00, 680.00,  100.00,  175.00, 19466.00},
   {680.00, 1400.00, 175.00,  200.00,  23575.00}, {680.00, 680.00,  200.00,  225.00, 30772.00},
   {680.00, 1400.00, 225.00,  350.00,  35045.00}, {680.00, 680.00,  350.00,  400.00, 40218.00}, 
   {680.00, 1400.00, 400.00,  475.00,  51951.00}, {680.00, 680.00,  475.00,  600.00, 58854.00},
   {680.00, 1400.00, 600.00,  750.00,  66279.00}
};

typedef enum
{
   EXT_COST = 0,
   INJ_COST = 1,
   LAB_COST = 2,
   NRG_COST = 3,
   ANA_COST = 4,
   DIS_COST = 5,
   MNT_COST = 6,
   DRL_COST = 7,
   PMP_COST = 8,
   TRC_COST = 9,
   TRO_COST = 10
} CostIdxEnum;

/******************************************************************************
GetResponseVarGroup()
******************************************************************************/
void * PATO::GetResponseVarGroup(void) 
{ 
	return (void *)m_pRespGroup; 
}/* end GetResponseVarGroup() */

/******************************************************************************
PATO::CTOR

Sets the observation group pointer.
******************************************************************************/
PATO::PATO(ParameterGroup * pParamGroup)
{
   int i;
   m_pObsGroup = NULL;
   m_pParamGroup = pParamGroup;
   strcpy(m_ObjFuncStr, "PATO");
   m_ObjType = PATO_OBJ_RATE;
   m_PenType = PEN_TYPE_MPM;
   m_pConstraints = NULL;
   m_pRespGroup = NULL;
   m_pPlumes = NULL;
   m_pWells = NULL;
   m_NumPlumes = 0;
   m_MaxNumWells = 0;
   m_pTbl = NULL;
   m_TblSize = 0;

   m_ExtRateCF = 0.00;
   m_InjRateCF = 0.00;
   m_FixWellCF = 0.00;
   m_VarWellCF = 0.00;
   m_MayerDrillCF = 0.00;
   m_MayerPumpCF = 0.00;
   m_RateUCF = 0.00;
   m_LiftUCF = 0.00;
   m_TimeFrame = 0.00;
   m_IntRate = 0.00;
   m_LaborRate = 0.00;
   m_ExtEnergyRate = 0.00;
   m_InjEnergyRate = 0.00;
   m_AnalyticRate = 0.00;
   m_SampleFreq = 0.00;
   m_DisposalRate = 0.00;
   m_MaintFactor = 0.00;
   m_TreatCapCoeff = 0.00;
   m_TreatCapExpon = 1.00;
   m_TreatOpCoeff = 0.00;
   m_TreatOpExpon = 1.00;
   m_RateThresh   = 0.00;

   for(i = 0; i < TRO_COST + 1; i++){ m_Costs[i] = 0.00;}

   InitFromFile();

   IncCtorCount();
}/* end PATO CTOR */

/******************************************************************************
PATO::WriteSetupToFile()

Output summary of setup.
******************************************************************************/
void PATO::WriteSetupToFile(FILE * pFile)
{
   int count;
   ConstraintABC * pCur;

   //count constraints
   count = 0;
   pCur = m_pConstraints;
   while(pCur != NULL){count++; pCur = pCur->GetNext();}
   
   fprintf(pFile, "Number of Resp. Vars   : %d\n", m_pRespGroup->GetNumRespVars());
   fprintf(pFile, "Number of Constraints  : %d\n", count);
   fprintf(pFile, "Max. Number of Wells   : %d\n", m_MaxNumWells);
   fprintf(pFile, "Number of Plumes       : %d\n", m_NumPlumes);
   fprintf(pFile, "Cost Function          : %s\n", CostFuncMap[m_ObjType]);
   fprintf(pFile, "Penalty Method         : %s\n", GetPenMethStr(m_PenType));
}/* end WriteSetupToFile() */

/******************************************************************************
PATO::InitFromFile()

Initialize the PATO classes by parsing the information in the input file.
******************************************************************************/
void PATO::InitFromFile(void)
{
   FILE * pFile;
   char * lineStr;
   char tmp1[DEF_STR_SZ], tmp2[DEF_STR_SZ];
   IroncladString fileName = GetInFileName();

   pFile = fopen(fileName, "r");
   
   if(pFile == NULL)
   {
      FileOpenFailure("PATO::InitFromFile", fileName);
   }/* end if() */

   //make sure correct tokens are present
   FindToken(pFile, "BeginPumpAndTreat", fileName);
   FindToken(pFile, "EndPumpAndTreat", fileName);
   rewind(pFile);

   FindToken(pFile, "BeginPumpAndTreat", fileName);
   lineStr = GetNxtDataLine(pFile, fileName);

   while(strstr(lineStr, "EndPumpAndTreat") == NULL)
   {      
      //read in objective function type
      if(strstr(lineStr, "CostFunction") != NULL)
      {
         sscanf(lineStr, "%s %s", tmp1, tmp2);
         MyStrLwr(tmp2);
         if(strcmp(tmp2, "pumprate") == 0)
            { m_ObjType = PATO_OBJ_RATE;}
         else if(strcmp(tmp2, "opcost") == 0)
            { m_ObjType = PATO_OBJ_OP;}
         else if(strcmp(tmp2, "cap&opcost") == 0)
            { m_ObjType = PATO_OBJ_CAP_OP;}
         else if(strcmp(tmp2, "mayer") == 0)
            { m_ObjType = PATO_OBJ_MAYER;}
         else if(strcmp(tmp2, "cap&op&treatcost") == 0)
            { m_ObjType = PATO_OBJ_CAP_OP_TRE;}
         else
         {
            sprintf(tmp1, "PATO::InitFromFile() invalid Cost Function: |%s|", tmp2);
            LogError(ERR_FILE_IO, tmp1);
         }
      }/* end if() --> Objective Function Type */

      //read in penalty function type
      else if(strstr(lineStr, "PenaltyFunction") != NULL)
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
            sprintf(tmp1, "PATO::InitFromFile() invalid Penalty Function: |%s|", tmp2);
            LogError(ERR_FILE_IO, tmp1);
         }         
      }/* end else if() --> Penalty Function Type */

      //read in on/off threshold
      else if(strstr(lineStr, "OnOffThreshold") != NULL)
      {
         sscanf(lineStr, "%s %lf", tmp1, &m_RateThresh);
      }/* end else if() --> On/Off Threshold */

      //read in extraction rate cost coefficient
      else if(strstr(lineStr, "ExtRateCF") != NULL)
      {
         sscanf(lineStr, "%s %lf", tmp1, &m_ExtRateCF);
      }/* end else if() --> Extraction Rate Cost Coeff */

      //read in injection rate cost coefficient
      else if(strstr(lineStr, "InjRateCF") != NULL)
      {
         sscanf(lineStr, "%s %lf", tmp1, &m_InjRateCF);
      }/* end else if() --> Injection Rate Cost Coeff */

      //read in well construction fixed cost coefficient
      else if(strstr(lineStr, "FixedWellCF") != NULL)
      {
         sscanf(lineStr, "%s %lf", tmp1, &m_FixWellCF);
      }/* end else if() --> Well Construction Fixed Cost Coefficient */

      //read in well construction depth-dependent cost coefficient
      else if(strstr(lineStr, "DepthDepWellCF") != NULL)
      {
         sscanf(lineStr, "%s %lf", tmp1, &m_VarWellCF);
      }/* end else if() --> Well Construction Dependent Cost Coefficient */

      //read in drilling cost factor (for Mayer formulation)
      else if(strstr(lineStr, "MayerDrillCF") != NULL)
      {
         sscanf(lineStr, "%s %lf", tmp1, &m_MayerDrillCF);
      }/* end else if() --> Drill Cost Factor */

      //read in pumping cost factor (for Mayer formulation)
      else if(strstr(lineStr, "MayerPumpCF") != NULL)
      {
         sscanf(lineStr, "%s %lf", tmp1, &m_MayerPumpCF);
      }/* end else if() --> Pump Cost Factor */

      //read in extraction rate unit conversion factor
      else if(strstr(lineStr, "RateUCF") != NULL)
      {
         sscanf(lineStr, "%s %lf", tmp1, &m_RateUCF);
      }/* end else if() --> Extraction Rate Unit COnversion Factor */

      //read in lift unit conversion factor
      else if(strstr(lineStr, "LiftUCF") != NULL)
      {
         sscanf(lineStr, "%s %lf", tmp1, &m_LiftUCF);
      }/* end else if() --> Lift Unit COnversion Factor */

      //read in extraction energy cost rate
      else if(strstr(lineStr, "ExtEnergyRate") != NULL)
      {
         sscanf(lineStr, "%s %lf", tmp1, &m_ExtEnergyRate);
      }/* end else if() --> Extraction Energy Cost Rate */

      //read in injection energy cost rate
      else if(strstr(lineStr, "InjEnergyRate") != NULL)
      {
         sscanf(lineStr, "%s %lf", tmp1, &m_InjEnergyRate);
      }/* end else if() --> Injection Energy Cost Rate */

      //read in Labor Cost Rate
      else if(strstr(lineStr, "LaborRate") != NULL)
      {
         sscanf(lineStr, "%s %lf", tmp1, &m_LaborRate);
      }/* end else if() --> Labor Cost Rate */

      //read in Analytic cost rate
      else if(strstr(lineStr, "AnalyticRate") != NULL)
      {
         sscanf(lineStr, "%s %lf", tmp1, &m_AnalyticRate);
      }/* end else if() --> Analytic Cost Rate */

      //read in sample frequency
      else if(strstr(lineStr, "SampleFreq") != NULL)
      {
         sscanf(lineStr, "%s %lf", tmp1, &m_SampleFreq);
      }/* end else if() --> Sample Frequency */

      //read in disposal cost rate
      else if(strstr(lineStr, "DisposalRate") != NULL)
      {
         sscanf(lineStr, "%s %lf", tmp1, &m_DisposalRate);
      }/* end else if() --> Disposal Cost Rate */

      //read in maintenance factor
      else if(strstr(lineStr, "MaintFactor") != NULL)
      {
         sscanf(lineStr, "%s %lf", tmp1, &m_MaintFactor);
      }/* end else if() --> Maintenance Factor */

      //read in remediation time frame
      else if(strstr(lineStr, "TimeFrame") != NULL)
      {
         sscanf(lineStr, "%s %lf", tmp1, &m_TimeFrame);
      }/* end else if() --> Remediation Time Frame */

      //read in interest rate
      else if(strstr(lineStr, "InterestRate") != NULL)
      {
         sscanf(lineStr, "%s %lf", tmp1, &m_IntRate);
      }/* end else if() --> Interest Rate */

      //read in treatment capital cost coeff
      else if(strstr(lineStr, "TreatCapCoeff") != NULL)
      {
         sscanf(lineStr, "%s %lf", tmp1, &m_TreatCapCoeff);
      }/* end else if() --> Treatment Capital Cost Coefficient */

      //read in treatment capital cost exponent
      else if(strstr(lineStr, "TreatCapExpon") != NULL)
      {
         sscanf(lineStr, "%s %lf", tmp1, &m_TreatCapExpon);
      }/* end else if() --> Treatment Capital Cost Exponent */

      //read in treatment operational cost coeff
      else if(strstr(lineStr, "TreatOpCoeff") != NULL)
      {
         sscanf(lineStr, "%s %lf", tmp1, &m_TreatOpCoeff);
      }/* end else if() --> Treatment Operational Cost Coefficient */

      //read in treatment operation cost exponent
      else if(strstr(lineStr, "TreatOpExpon") != NULL)
      {
         sscanf(lineStr, "%s %lf", tmp1, &m_TreatOpExpon);
      }/* end else if() --> Treatment Operational Cost Exponent */

      else
      {
         sprintf(tmp1, "PATO::InitFromFile(): unknown token |%s|", lineStr);
         LogError(ERR_FILE_IO, tmp1);
      }

      lineStr = GetNxtDataLine(pFile, fileName);
   }/* end while() */

   fclose(pFile);

   //read in plumes
   InitPlumes();

   //read in response variables
   InitResponseVars();

   //read in constraints
   InitConstraints();

   InitLookupTable();
   
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

   //assemble array of wells structures
   InitWells();

   /* Display well information 
   WriteWells(stdout, WRITE_BNR);
   WriteWells(stdout, WRITE_SCI); */
}/* end InitFromFile()*/

/******************************************************************************
PATO::InitPlumes()

Initialize all plume geometries by reading them from the input file.
******************************************************************************/
void PATO::InitPlumes(void)
{
   FILE * pFile;
   char * lineStr;
   char tmp1[DEF_STR_SZ], tmp2[DEF_STR_SZ];
   IroncladString fileName = GetInFileName();
   int state, i, len;
   double tmpx, tmpy;

   pFile = fopen(fileName, "r");

   if(pFile == NULL)
   {
      FileOpenFailure("PATO::InitPlumes", fileName);
   }/* end if() */

   if(CheckToken(pFile, "BeginPlumeGeometry", fileName) == true)
   {
      FindToken(pFile, "EndPlumeGeometry", fileName);
      rewind(pFile);

      /* ----------------------------------------------------
      Count the number of plume names and size the plume 
      array accordingly.
      -----------------------------------------------------*/
      FindToken(pFile, "BeginPlumeGeometry", fileName);
      lineStr = GetNxtDataLine(pFile, fileName);
      m_NumPlumes = 0;

      while(strcmp(lineStr, "EndPlumeGeometry") != 0)
      {    
         sscanf(lineStr, "%s %s", tmp1, tmp2);
         if(strcmp(tmp1, "PlumeName") == 0){ m_NumPlumes++;}
         lineStr = GetNxtDataLine(pFile, fileName);
      }/* end while() --> Count Number of Plumes */

      NEW_PRINT("Plume2D", m_NumPlumes);
      m_pPlumes = new Plume2D[m_NumPlumes];
      MEM_CHECK(m_pPlumes);
      for(i = 0; i < m_NumPlumes; i++)
      {
         m_pPlumes[i].name = NULL;
         m_pPlumes[i].nv = 0;
         m_pPlumes[i].poly = NULL;
      }/* end for() */

      rewind(pFile);
      FindToken(pFile, "BeginPlumeGeometry", fileName);
      lineStr = GetNxtDataLine(pFile, fileName);
      state = 0; //expecting plume name
      i = 0;

      /*--------------------------------------------------------------------------
      The while loop processes a state macine for reading in plume geometries. This
      allows multiple plumes to be embedded into the PlumeGeometry section. The 
      ResizePlume() routine is used to dynamically resize and insert polygon 
      vertices.
      --------------------------------------------------------------------------*/
      while(strcmp(lineStr, "EndPlumeGeometry") != 0)
      {         
         switch(state)
         {
            case(0) : //expecting PlumeName
            {
               strcpy(tmp1, "");
               strcpy(tmp2, "");
               sscanf(lineStr, "%s %s", tmp1, tmp2);

               if(strcmp(tmp1, "PlumeName") == 0)
               {
                  len = (int)strlen(tmp2)+1;
                  NEW_PRINT("char", len);
                  m_pPlumes[i].name = new char[len];
                  MEM_CHECK(m_pPlumes[i].name);
                  strcpy(m_pPlumes[i].name, tmp2);
                  state = 1;
               }
               else
               {
                  sprintf(tmp2, "PATO::InitPlumes() expected PlumeName, got |%s|", tmp1);
                  LogError(ERR_FILE_IO, tmp2);
               }
               break;
            }/* end case(0) --> PlumeName */
            case(1) : //expecting BeginPlumeCoords
            {
               if(strcmp(lineStr, "BeginPlumeCoords") == 0){ state = 2;}
               else
               {
                  sprintf(tmp1, "PATO::InitPlumes() expected BeginPlumeCoords, got |%s|", lineStr);
                  LogError(ERR_FILE_IO, tmp1);
               }
               break;
            }/* end case(1) --> BeginPlumeCoords */
            case(2) :
            {
               if(strcmp(lineStr, "EndPlumeCoords") == 0)
               { 
                  i++;
                  state = 0;
               }
               else //must be a coordinate pair
               {
                  sscanf(lineStr, "%lf %lf", &tmpx, &tmpy);
                  ResizePlume(tmpx, tmpy, &(m_pPlumes[i]));
               }
               break;
            }/* end case(2) --> EndPlumeCoords */
         }/* end switch(state) */

         lineStr = GetNxtDataLine(pFile, fileName);
      }/* end while() --> Reading in plume gemoetries */
   } /* end if() */
   else
   {
      LogError(ERR_FILE_IO, "No plume geometry specified.");
   }/* end else() */

   fclose(pFile);
}/* end InitPlumes()*/

/******************************************************************************
PATO::ResizePlume()

Add a coordinate pair to the plume polygon.
******************************************************************************/
void PATO::ResizePlume(double x, double y, Plume2D * pPlume)
{
   Point2D * pTmp;
   int nv, i;

   if(pPlume->nv == 0)
   {
      pPlume->nv = 1;
      NEW_PRINT("Point2D", 1);
      pPlume->poly = new Point2D[1];
      MEM_CHECK(pPlume->poly);
      pPlume->poly[0].x = x;
      pPlume->poly[0].y = y;
   }
   else
   {
      nv = (pPlume->nv + 1);
      NEW_PRINT("Point2D", nv);
      pTmp = new Point2D[nv];
      MEM_CHECK(pTmp);
      //copy previous coords
      for(i = 0; i < (nv-1); i++)
      {
         pTmp[i].x = pPlume->poly[i].x;
         pTmp[i].y = pPlume->poly[i].y;
      }/* end for() */
      //add coord
      pTmp[i].x = x;
      pTmp[i].y = y;
      //replace plume ploygon
      delete [] (pPlume->poly);
      pPlume->poly = pTmp;
      pPlume->nv = nv;
   }/* end else() */
}/* end ResizePlume() */

/******************************************************************************
PATO::InitResponseVars()

Initialize all response variables, which are the basis for the constraints.
******************************************************************************/
void PATO::InitResponseVars(void)
{
   NEW_PRINT("ResponseVarGroup", 1);
   m_pRespGroup = new ResponseVarGroup();
   MEM_CHECK(m_pRespGroup);
}/* end InitResponseVars()*/

/******************************************************************************
PATO::InitConstraints()

Initialize all constraints by parsing the information in the "Constraints" 
section of the input file.
******************************************************************************/
void PATO::InitConstraints(void)
{
   FILE * pFile;
   StringType * names;
   Point2D * pPlume;
   RespVarABC * pLoc1, * pLoc2;
   char * pTok, * pOld;
   int i, j, np, nv, len;
   double conv, lwr, upr;
   char * lineStr, nameStr[DEF_STR_SZ], typeStr[DEF_STR_SZ];
   char tmp1[DEF_STR_SZ], tmp2[DEF_STR_SZ], tmp3[DEF_STR_SZ];
   IroncladString fileName = GetInFileName();

   pFile = fopen(fileName, "r");

   if(pFile == NULL)
   {
      FileOpenFailure("PATO::InitConstraints", fileName);
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

         if(strcmp(typeStr, "capacity") == 0)
         {
            pTok = lineStr;
            //extract name of constraint (no spaces allowed)
            j = ExtractString(pTok, nameStr);
            j = ValidateExtraction(j, 1, 1, "PATO::InitConstraints()");
            pTok += j;
            //extract type  
            j = ExtractString(pTok, typeStr);
            j = ValidateExtraction(j, 1, 1, "PATO::InitConstraints()");
            pTok += j;
            //extract conversion factor
            j = ExtractString(pTok, tmp1);
            j = ValidateExtraction(j, 1, 1, "PATO::InitConstraints()");
            conv = atof(tmp1);
            pTok += j;
            //extract lower bound
            j = ExtractString(pTok, tmp1);
            j = ValidateExtraction(j, 1, 1, "PATO::InitConstraints()");
            lwr = atof(tmp1);
            pTok += j;
            //extract upper bound
            j = ExtractString(pTok, tmp1);
            j = ValidateExtraction(j, 1, 1, "PATO::InitConstraints()");
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
               j = ValidateExtraction(j, i, np, "PATO::InitConstraints()");
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
         else if (strcmp(typeStr, "drawdown") == 0)
         {
            sscanf(lineStr, "%s %s %lf %lf %lf %s", nameStr, typeStr, &conv, &lwr, &upr, tmp1);
            pLoc1 = m_pRespGroup->GetRespVarPtr(tmp1);
            if(pLoc1 == NULL)
            {
               sprintf(typeStr, "PATO::InitConstraints() unknown response variable |%s|", tmp1);
               LogError(ERR_FILE_IO, typeStr);
               ExitProgram(1);
            }/* end if() */
            //create constraint
            NEW_PRINT("DrawdownConstraint", 1);
            DrawdownConstraint * pNewDD = new DrawdownConstraint(nameStr, pLoc1, lwr, upr, conv);
            MEM_CHECK(pNewDD);
            //insert into linked list
            if(m_pConstraints == NULL){ m_pConstraints = pNewDD;}
            else{ m_pConstraints->AddConstraint(pNewDD);}
         }/* end else if() ---> Drawdown Constraint */
         else if (strcmp(typeStr, "general") == 0)
         {
            sscanf(lineStr, "%s %s %lf %lf %lf %s", nameStr, typeStr, &conv, &lwr, &upr, tmp1);
            pLoc1 = m_pRespGroup->GetRespVarPtr(tmp1);
            if(pLoc1 == NULL)
            {
               sprintf(typeStr, "PATO::InitConstraints() unknown response variable |%s|", tmp1);
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
         else if (strcmp(typeStr, "hydgrad") == 0)
         {            
            sscanf(lineStr, "%s %s %lf %lf %lf %s %s", nameStr, typeStr, &conv, &lwr, &upr, tmp1, tmp2);
            pLoc1 = m_pRespGroup->GetRespVarPtr(tmp1);
            pLoc2 = m_pRespGroup->GetRespVarPtr(tmp2);
            if(pLoc1 == NULL)
            {
               sprintf(typeStr, "PATO::InitConstraints() unknown response variable |%s|", tmp1);
               LogError(ERR_FILE_IO, typeStr);
               ExitProgram(1);
            }/* end if() */
            if(pLoc2 == NULL)
            {
               sprintf(typeStr, "PATO::InitConstraints() unknown response variable |%s|", tmp2);
               LogError(ERR_FILE_IO, typeStr);
               ExitProgram(1);
            }/* end if() */

            //create constraint
            NEW_PRINT("HydGradConstraint", 1);
            HydGradConstraint * pNewHG = new HydGradConstraint(nameStr, pLoc1, pLoc2, lwr, upr, conv);
            MEM_CHECK(pNewHG);
            //insert into linked list
            if(m_pConstraints == NULL){ m_pConstraints = pNewHG;}
            else{ m_pConstraints->AddConstraint(pNewHG);}
         }/* end else if() ---> Hydraulic Gradient Constraint */
         else if (strcmp(typeStr, "partcap") == 0)
         {            
            sscanf(lineStr, "%s %s %lf %s %s %s", nameStr, typeStr, &conv, tmp1, tmp2, tmp3);
            pLoc1 = m_pRespGroup->GetRespVarPtr(tmp1);
            pLoc2 = m_pRespGroup->GetRespVarPtr(tmp2);
            if(pLoc1 == NULL)
            {
               sprintf(typeStr, "PATO::InitConstraints() unknown response variable |%s|", tmp1);
               LogError(ERR_FILE_IO, typeStr);
               ExitProgram(1);
            }/* end if() */
            if(pLoc2 == NULL)
            {
               sprintf(typeStr, "PATO::InitConstraints() unknown response variable |%s|", tmp2);
               LogError(ERR_FILE_IO, typeStr);
               ExitProgram(1);
            }/* end if() */

            pPlume = NULL;
            for(i = 0; i < m_NumPlumes; i++)
            {
               if(strcmp(m_pPlumes[i].name, tmp3) == 0)
               {
                  pPlume = m_pPlumes[i].poly;
                  nv = m_pPlumes[i].nv;                  
               }
            }/* end for() */

            if(pPlume == NULL)
            {
               sprintf(typeStr, "PATO::InitConstraints() unknown plume name |%s|", tmp3);
               LogError(ERR_FILE_IO, typeStr);
               ExitProgram(1);
            }/* end if() */

            //create constraint
            NEW_PRINT("ParticleCaptureConstraint", 1);
            ParticleCaptureConstraint * pNewPC = new ParticleCaptureConstraint(nameStr, pLoc1, pLoc2, pPlume, nv, conv);
            MEM_CHECK(pNewPC);
            //insert into linked list
            if(m_pConstraints == NULL){ m_pConstraints = pNewPC;}
            else{ m_pConstraints->AddConstraint(pNewPC);}
         }/* end else if() ---> Particle Capture Constraint */
         else
         {
            sprintf(nameStr, "PATO::InitConstrints() unknown type |%s|", typeStr);
            LogError(ERR_FILE_IO, nameStr);
         }         
         lineStr = GetNxtDataLine(pFile, fileName);
      }/* end while() --> Count Number of Plumes */
   } /* end if() */
   else
   {
      LogError(ERR_FILE_IO, "No constraints specified.");
   }/* end else() */

   fclose(pFile);
}/* end InitConstraints()*/

/******************************************************************************
PATO::InitWells()

Initialize the list of candidate wells. This is a listing of m_MaxNumWells
and contains pointers to the pumping rate, x-location and y-location parameters
(stored in m_pParamGroup), the response variables used in computation of lift
(head and surface elevation at well), and the name given to the well.
******************************************************************************/
void PATO::InitWells(void)
{
   FILE * pFile;
   char * lineStr;
   char name[DEF_STR_SZ], xloc[DEF_STR_SZ], yloc[DEF_STR_SZ], rate[DEF_STR_SZ];
   char head[DEF_STR_SZ], topo[DEF_STR_SZ], base[DEF_STR_SZ], msg[DEF_STR_SZ];

   IroncladString fileName = GetInFileName();
   int i, len;

   pFile = fopen(fileName, "r");

   if(pFile == NULL)
   {
      FileOpenFailure("PATO::InitWells", fileName);
   }/* end if() */

   if(CheckToken(pFile, "BeginCandidateWells", fileName) == true)
   {
      FindToken(pFile, "EndCandidateWells", fileName);
      rewind(pFile);

      /* ----------------------------------------------------
      Count the number of wells and size the well array accordingly.
      -----------------------------------------------------*/
      FindToken(pFile, "BeginCandidateWells", fileName);
      lineStr = GetNxtDataLine(pFile, fileName);
      m_MaxNumWells = 0;

      while(strcmp(lineStr, "EndCandidateWells") != 0)
      {    
         m_MaxNumWells++;
         lineStr = GetNxtDataLine(pFile, fileName);
      }/* end while() --> Count Number of Wells */

      NEW_PRINT("WellStruct", m_MaxNumWells);
      m_pWells = new WellStruct[m_MaxNumWells];
      MEM_CHECK(m_pWells);
      for(i = 0; i < m_MaxNumWells; i++)
      {
         m_pWells[i].name   = NULL;
         m_pWells[i].pQ     = NULL;
         m_pWells[i].pXloc  = NULL;
         m_pWells[i].pYloc  = NULL;
         m_pWells[i].pHead  = NULL;
         m_pWells[i].pTopo  = NULL;
         m_pWells[i].Topo   = 0.00;
         m_pWells[i].pBase  = NULL;
         m_pWells[i].Base   = 0.00;
         m_pWells[i].Cdrill = 0.00;
         m_pWells[i].Cpump  = 0.00;
         m_pWells[i].Cnrg   = 0.00;
         m_pWells[i].Ctot   = 0.00;
      }/* end for() */

      rewind(pFile);
      FindToken(pFile, "BeginCandidateWells", fileName);
      lineStr = GetNxtDataLine(pFile, fileName);
      i = 0;

      //read in well information
      while(strcmp(lineStr, "EndCandidateWells") != 0)
      {
         //init. tmp vars, in case sscanf() fails....
         name[0] = 0; xloc[0] = 0; yloc[0] = 0; 
         rate[0] = 0; head[0] = 0; topo[0] = 0; base[0] = 0;

         sscanf(lineStr, "%s %s %s %s %s %s %s", name, xloc, yloc, rate, head, topo, base);

         //name
         len = (int)strlen(name)+1;
         NEW_PRINT("char", len);
         m_pWells[i].name = new char[len];
         MEM_CHECK(m_pWells[i].name);
         strcpy(m_pWells[i].name, name);

         //x-coord
         m_pWells[i].pXloc = m_pParamGroup->GetParamPtr(xloc);
         //if name not found in list, must abort program
         if( m_pWells[i].pXloc == NULL)
         {
            sprintf(msg, "PATO::InitWells(), unknown parameter : |%s|", xloc);
            LogError(ERR_FILE_IO, msg);
            ExitProgram(1);         
         }/* end if() */

         //y-coord
         m_pWells[i].pYloc = m_pParamGroup->GetParamPtr(yloc);
         //if name not found in list, must abort program
         if(m_pWells[i].pYloc == NULL)
         {
            sprintf(msg, "PATO::InitWells(), unknown parameter : |%s|", yloc);
            LogError(ERR_FILE_IO, msg);
            ExitProgram(1);         
         }/* end if() */

         //rate
         m_pWells[i].pQ = m_pParamGroup->GetParamPtr(rate);
         //if name not found in list, must abort program
         if(m_pWells[i].pQ == NULL)
         {
            sprintf(msg, "PATO::InitWells(), unknown parameter : |%s|", rate);
            LogError(ERR_FILE_IO, msg);
            ExitProgram(1);         
         }/* end if() */
         m_pWells[i].pQ->SetThreshVal(-m_RateThresh, m_RateThresh, 0.00);
         
         if(m_ObjType != PATO_OBJ_RATE)
         {
            //head at well (water table elevation)
            m_pWells[i].pHead = m_pRespGroup->GetRespVarPtr(head);
            if(m_pWells[i].pHead == NULL)
            {
               sprintf(msg, "PATO::InitWells() unknown response variable |%s|", head);
               LogError(ERR_FILE_IO, msg);
               ExitProgram(1);
            }/* end if() */

            //topography at well (surface elevation)
            m_pWells[i].pTopo = m_pRespGroup->GetRespVarPtr(topo);
            if(m_pWells[i].pTopo == NULL)
            {
               //topo string might refer to a constant....
               for(int c = 0; c < (int)strlen(topo); c++)
               {
                  if(((topo[c] < '0') || (topo[c] > '9')) && (topo[c] != '.') && 
                     (topo[c] != 'e') && (topo[c] != 'E') && (topo[c] != '+') &&
                     (topo[c] != '-'))
                  {
                     //not a valid number format
                     sprintf(msg, "PATO::InitWells() unknown response variable or invalid number format |%s|", topo);
                     LogError(ERR_FILE_IO, msg);
                     ExitProgram(1);
                  }
               }
               m_pWells[i].Topo = atof(topo);               
            }/* end if() */

            //aquifer base at well
            m_pWells[i].pBase = m_pRespGroup->GetRespVarPtr(base);
            if(m_pWells[i].pBase == NULL)
            {
               //base string might refer to a constant....
               for(int c = 0; c < (int)strlen(base); c++)
               {
                  if(((base[c] < '0') || (base[c] > '9')) && (base[c] != '.') && 
                     (base[c] != 'e') && (base[c] != 'E') && (base[c] != '+') &&
                     (base[c] != '-'))
                  {
                     //not a valid number format
                     sprintf(msg, "PATO::InitWells() unknown response variable or invalid number format |%s|", base);
                     LogError(ERR_FILE_IO, msg);
                     ExitProgram(1);
                  }
               }
               m_pWells[i].Base = atof(base);
            }/* end if() */
         }/* end if() */
         
         lineStr = GetNxtDataLine(pFile, fileName);
         i++;
      }/* end while() --> Reading in wells */
   } /* end if() */
   else
   {
      LogError(ERR_FILE_IO, "No wells specified.");
   }/* end else() */

   fclose(pFile);  
}/* end InitWells()*/

/******************************************************************************
PATO::InitLookupTable()

Initialize the pump cost lookup table. Either the data is read from a file, or 
it is initialized to a default lookup table transcribed from the RS Means 
catalog.
******************************************************************************/
void PATO::InitLookupTable(void)
{
   FILE * pFile;
   char * pTok;
   char * lineStr, tmpStr[DEF_STR_SZ];

   IroncladString fileName = GetInFileName();
   int i, j;

   pFile = fopen(fileName, "r");

   if(pFile == NULL)
   {
      FileOpenFailure("PATO::InitLookupTable", fileName);
   }/* end if() */

   if(CheckToken(pFile, "BeginLookupTable", fileName) == true)
   {
      FindToken(pFile, "EndLookupTable", fileName);
      rewind(pFile);

      /* ----------------------------------------------------
      Count the number of table entries and size the table 
      array accordingly.
      -----------------------------------------------------*/
      FindToken(pFile, "BeginLookupTable", fileName);
      lineStr = GetNxtDataLine(pFile, fileName);
      m_TblSize = 0;

      while(strcmp(lineStr, "EndLookupTable") != 0)
      {    
         m_TblSize++;
         lineStr = GetNxtDataLine(pFile, fileName);
      }/* end while() --> Count Number of Wells */

      NEW_PRINT("PumpLkupTableStruct", m_TblSize);
      m_pTbl = new PumpLkupTableStruct[m_TblSize];
      MEM_CHECK(m_pTbl);
      for(i = 0; i < m_TblSize; i++)
      {
         m_pTbl[i].cost  = 0.00;
         m_pTbl[i].Lmax  = 0.00;
         m_pTbl[i].Lmin  = 0.00;
         m_pTbl[i].Qmax  = 0.00;
         m_pTbl[i].Qmin  = 0.00;
      }/* end for() */

      rewind(pFile);
      FindToken(pFile, "BeginLookupTable", fileName);
      lineStr = GetNxtDataLine(pFile, fileName);
      i = 0;

      //read in table information
      while(strcmp(lineStr, "EndLookupTable") != 0)
      {
         pTok = lineStr;
         //extract Qmin
         j = ExtractString(pTok, tmpStr);
         j = ValidateExtraction(j, 1, 1, "PATO::InitLookupTable()");
         pTok += j;
         m_pTbl[i].Qmin = atof(tmpStr);
         //extract Qmax
         j = ExtractString(pTok, tmpStr);
         j = ValidateExtraction(j, 1, 1, "PATO::InitLookupTable()");
         pTok += j;
         m_pTbl[i].Qmax = atof(tmpStr);
         //extract Lmin
         j = ExtractString(pTok, tmpStr);
         j = ValidateExtraction(j, 1, 1, "PATO::InitLookupTable()");
         pTok += j;
         m_pTbl[i].Lmin = atof(tmpStr);
         //extract Lmax
         j = ExtractString(pTok, tmpStr);
         j = ValidateExtraction(j, 1, 1, "PATO::InitLookupTable()");
         pTok += j;
         m_pTbl[i].Lmax = atof(tmpStr);
         //extract cost
         j = ExtractString(pTok, tmpStr);         
         pTok += j;
         m_pTbl[i].cost = atof(tmpStr);
         
         lineStr = GetNxtDataLine(pFile, fileName);
         i++;
      }/* end while() --> Reading in wells */
   } /* end if() */
   else
   {
      LogError(ERR_FILE_IO, "Using default lookup table for pump costs.");
      m_TblSize = 77;
      NEW_PRINT("PumpLkupTableStruct", m_TblSize);
      m_pTbl = new PumpLkupTableStruct[m_TblSize];
      MEM_CHECK(m_pTbl);
      /*--------------------------------------------------------------------
      Fill in values in GPM, feet and US dollars read from the RS Means 2004 
      reference.
      ---------------------------------------------------------------------*/
      for(i = 0; i < m_TblSize; i++)
      {         
         m_pTbl[i].Qmin  = gPumpCosts[i][0];
         m_pTbl[i].Qmax  = gPumpCosts[i][1];
         m_pTbl[i].Lmin  = gPumpCosts[i][2];
         m_pTbl[i].Lmax  = gPumpCosts[i][3];
         m_pTbl[i].cost  = gPumpCosts[i][4];
      }/* end for() */
      //convert values to m3/day and meters
      for(i = 0; i < m_TblSize; i++)
      {
         m_pTbl[i].Lmax  *= 0.3048;
         m_pTbl[i].Lmin  *= 0.3048;
         m_pTbl[i].Qmax  *= 5.4496;
         m_pTbl[i].Qmin  *= 5.4496;
      }/* end for() */
   }/* end else() */

   fclose(pFile);  
}/* end InitLookupTable() */

/******************************************************************************
PATO::Destroy()
******************************************************************************/
void PATO::Destroy(void)
{
   int i;

   delete m_pConstraints;

   for(i = 0; i < m_MaxNumWells; i++)
   {
      delete [] m_pWells[i].name;
   }
   delete [] m_pWells;

   delete [] m_pTbl;

  //destroy plumes
   for(i = 0; i < m_NumPlumes; i++)
   {
      delete [] m_pPlumes[i].name;
      delete [] m_pPlumes[i].poly;
   }
   delete [] m_pPlumes;

   delete m_pRespGroup;

   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
PATO::CalcObjFunc()

Computes the objective function and returns the result.
******************************************************************************/
double PATO::CalcObjFunc(void)
{
   static bool firstTime = true;
   FILE * pFile;
   double cost = 0.00;
   double penalty = 0.00;
   int i;
  
   ConstraintABC * pCur;
   
   m_pRespGroup->ExtractVals();
   
   if(firstTime == true)
   {
      pFile = fopen("OstPatoOut.txt", "w");
      fprintf(pFile, "True Cost \tPenalty \tAdjusted Cost\n");
      firstTime = false;
   }
   else
   { 
      pFile = fopen("OstPatoOut.txt", "a+");
   }

   //reset total cost
   for(i = 0; i < m_MaxNumWells; i++){ m_pWells[i].Ctot = 0;}

   switch(m_ObjType)
   {
      case(PATO_OBJ_RATE) :
      {
         cost = CalcPumpingRate();
         break;
      }
      case(PATO_OBJ_OP) :
      {
         cost = CalcOperationCost();
         break;
      } 
      case(PATO_OBJ_CAP_OP) :
      {
         cost  = CalcCapitalCost();
         cost += CalcOperationCost();
         break;
      } 
      case(PATO_OBJ_MAYER) :
      {
         cost  = CalcCapitalCost();
         cost += CalcOperationCost();
         break;
      } 
      case(PATO_OBJ_CAP_OP_TRE) :
      {
         cost  = CalcCapitalCost();
         cost += CalcOperationCost();
         cost += CalcTreatmentCost();
         break;
      } 
   }/* end switch() */

   //compute total cost for each well
   for(i = 0; i < m_MaxNumWells; i++)
   {
      m_pWells[i].Ctot += m_pWells[i].Cpump;
      m_pWells[i].Ctot += m_pWells[i].Cdrill;
      m_pWells[i].Ctot += m_pWells[i].Cnrg;
   }/* end for() */
   

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
} /* end PATO::CalcObjFunc() */

/******************************************************************************
PATO::CalcPumpingRate()

Computes the total pumping rate as the sum of the recharge and discharge rates.
******************************************************************************/
double PATO::CalcPumpingRate(void)
{
   int i;
   double pump_sum, inj_sum, sum, Qi;
   
   pump_sum = inj_sum = 0.00;
   for(i = 0; i < m_MaxNumWells; i++)
   {
      Qi = m_pWells[i].pQ->GetEstVal();
      //pumping wells
      if(Qi > m_RateThresh)
      {
         pump_sum += Qi;
         m_pWells[i].Ctot = Qi*m_ExtRateCF;
      }
      else if(Qi < -m_RateThresh) //injection wells
      {
         inj_sum += fabs(Qi);
         m_pWells[i].Ctot = fabs(Qi)*m_InjRateCF;
      }
      else{} //inactive well
   }/* end for() */
   m_Costs[EXT_COST] = pump_sum*m_ExtRateCF;
   m_Costs[INJ_COST] = inj_sum*m_InjRateCF;

   sum = pump_sum*m_ExtRateCF + inj_sum*m_InjRateCF;
   return (sum);
}/* end CalcPumpingRate() */

/******************************************************************************
PATO::CalcOperationCost()

Computes the total operational cost as a function of pumping rate and hydraulic 
lift.
******************************************************************************/
double PATO::CalcOperationCost(void)
{
   int i, max;
   double sum1, sum2, T, NW;
   double Qk, Hk, Tk; //rate, head and surface elevation at well k
   double CL, CE, CA, CD, CM, cap;
   double InjTot, ExtTot;
         
   InjTot = ExtTot = NW = CL = CA = CE = CD = CM = 0.00;
   
   if(m_IntRate == 0.00){ T = m_TimeFrame;}
   else{ T = (1.00 - pow((1 - m_IntRate), -m_TimeFrame)) / m_IntRate;}

   //energy cost due to injection (depends on injection rate)
   sum1 = 0.00;
   for(i = 0; i < m_MaxNumWells; i++)
   {
      Qk = m_pWells[i].pQ->GetEstVal();
      //injection wells
      if(Qk < -m_RateThresh)
      { 
         sum1 += fabs(Qk); 
         m_pWells[i].Cnrg = fabs(Qk)*m_InjEnergyRate*T;
      }
      
      //count total number of active wells
      if((Qk > m_RateThresh) || (Qk < -m_RateThresh))
      {  
         NW = NW + 1.00;
      }
   }/* end for() */
   InjTot = sum1;
   sum1 *= m_InjEnergyRate;

   //energy cost due to extraction (depends on lift)
   sum2 = 0.00;   
   for(i = 0; i < m_MaxNumWells; i++)
   {
      Qk = m_pWells[i].pQ->GetEstVal();
      //only compute lift for pumping wells
      if(Qk > m_RateThresh)
      {
         ExtTot += Qk;
         Hk = m_pWells[i].pHead->GetCurrentVal();
         if(m_pWells[i].pTopo != NULL){ Tk = m_pWells[i].pTopo->GetCurrentVal();}
         else{ Tk = m_pWells[i].Topo;}
         if(Tk > Hk)
         { 
            sum2 += ((Tk-Hk)*Qk);
            m_pWells[i].Cnrg = (Tk-Hk)*Qk*m_ExtEnergyRate*T;
         }
      }/* end if() */
   }/* end for() */
   sum2 *= m_ExtEnergyRate;

   CE = sum1+sum2;

   //labor and analytic costs
   CL = m_LaborRate*110.00*sqrt(NW/3.00);
   CA = m_AnalyticRate*m_SampleFreq*10.00*NW;
   //disposal costs
   CD = m_DisposalRate*(ExtTot - InjTot);
   if(CD < 0.00){ CD = 0.00;}
   //maintenance total
   if(m_ObjType != PATO_OBJ_OP)
   {
      if(m_TimeFrame < 30.00){ max = (int)m_TimeFrame;}
      else{max = 30;}
   
      for(i = 0; i < max; i++)
      {
         CM += gAdjFact[i];
      }/* end for() */
      cap = CalcCapitalCost();
      CM *= (m_MaintFactor*cap);
   }

   m_Costs[LAB_COST] = CL*T;
   m_Costs[NRG_COST] = CE*T;
   m_Costs[ANA_COST] = CA*T;
   m_Costs[DIS_COST] = CD*T;
   m_Costs[MNT_COST] = CM*T;

   return ((CL+CE+CA+CD+CM)*T);
}/* end CalcOperationCost() */

/******************************************************************************
PATO::CalcCapitalCost()

Computes the total capital cost as a function of the number of wells, surface 
elevation and pumping rate.
******************************************************************************/
double PATO::CalcCapitalCost(void)
{
   if(m_ObjType == PATO_OBJ_MAYER) return CalcMayerCost();

   /* RS Means calculation of capital costs */
   int i;
   double sum1, sum2, Qi, di, Bi, Ti, Hi, sum_fix;
   
   sum_fix = sum1 = sum2 = 0.00;

   //construction costs for active wells (depends on depth to well)
   for(i = 0; i < m_MaxNumWells; i++)
   {
      Qi = m_pWells[i].pQ->GetEstVal();
      if((Qi > m_RateThresh) || (Qi < -m_RateThresh))
      { 
         /* Depth of the well is from surface to aquifer base */
         if(m_pWells[i].pBase != NULL){ Bi = m_pWells[i].pBase->GetCurrentVal();}
         else{ Bi = m_pWells[i].Base;}

         if(m_pWells[i].pTopo != NULL){ Ti = m_pWells[i].pTopo->GetCurrentVal();}
         else{ Ti = m_pWells[i].Topo;}

         if(Ti > Bi){ di = Ti-Bi;}
         else{ di = 0.00;}

         sum1 += di;
         sum_fix += m_FixWellCF;
         m_pWells[i].Cdrill = (di * m_VarWellCF)+m_FixWellCF;
      }
   }/* end for() */
   sum1 *= m_VarWellCF;

   //pump costs (depends on rate and depth to well)
   for(i = 0; i < m_MaxNumWells; i++)
   {
      Qi = m_pWells[i].pQ->GetEstVal();
      if(Qi > m_RateThresh) //extraction wells only
      { 
         Hi = m_pWells[i].pHead->GetCurrentVal();
         if(m_pWells[i].pTopo != NULL){ Ti = m_pWells[i].pTopo->GetCurrentVal();}
         else{ Ti = m_pWells[i].Topo;}
         if(Ti > Hi){ di = Ti-Hi;}
         else{ di = 0.00;}
               
         m_pWells[i].Cpump = LookupPumpCost((m_RateUCF*Qi), (m_LiftUCF*di));
         sum2 += m_pWells[i].Cpump;
      }
   }/* end for() */
   
   m_Costs[DRL_COST] = sum1+sum_fix;
   m_Costs[PMP_COST] = sum2;

   return (sum1 + sum2 + sum_fix);
}/* end CalcCapitalCost() */

/******************************************************************************
PATO::CalcMayerCost()

Computes the total capital cost as a function of the number of wells, surface 
elevation and pumping rate, based on Mayer's community problem formulation.
******************************************************************************/
double PATO::CalcMayerCost(void)
{
   int i;
   double sum1, sum2, Qi;
   
   sum1 = sum2 = 0.00;

   //drlling costs for active wells (depends on depth to well)
   for(i = 0; i < m_MaxNumWells; i++)
   {
      Qi = m_pWells[i].pQ->GetEstVal();
      if((Qi > m_RateThresh) || (Qi < -m_RateThresh))
      { 
         sum1 += m_MayerDrillCF;
         m_pWells[i].Cdrill = m_MayerDrillCF;
      }
   }/* end for() */

   //pump costs (depends on number of pumping wells)
   for(i = 0; i < m_MaxNumWells; i++)
   {
      Qi = m_pWells[i].pQ->GetEstVal();
      if(Qi > m_RateThresh) //extraction wells only
      {             
         m_pWells[i].Cpump = m_MayerPumpCF;
         sum2 += m_MayerPumpCF;
      }
   }/* end for() */
   
   m_Costs[DRL_COST] = sum1;
   m_Costs[PMP_COST] = sum2;

   return (sum1 + sum2);
}/* end CalcMayerCost() */

/******************************************************************************
PATO::CalcTreatmentCost()

Computes the total treament cost as a function of the total pumping rate.
******************************************************************************/
double PATO::CalcTreatmentCost(void)
{
   int i;
   double Qi, Qtot, sum1, sum2, T;
   
   Qtot = sum1 = sum2 = 0.00;

   //compute total extraction rate
   for(i = 0; i < m_MaxNumWells; i++)
   {
      Qi = m_pWells[i].pQ->GetEstVal();
      if(Qi > m_RateThresh){ Qtot += Qi;}
   }/* end for() */

   //compute capital costs
   sum1 = m_TreatCapCoeff * pow(Qtot, m_TreatCapExpon);

   //compute operational costs
   if(m_IntRate == 0.00){ T = m_TimeFrame;}
   else{ T = (1.00 - pow((1 - m_IntRate), -m_TimeFrame)) / m_IntRate;}

   sum2 = (m_TreatOpCoeff * pow(Qtot, m_TreatOpExpon))*T;

   m_Costs[TRC_COST] = sum1;
   m_Costs[TRO_COST] = sum2;

   return (sum1 + sum2);
}/* end CalcTreatmentCost() */

/******************************************************************************
PATO::WriteCost()

Display cost information.
******************************************************************************/
void PATO::WriteCost(FILE * pFile, int type)
{
   double Cext, Cinj, Clab, Cnrg, Cana, Cdis, Cmnt, Cdrl, Cpmp, Ctrc, Ctro;

   Cext = m_Costs[EXT_COST];
   Cinj = m_Costs[INJ_COST];
   Clab = m_Costs[LAB_COST];
   Cnrg = m_Costs[NRG_COST];
   Cana = m_Costs[ANA_COST];
   Cdis = m_Costs[DIS_COST];
   Cmnt = m_Costs[MNT_COST];
   Cdrl = m_Costs[DRL_COST];
   Cpmp = m_Costs[PMP_COST];
   Ctrc = m_Costs[TRC_COST];
   Ctro = m_Costs[TRO_COST];

   switch(type)
   {
      case(WRITE_SCI) : 
      {
         switch(m_ObjType)
         {
            case(PATO_OBJ_RATE) :
            {
               fprintf(pFile, "Extraction Cost : %E\n", Cext);
               fprintf(pFile, "Injection Cost  : %E\n", Cinj);
               fprintf(pFile, "Total Cost      : %E\n", Cext+Cinj);
               break;
            }/* end case() */
            case(PATO_OBJ_OP) :
            {
               fprintf(pFile, "Labor Cost       : %E\n", Clab);
               fprintf(pFile, "Energy Cost      : %E\n", Cnrg);
               fprintf(pFile, "Analytic Cost    : %E\n", Cana);
               fprintf(pFile, "Disposal Cost    : %E\n", Cdis);
               fprintf(pFile, "Maintenance Cost : %E\n", Cmnt);
               fprintf(pFile, "Total Cost       : %E\n", Clab+Cnrg+Cana+Cdis+Cmnt);
               break;
            }/* end case() */
            case(PATO_OBJ_CAP_OP) :
            case(PATO_OBJ_MAYER)  :
            {
               fprintf(pFile, "Labor Cost       : %E\n", Clab);
               fprintf(pFile, "Energy Cost      : %E\n", Cnrg);
               fprintf(pFile, "Analytic Cost    : %E\n", Cana);
               fprintf(pFile, "Disposal Cost    : %E\n", Cdis);
               fprintf(pFile, "Maintenance Cost : %E\n", Cmnt);
               fprintf(pFile, "Drilling Cost    : %E\n", Cdrl);
               fprintf(pFile, "Cost of Pumps    : %E\n", Cpmp);
               fprintf(pFile, "Total Cost       : %E\n", Clab+Cnrg+Cana+Cdis+Cmnt+Cdrl+Cpmp);
               break;
            }/* end case() */
            case(PATO_OBJ_CAP_OP_TRE) :
            {
               fprintf(pFile, "Labor Cost                 : %E\n", Clab);
               fprintf(pFile, "Energy Cost                : %E\n", Cnrg);
               fprintf(pFile, "Analytic Cost              : %E\n", Cana);
               fprintf(pFile, "Disposal Cost              : %E\n", Cdis);
               fprintf(pFile, "Maintenance Cost           : %E\n", Cmnt);
               fprintf(pFile, "Drilling Cost              : %E\n", Cdrl);
               fprintf(pFile, "Cost of Pumps              : %E\n", Cpmp);
               fprintf(pFile, "Treatment Capital Cost     : %E\n", Ctrc);
               fprintf(pFile, "Treatment Operational Cost : %E\n", Ctro);
               fprintf(pFile, "Total Cost       : %E\n", Clab+Cnrg+Cana+Cdis+Cmnt+Cdrl+Cpmp);
               break;
            }/* end case() */
         }/* end switch() */
         break;
      }/* end case(WRITE_SCI) */
      case(WRITE_BNR) :
      {
         break;
      }/* end case(WRITE_BNR) */
      default:
      case(WRITE_DBG) :
      case(WRITE_DEC) :
      {
         switch(m_ObjType)
         {
            case(PATO_OBJ_RATE) :
            {
               fprintf(pFile, "Extraction Cost : %10.2lf\n", Cext);
               fprintf(pFile, "Injection Cost  : %10.2lf\n", Cinj);
               fprintf(pFile, "Total Cost      : %10.2lf\n", Cext+Cinj);
               break;
            }/* end case() */
            case(PATO_OBJ_OP) :
            {
               fprintf(pFile, "Labor Cost       : %10.2lf\n", Clab);
               fprintf(pFile, "Energy Cost      : %10.2lf\n", Cnrg);
               fprintf(pFile, "Analytic Cost    : %10.2lf\n", Cana);
               fprintf(pFile, "Disposal Cost    : %10.2lf\n", Cdis);
               fprintf(pFile, "Maintenance Cost : %10.2lf\n", Cmnt);
               fprintf(pFile, "Total Cost       : %10.2lf\n", Clab+Cnrg+Cana+Cdis+Cmnt);
               break;
            }/* end case() */
            case(PATO_OBJ_CAP_OP) :
            case(PATO_OBJ_MAYER)  :
            {
               fprintf(pFile, "Labor Cost       : %10.2lf\n", Clab);
               fprintf(pFile, "Energy Cost      : %10.2lf\n", Cnrg);
               fprintf(pFile, "Analytic Cost    : %10.2lf\n", Cana);
               fprintf(pFile, "Disposal Cost    : %10.2lf\n", Cdis);
               fprintf(pFile, "Maintenance Cost : %10.2lf\n", Cmnt);
               fprintf(pFile, "Drilling Cost    : %10.2lf\n", Cdrl);
               fprintf(pFile, "Cost of Pumps    : %10.2lf\n", Cpmp);
               fprintf(pFile, "Total Cost       : %10.2lf\n", Clab+Cnrg+Cana+Cdis+Cmnt+Cdrl+Cpmp);
               break;
            }/* end case() */
            case(PATO_OBJ_CAP_OP_TRE) :
            {
               fprintf(pFile, "Labor Cost                 : %10.2lf\n", Clab);
               fprintf(pFile, "Energy Cost                : %10.2lf\n", Cnrg);
               fprintf(pFile, "Analytic Cost              : %10.2lf\n", Cana);
               fprintf(pFile, "Disposal Cost              : %10.2lf\n", Cdis);
               fprintf(pFile, "Maintenance Cost           : %10.2lf\n", Cmnt);
               fprintf(pFile, "Drilling Cost              : %10.2lf\n", Cdrl);
               fprintf(pFile, "Cost of Pumps              : %10.2lf\n", Cpmp);
               fprintf(pFile, "Treatment Capital Cost     : %10.2lf\n", Ctrc);
               fprintf(pFile, "Treatment Operational Cost : %10.2lf\n", Ctro);
               fprintf(pFile, "Total Cost       : %10.2lf\n", Clab+Cnrg+Cana+Cdis+Cmnt+Cdrl+Cpmp);
               break;
            }/* end case() */
         }/* end switch() */
         break;
      }/* end case() */
   }/* end switch() */
}/* end WriteCost() */

/******************************************************************************
PATO::WriteConstraints()

Display constraint information.
******************************************************************************/
void PATO::WriteConstraints(FILE * pFile, int type)
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
PATO::WriteWells()

Display well information.
******************************************************************************/
void PATO::WriteWells(FILE * pFile, int type)
{
   int i;
   double x, y, q, lift, depth, H, T, B;
   char activeStr[10];
   
   switch(type)
   {
      case(WRITE_SCI) : 
      {
         for(i = 0; i < m_MaxNumWells; i++)
         {
            x = m_pWells[i].pXloc->GetEstVal();
            y = m_pWells[i].pYloc->GetEstVal();
            q = m_pWells[i].pQ->GetEstVal();
            if((q < -m_RateThresh) || (q > m_RateThresh)){strcpy(activeStr, "YES");}
            else{strcpy(activeStr, "NO");}
            if(m_pWells[i].pHead != NULL)
            {
               H = m_pWells[i].pHead->GetCurrentVal();

               if(m_pWells[i].pTopo != NULL){ T = m_pWells[i].pTopo->GetCurrentVal();}
               else{ T = m_pWells[i].Topo;}

               if(m_pWells[i].pBase != NULL){ B = m_pWells[i].pBase->GetCurrentVal();}
               else{ B = m_pWells[i].Base;}

               lift = T - H; 
               depth = T - B;            
               fprintf(pFile, "%-16s  %-6s  %E  %E  %E  ", 
                       m_pWells[i].name, activeStr, x, y, q);
               fprintf(pFile, "%E  %E  %E  %E  %E  ",
                       H, T, B, lift, depth);
            }
            else
            {
               fprintf(pFile, "%-12s  %-6s  %E  %E  %E  ", 
                       m_pWells[i].name, activeStr, x, y, q);
               fprintf(pFile, "n/a         n/a         n/a         n/a         n/a           ");
            }
            //cost information depends on type of cost function
            if(m_ObjType == PATO_OBJ_RATE)
            {
               fprintf(pFile, "n/a         n/a         n/a         %E\n", m_pWells[i].Ctot);
            }
            else if(m_ObjType == PATO_OBJ_OP)
            {
               fprintf(pFile, "n/a         n/a         %E  %E\n", m_pWells[i].Cnrg, m_pWells[i].Ctot);
            }
            else
            {
               fprintf(pFile, "%E  %E  %E  %E\n", m_pWells[i].Cdrill, m_pWells[i].Cpump, 
                       m_pWells[i].Cnrg, m_pWells[i].Ctot);
            }
         }/* end for() */
         break;
      }/* end case(WRITE_SCI) */
      case(WRITE_DEC) :
      {
         for(i = 0; i < m_MaxNumWells; i++)
         {
            x = m_pWells[i].pXloc->GetEstVal();
            y = m_pWells[i].pYloc->GetEstVal();
            q = m_pWells[i].pQ->GetEstVal();
            if((q < -m_RateThresh) || (q > m_RateThresh)){strcpy(activeStr, "YES");}
            else{strcpy(activeStr, "NO");}
            if(m_pWells[i].pHead != NULL)
            {
               H = m_pWells[i].pHead->GetCurrentVal();

               if(m_pWells[i].pTopo != NULL){ T = m_pWells[i].pTopo->GetCurrentVal();}
               else{ T = m_pWells[i].Topo;}

               if(m_pWells[i].pBase != NULL){ B = m_pWells[i].pBase->GetCurrentVal();}
               else{ B = m_pWells[i].Base;}

               lift = T - H; 
               depth = T - B;

               fprintf(pFile, "%-12s  %-6s  %-10.3lf  %-10.3lf  %-10.3lf  ", 
                       m_pWells[i].name, activeStr, x, y, q);
               fprintf(pFile, "%-10.3lf  %-10.3lf  %-10.3lf  %-10.3lf  %-10.3lf  ",
                       H, T, B, lift, depth);
            }
            else
            {
               fprintf(pFile, "%-12s  %-6s  %-10.3lf  %-10.3lf  %-10.3lf  ", 
                       m_pWells[i].name, activeStr, x, y, q);
               fprintf(pFile, "n/a         n/a         n/a         n/a         n/a           ");
            }
            //cost information depends on type of cost function
            if(m_ObjType == PATO_OBJ_RATE)
            {
               fprintf(pFile, "n/a         n/a         n/a         %-10.2lf\n", 
                       m_pWells[i].Ctot);
            }
            else if(m_ObjType == PATO_OBJ_OP)
            {
               fprintf(pFile, "n/a         n/a         %-10.2lf  %-10.2lf\n", 
                       m_pWells[i].Cnrg, m_pWells[i].Ctot);
            }
            else
            {
               fprintf(pFile, "%-10.2lf  %-10.2lf  %-10.2lf  %-10.2lf\n", 
                       m_pWells[i].Cdrill, m_pWells[i].Cpump, m_pWells[i].Cnrg, 
                       m_pWells[i].Ctot);
            }
         }/* end for() */
         break;
      }/* end case(WRITE_DEC) */
      case(WRITE_BNR) :
      {
	      fprintf(pFile, "Name           Active?  X-loc       Y-loc       Rate        ");
         fprintf(pFile, "Head        Surface     Base        Lift        Depth       ");
         fprintf(pFile, "Drill Cost  Pump Cost   Energy      Total Cost\n");
         break;
      }/* end case(WRITE_BNR) */
      default:
      case(WRITE_DBG) :
      {
         for(i = 0; i < m_MaxNumWells; i++)
         {
            q = m_pWells[i].pQ->GetEstVal();
            if((q < -m_RateThresh) || (q > m_RateThresh)){strcpy(activeStr, "YES");}
            else{strcpy(activeStr, "NO");}
            fprintf(pFile, "\n***** Well[%d] Information *****\n", i);
            fprintf(pFile, "Name : %s\n", m_pWells[i].name);
            fprintf(pFile, "---------- x-location ----------\n");
            m_pWells[i].pXloc->Write(pFile, type);
            fprintf(pFile, "\n---------- y-location ----------\n");
            m_pWells[i].pYloc->Write(pFile, type);   
            fprintf(pFile, "\n---------- rate       ----------\n");
            m_pWells[i].pQ->Write(pFile, type);
            fprintf(pFile, "Active? %s\n", activeStr);
            if(m_pWells[i].pHead != NULL)
            {
               fprintf(pFile, "\n---------- head ----------\n");
               m_pWells[i].pHead->Write(pFile, type);
               fprintf(pFile, "\n---- surface elevation ----\n");
               if(m_pWells[i].pTopo != NULL){ m_pWells[i].pTopo->Write(pFile, type);}
               else{ fprintf(pFile, "%lf\n", m_pWells[i].Topo);}
               fprintf(pFile, "\n---- aquifer base ----\n");
               if(m_pWells[i].pBase != NULL){ m_pWells[i].pBase->Write(pFile, type);}
               else{ fprintf(pFile, "%lf\n", m_pWells[i].Base);}
            }
            fprintf(pFile, "Drill Cost  : %lf\n", m_pWells[i].Cdrill); 
            fprintf(pFile, "Pump  Cost  : %lf\n", m_pWells[i].Cpump); 
            fprintf(pFile, "Energy Cost : %lf\n", m_pWells[i].Cnrg); 
            fprintf(pFile, "Total  Cost : %lf\n", m_pWells[i].Ctot); 
         }
         break;
      }/* end case(WRITE_DBG) */
   }/* end switch() */
}/* end WriteWells() */

/******************************************************************************
PATO::LookupPumpCost()

Find the cheapest pump that will operate at the desired rate and lift and 
return its cost.
******************************************************************************/
double PATO::LookupPumpCost(double rate, double lift)
{
   double cur, min, rateLow, rateHi, liftLow, liftHi;
   int i;

   if(m_pTbl == NULL){ return 0.00;}

   //init. min to be max cost
   min = 0.00;
   for(i = 0; i < m_TblSize; i++)
   { 
      cur = m_pTbl[i].cost;
      if(cur > min){ min = cur;}
   }

   //find cheapest match in the table
   for(i = 0; i < m_TblSize; i++)
   {
      cur     = m_pTbl[i].cost;
      rateLow = m_pTbl[i].Qmin;
      rateHi  = m_pTbl[i].Qmax;
      liftLow = m_pTbl[i].Lmin;
      liftHi  = m_pTbl[i].Lmax;

      if((rate >= rateLow) && (rate <= rateHi) && 
         (lift >= liftLow) && (lift <= liftHi) &&
         (cur < min))
      {
         min = cur;
      }
   }
   return min;
}/* end LookupPumpCost() */


/******************************************************************************
PATO::GetConstraintPtr()

Retrieve constraint associated with pName
******************************************************************************/
ConstraintABC * PATO::GetConstraintPtr(IroncladString pName)
{
   ConstraintABC * pCur;

   for(pCur = m_pConstraints; pCur != NULL; pCur = pCur->GetNext());
   {
      if(strcmp(pCur->GetName(), pName) == 0) return pCur;
   }	
   return NULL;
}/* end GetConstraintPtr() */

