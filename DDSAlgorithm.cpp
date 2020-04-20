/******************************************************************************
File     : DDSAlgorithm.cpp
Author   : James Craig and L. Shawn Matott
Copyright: 2009, James Craig and L. Shawn Matott

An implementation of the DDS algorithm.

Version History
09-09-09    lsm   Created
******************************************************************************/
#include <math.h>
#include <string.h>

#include "DDSAlgorithm.h"
#include "Model.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "StatsClass.h"

#include "Utility.h"
#include "WriteUtility.h"
#include "Exception.h"

/******************************************************************************
WarmStart()

Read the best solution from a previous run.
******************************************************************************/
void DDSAlgorithm::WarmStart(void)
{
   int np = m_pModel->GetParamGroupPtr()->GetNumParams();
   double * pbest = new double[np+1];
   int newcount = SimpleWarmStart(np, pbest);
   m_pModel->GetParamGroupPtr()->WriteParams(pbest);
   ((Model *)m_pModel)->SetCounter(newcount);
   delete [] pbest;
}/* end WarmStart() */

/**********************************************************************
CTOR
	Allocates memory for DDS class members
**********************************************************************/
DDSAlgorithm::DDSAlgorithm(ModelABC *pModel)
{
	FILE * inFile;
   char * line;
   char tmp[DEF_STR_SZ];
   char begin_token[DEF_STR_SZ];
   char end_token[DEF_STR_SZ];
	
	RegisterAlgPtr(this);

	m_pModel = pModel;
  m_pStats = NULL;
	
	//init. everything to reasonable defaults
   m_r_val = 0.2;
	m_UserSeed = GetRandomSeed();
	m_MaxIter = 100;
	m_UserSuppliedInit = false;

	//Read Data From Algorithm Input file 
	IroncladString pFileName = GetInFileName();
	inFile = fopen(pFileName, "r");  
  if(inFile == NULL) 
  {
    FileOpenFailure("DDSAlgorithm::CTOR", pFileName);
  }/* end if() */ 

  //accept multiple section headings
  if(CheckToken(inFile, "BeginDDSAlg", pFileName) == true)
  {
     strcpy(begin_token, "BeginDDSAlg");
     strcpy(end_token, "EndDDSAlg");
  }
  else
  {
     rewind(inFile);
     if(CheckToken(inFile, "BeginDDS", pFileName) == true)
     {
        strcpy(begin_token, "BeginDDS");
        strcpy(end_token, "EndDDS");
     }
     else //default
     {
        strcpy(begin_token, "BeginDDS");
        strcpy(end_token, "EndDDS");
     }
   }
   rewind(inFile);

  if(CheckToken(inFile, begin_token, pFileName) == true)
  {
      FindToken(inFile, end_token, pFileName);
      rewind(inFile);

      FindToken(inFile, begin_token, pFileName);
      line = GetNxtDataLine(inFile, pFileName);
      
      while(strstr(line, end_token) == NULL)
      {
         if (strstr(line, "PerturbationValue") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_r_val);
         }
         else if(strstr(line, "MaxIterations") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_MaxIter);
         }
         else if(strstr(line, "UseInitialParamValues") != NULL)
         {
            m_UserSuppliedInit=true;
         }
         else if(strstr(line, "UseRandomParamValues") != NULL)
         {
            m_UserSuppliedInit=false;
         }
         line = GetNxtDataLine(inFile, pFileName);
      }/* end while() */
   } /* end if() */
   else
   {
      if(GetProgramType() != DDSAU_PROGRAM)
      {
         LogError(ERR_FILE_IO, "Using default DDS algorithm setup.");
      }
   }/* end else() */

	 if ((m_r_val<0) || (m_r_val>1)){
     LogError(ERR_FILE_IO,"Bad Perturbation value specified for DDS Algorithm");
		 ExitProgram(1);
	 }
	 if (m_MaxIter<1){
     LogError(ERR_FILE_IO,"Maximum DDS Algorithm iterations must be >0");
		 ExitProgram(1);
	 }
   fclose(inFile);

  IncCtorCount();
}/* end CTOR() */

/**********************************************************************
Destroy

  Free up memory.
**********************************************************************/
void   DDSAlgorithm::Destroy()
{
  delete m_pStats;
  IncDtorCount();
}/* end Destroy() */

/**********************************************************************
ResetUserSeed()

**********************************************************************/
void DDSAlgorithm::ResetUserSeed(int seed)
{
   m_UserSeed = seed; 
   ResetRandomSeed(seed);
}/* end ResetUserSeed() */

/**********************************************************************
Optimize
**********************************************************************/
void   DDSAlgorithm::Optimize() 
{
	//------------------------------------------------------------
	//				VARIABLE DECLARATION
	//------------------------------------------------------------
	double	Ftest,Fbest, * Cbest;
	double  Pn,convergence;
	int     dv,dvn_count;			//index, counter for how many dec vars vary in neighbourhood
	int     iters_remaining,InitFunctEvals;
	int	  i,j, k;							//counter variables
	int     NumParams;
	double *BestParams,*TestParams;
   bool bBanner, bWarmStart;
   double a = 0.00; //the fraction of elapsed budget
   FILE * pPnFile;

	StatusStruct    pStatus;
	ParameterABC   *pParam;
	ParameterGroup *pParamGroup;

	pParamGroup = m_pModel->GetParamGroupPtr(); 
	NumParams = pParamGroup->GetNumParams();
   int nSpecial = pParamGroup->GetNumSpecialParams();
   Cbest =  new double[nSpecial];
	
	NEW_PRINT("DDSMembers", NumParams);
	BestParams		=new double [NumParams];
	TestParams		=new double [NumParams];
	MEM_CHECK(TestParams);

	for (k=0; k<NumParams;k++)
	{
		BestParams[k]=TestParams[k]=pParamGroup->GetParamPtr(k)->GetEstVal();		
	}

   //write setup
   WriteSetup(m_pModel, "Dynamically Dimensioned Search Algorithm (DDS)");   
   //write banner
   WriteBanner(m_pModel, "trial    best fitness   ", " trials remaining");
	pStatus.maxIter =m_MaxIter;

	//------------------------------------------------------------
	//				INITIALIZATION
	//------------------------------------------------------------

   //read in best result from previous run, if desired
   bWarmStart = m_pModel->CheckWarmStart();
   if(bWarmStart == true)
   {
      InitFunctEvals=1;
      WarmStart();
		for(k=0;k<NumParams;k++)
      {
         //Estimated value=warm start solution
			TestParams[k]=pParamGroup->GetParamPtr(k)->GetEstVal(); 
		}
   }
	//user supplied fixed initial solution------------------------
	else if (m_UserSuppliedInit)
	{
		InitFunctEvals=1;
		for(k=0;k<NumParams;k++)
      {
         //Estimated value=initial solution
			TestParams[k]=pParamGroup->GetParamPtr(k)->GetEstVal(); 
		}
		// printf("Evaluating user supplied initial solution....  \n");
	}   
   //standard random initialization------------------------------
	else                 
	{
		//Calculate # function evals to use for random sampling to initialize DDS 
		InitFunctEvals=iMax(5,(int)(0.005*double(m_MaxIter)));
      if(InitFunctEvals >= m_MaxIter)
      {
         InitFunctEvals = 1;
      }
		// printf("Sampling for initial DDS solution....   \n");
	}

	iters_remaining=(m_MaxIter-InitFunctEvals); // reduce number of fevals in DDS loop
	if (iters_remaining<=0) 
   {
		LogError(ERR_FILE_IO,"DDSAlgorithm: # of Initialization samples >= Max # func evaluations");
		ExitProgram(1);
	}

	//Use random sampling (rs) algorithm to initialize DDS ------------------------------------
	//note:  The only way to start runs from same solution for a different
	//       m_MaxIter limit is to read in the initial solution matrix from an input file
   j = 0;
   m_CurIter = 0;
	for(i=1;i<=InitFunctEvals;i++)
	{
		pStatus.curIter = i;      

      if(IsQuit() == true){ break;}

      // sample an initial solution candidate-overwrites initialization of TestParams
		if((m_UserSuppliedInit == false) && (bWarmStart == false))
		{
			for(k=0;k<NumParams;k++)
			{
				pParam  = pParamGroup->GetParamPtr(k); 
				TestParams[k]=UniformRandom()*(pParam->GetUprBnd()-pParam->GetLwrBnd())+pParam->GetLwrBnd();
			}
		} 

		// Evaluate obj function for initial solution guess   	
      MakeParameterCorrections(TestParams, BestParams, NumParams, a);		
		m_pModel->Execute();
      a += 1.00/(double)m_MaxIter;
		Ftest=m_pModel->GetObjFuncVal();
      m_CurIter++;

	  if ((i==1) || (Ftest<=Fbest))  // if initializing or solution is best, update solution
	  {
         if(i != 1)
         {
            WriteInnerEval(++j, 0, '.');
            WriteInnerEval(WRITE_ENDED, 0, '.');
         }
         else
         {
            /* --------------------------------------------
            enable special parameters now that the best 
            parameter set is initialized.
            -------------------------------------------- */
            pParamGroup->EnableSpecialParams();
         }
			Fbest = Ftest;
			for (k=0;k<NumParams;k++){BestParams[k]=TestParams[k];}
			convergence=(double)(m_MaxIter-i);//TMP DEBUG - should be relative convergence
			WriteRecord(m_pModel, i, Fbest, convergence);
         bBanner = true;
         m_pModel->SaveBest(0); //save the input and output files of the best configuration

         //update special parameters (i.e. pre-emption)
         pParamGroup->GetSpecialConstraints(Cbest);
			pParamGroup->ConfigureSpecialParams(Fbest, Cbest);
		}
      else
      {
         if(bBanner == true)
         {
            WriteInnerEval(WRITE_DDS, 0, '.');
            bBanner = false;
            j=0;
         }
         WriteInnerEval(++j, 0, '.');
      }
	}	// end for(i=1;i<=InitFunctEvals;i++) (initialization rs sampling loop)

  // printf("Done DDS initialization. DDS running...");

	//------------------------------------------------------------
	//				MAIN DDS LOOP
	//------------------------------------------------------------

   pPnFile = fopen("OstDDSPn.txt", "w");
   fprintf(pPnFile, "EVAL  Pn\n");
   fclose(pPnFile);

	for(i=1;i<=iters_remaining;i++) 
	{
		pStatus.curIter = i+InitFunctEvals;      

      if(IsQuit() == true){ break;}

		// Determine variable selected as neighbour 
		Pn=1.0-log(double(i))/log(double(iters_remaining)); 
      pPnFile = fopen("OstDDSPn.txt", "a");
      fprintf(pPnFile, "%04d  %E\n", i, Pn);
      fclose(pPnFile);

		dvn_count=0; 

		// define TestParams initially as best current solution
		for (k=0;k<NumParams;k++)
      {
			TestParams[k]=BestParams[k];
		}

		for(k=0;k<NumParams;k++)// randomly select DVs in neighbourhood to perturb 
		{	 
			if (UniformRandom()<Pn)
         { 
            dvn_count=dvn_count+1;
				TestParams[k]=PerturbParam(BestParams[k],pParamGroup->GetParamPtr(k));
			}
		}
		if (dvn_count==0)			  // if no DVs selected at random, select ONE
		{ 
				dv=(int)(ceil((double)(NumParams)*UniformRandom()))-1; // index for one DV 
				TestParams[dv]=PerturbParam(BestParams[dv],pParamGroup->GetParamPtr(dv));
		}
		
		// Evaluate obj function for current set of DVs
      MakeParameterCorrections(TestParams, BestParams, NumParams, a);
		m_pModel->Execute();
      a += 1.00/(double)m_MaxIter;
		Ftest=m_pModel->GetObjFuncVal();
      m_CurIter++;

		if (Ftest<=Fbest) // update current (best) solution
		{							 
			Fbest = Ftest;
			for (k=0;k<NumParams;k++){BestParams[k]=TestParams[k];}
			//	cout<< f_count[ind1]<<" "<< to_max*Fbest<<endl;
   
		   //write results
         WriteInnerEval(++j, 0, '.');
         WriteInnerEval(WRITE_ENDED, 0, '.');
		   convergence=(double)(m_MaxIter-(i+InitFunctEvals));//TMP DEBUG - should be relative convergence
         pParamGroup->WriteParams(BestParams);
		   WriteRecord(m_pModel, i+InitFunctEvals, Fbest, convergence);
		   pStatus.pct  = ((float)(100)*(float)(i+InitFunctEvals))/(float)(m_MaxIter);
		   pStatus.numRuns = m_pModel->GetCounter();
		   WriteStatus(&pStatus);
         bBanner = true;
         m_pModel->SaveBest(0); //save the input and output files of the best configuration

         //update special parameters (i.e. pre-emption)
         pParamGroup->GetSpecialConstraints(Cbest);
			pParamGroup->ConfigureSpecialParams(Fbest, Cbest);
      }
      else
      {
         if(bBanner == true)
         {
            WriteInnerEval(WRITE_DDS, 0, '.');
            bBanner = false;
            j=0;
         }
         WriteInnerEval(++j, 0, '.');
      }

      if(i == iters_remaining)
      {
         convergence=(double)(m_MaxIter-(i+InitFunctEvals));//TMP DEBUG - should be relative convergence
         WriteInnerEval(WRITE_ENDED, 0, '.');
         pParamGroup->WriteParams(BestParams);
		   WriteRecord(m_pModel, i+InitFunctEvals, Fbest, convergence);
      }
	}/* end for(i=1;i<=iters_remaining;i++)  (main DDS loop)*/

	pParamGroup->WriteParams(BestParams); 
   m_pModel->Execute();
   WriteOptimal(m_pModel, Fbest);
   m_pModel->SaveBest(0); //save the input and output files of the best configuration

	pStatus.pct = 100.0;
   pStatus.numRuns = m_pModel->GetCounter();
   WriteStatus(&pStatus);
   WriteAlgMetrics(this);

	delete [] BestParams;
	delete [] TestParams;
   delete [] Cbest;
}

/*************************************************************************************
MakeParameerCorrections()
*************************************************************************************/
void DDSAlgorithm::MakeParameterCorrections(double * x, double * xb, int n, double a)
{
   double lwr, upr;
   ParameterGroup * pParamGroup;
	pParamGroup = m_pModel->GetParamGroupPtr(); 

   for(int k = 0; k < n; k++)
   {
      lwr=m_pModel->GetParamGroupPtr()->GetParamPtr(k)->GetLwrBnd();
      upr=m_pModel->GetParamGroupPtr()->GetParamPtr(k)->GetUprBnd();
      x[k]=TelescopicCorrection(lwr, upr, xb[k], a, x[k]);
   }
   pParamGroup->WriteParams(x); 		

   //inerface with expert judgement module
   m_pModel->PerformParameterCorrections();
   for(int i = 0; i < n; i++)
   {
      x[i] = pParamGroup->GetParamPtr(i)->GetEstVal();
   }/* end for() */

}/* enad MakeParameterCorrections() */

/**********************************************************************
		PerturbParam
-----------------------------------------------------------------------
	Generates a neighboring decision variable value for a single
	decision variable value being perturbed by the DDS optimization algorithm.
	New DV value respects the upper and lower DV bounds.
	
	Translated by James Craig July 2006 from Bryan Tolson's Fortran DDS code
 
  returns  new decision variable value (within specified min and max)
**********************************************************************/
double DDSAlgorithm::PerturbParam(const double &x_best, //current best decision variable (DV) value
															    ParameterABC * pParam)     
{

	double x_new, r;

	double x_max = pParam->GetUprBnd();
	double x_min = pParam->GetLwrBnd();

   if(strcmp(pParam->GetType(), "real") == 0)
   {
	  x_new=x_best+GaussRandom()*m_r_val*(x_max-x_min);
   }
   else if(strcmp(pParam->GetType(), "integer") == 0)
   {
     r = 2.00*UniformRandom()-1.00; //-1 < r < +1
     x_new=x_best+r*m_r_val*(x_max-x_min);
     //check that a discrete change has been made.
     if((int)x_new == (int)x_best)
     {
        if(x_new<x_best) x_new = x_new-1.00;  //decrement
        else x_new = x_new+1.00; //increment
     }/* end if() */
     //bounds check
     if((int)x_new > (int)x_max) x_new = x_best - 1.00;
     if((int)x_new < (int)x_min) x_new = x_best + 1.00;
   }
   else
   {
     x_new=x_best+GaussRandom()*m_r_val*(x_max-x_min);
   }

	// need if statements to check within DV bounds.  If not, bounds are reflecting.
	if (x_new<x_min) 
	{
		x_new=x_min+(x_min-x_new); //reflect
		// if reflection goes past x_max then value should be x_min since 
		// without reflection the approach goes way past lower bound.  
		// This keeps x close to lower bound when x_best is close to lower bound
		if (x_new>x_max) {x_new=x_min;}
	}
	else if (x_new>x_max) 
	{
		x_new=x_max-(x_new-x_max);//reflect

		if (x_new<x_min){x_new=x_max;}
	}

	return x_new;
}
/******************************************************************************
Calibrate()

Solve the Least-Squares minimization problem using the DDS
******************************************************************************/
void DDSAlgorithm::Calibrate(void)
{ 
	char fileName[DEF_STR_SZ];
	FILE * pFile;
   
   NEW_PRINT("StatsClass", 1);
   m_pStats = new StatsClass(m_pModel);
   MEM_CHECK(m_pStats);
   RegisterStatsPtr(m_pStats);

	Optimize();
	
  sprintf(fileName, "OstOutput%d.txt", 0);

  //compute statistics (variance and covariance)
  m_pStats->CalcStats();

  //write statistics of best parameter set to output file
  pFile = fopen(fileName, "a");   
  m_pStats->WriteStats(pFile);
  fclose(pFile);

  //write statistics of best parameter set to output file
  m_pStats->WriteStats(stdout);

} /* end Calibrate() */
/**********************************************************************
		WriteMetrics
-----------------------------------------------------------------------
**********************************************************************/
void DDSAlgorithm::WriteMetrics(FILE * pFile)
{
   fprintf(pFile, "\nAlgorithm Metrics\n");
   fprintf(pFile, "Algorithm               : Dynamically-Dimensioned Search Algorithm (DDS)\n");
   fprintf(pFile, "Desired Convergence Val : N/A\n");
   fprintf(pFile, "Actual Convergence Val  : N/A\n");
   fprintf(pFile, "Max Generations         : %d\n", m_MaxIter);
   fprintf(pFile, "Actual Generations      : %d\n", m_MaxIter);
   fprintf(pFile, "Peterbation Value       : %lf\n", m_r_val);
   m_pModel->WriteMetrics(pFile);
   fprintf(pFile, "Algorithm successfully converged on a solution, however more runs may be needed\n");
}/* end WriteMetrics() */

/******************************************************************************
DDS_Program()
Calibrate the model using DDS.
******************************************************************************/
void DDS_Program(int argC, StringType argV[])
{
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model();

   NEW_PRINT("DDS", 1);
   DDSAlgorithm * DDS = new DDSAlgorithm(model);

   MEM_CHECK(DDS);
   
   if(model->GetObjFuncId() == OBJ_FUNC_WSSE) { DDS->Calibrate(); }
   else { DDS->Optimize(); }

   delete DDS;
   delete model;
} /* end DDS_Program() */

