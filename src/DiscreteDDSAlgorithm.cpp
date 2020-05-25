/******************************************************************************
File     : DiscreteDDSAlgorithm.cpp
Author   : L. Shawn Matott
Copyright: 2011, L. Shawn Matott

An implementation of the discrete DDS algorithm.

Version History
11-11-11    lsm   Created
******************************************************************************/
#include <math.h>
#include <string.h>

#include "DiscreteDDSAlgorithm.h"
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
void DiscreteDDSAlgorithm::WarmStart(void)
{
   int np = m_pModel->GetParamGroupPtr()->GetNumParams();
   double * pbest = new double[np+1];
   int newcount = SimpleWarmStart(np, pbest);
   m_pModel->GetParamGroupPtr()->WriteParams(pbest);
   ((Model *)m_pModel)->SetCounter(newcount);
   delete [] pbest;
}/* end WarmStart() */

/**********************************************************************
		DDS_Initialize
-----------------------------------------------------------------------
	Allocates memory for DDS class members
**********************************************************************/
DiscreteDDSAlgorithm::DiscreteDDSAlgorithm(ModelABC *pModel)
{
	FILE * inFile;
  char * line;
  char tmp[DEF_STR_SZ];
	
	RegisterAlgPtr(this);

	m_pModel = pModel;
  m_pStats = NULL;
	
	//init. everything to reasonable defaults
   m_r_val						= 0.2;
	m_UserSeed				   = GetRandomSeed();
	m_MaxIter					= 100;
	m_UserSuppliedInit = false;
   m_nCorr = 0;

	//Read Data From Algorithm Input file 
	IroncladString pFileName = GetInFileName();
	inFile = fopen(pFileName, "r");  
  if(inFile == NULL) 
  {
    FileOpenFailure("DDSAlgorithm::CTOR", pFileName);
  }/* end if() */ 

  if(CheckToken(inFile, "BeginDiscreteDDSAlg", pFileName) == true)
  {
      FindToken(inFile, "EndDiscreteDDSAlg", pFileName);
      rewind(inFile);

      FindToken(inFile, "BeginDiscreteDDSAlg", pFileName);
      line = GetNxtDataLine(inFile, pFileName);
      
      while(strstr(line, "EndDiscreteDDSAlg") == NULL)
      {
         if     (strstr(line, "PerturbationValue") != NULL)
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
      LogError(ERR_FILE_IO, "Using default DDS algorithm setup.");
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
}

/**********************************************************************
		Destroy
**********************************************************************/
void   DiscreteDDSAlgorithm::Destroy()
{
  delete m_pStats;
  IncDtorCount();
}
/**********************************************************************
		Optimize
-----------------------------------------------------------------------
**********************************************************************/
void   DiscreteDDSAlgorithm::Optimize() 
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

	StatusStruct    pStatus;
	ParameterABC   *pParam;
	ParameterGroup *pParamGroup;

	pParamGroup=m_pModel->GetParamGroupPtr(); 
	NumParams  =pParamGroup->GetNumParams();
   int nSpecial = pParamGroup->GetNumSpecialParams();
   Cbest = new double[nSpecial];
	
	NEW_PRINT("DDSMembers", NumParams);
	BestParams		=new double [NumParams];
	TestParams		=new double [NumParams];
	MEM_CHECK(TestParams);

	for (k=0; k<NumParams;k++)
	{
		BestParams[k]=TestParams[k]=pParamGroup->GetParamPtr(k)->GetEstVal();		
	}

  //write setup
  WriteSetup(m_pModel, "Discrete-valued Dynamically Dimensioned Search Algorithm (DDDS)");   
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
		// printf("Sampling for initial DDS solution....   \n");
	}

	iters_remaining=(m_MaxIter-InitFunctEvals); // reduce number of fevals in DDS loop
	if (iters_remaining<=0) {
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
		pParamGroup->WriteParams(TestParams); 
		m_pModel->Execute();
		Ftest=m_pModel->GetObjFuncVal();
      m_CurIter++;

	  if ((i==1) || (Ftest<=Fbest))  // if initializing or solution is best, update solution
	  {
         if(i != 1)
         {
            WriteInnerEval(++j, 0, '.');
            WriteInnerEval(WRITE_ENDED, 0, '.');
         }
			Fbest = Ftest;
			for (k=0;k<NumParams;k++){BestParams[k]=TestParams[k];}
			convergence=(double)(m_MaxIter-i);//TMP DEBUG - should be relative convergence
			WriteRecord(m_pModel, i, Fbest, convergence);
         bBanner = true;
         pParamGroup->EnableSpecialParams();
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

	for(i=1;i<=iters_remaining;i++) 
	{
		pStatus.curIter = i+InitFunctEvals;      
    if(IsQuit() == true){ break;}

		// Determine variable selected as neighbour 
		Pn=1.0-log(double(i))/log(double(iters_remaining)); 
		dvn_count=0; 

		// define TestParams initially as best current solution
		for (k=0;k<NumParams;k++){
			TestParams[k]=BestParams[k];
		}

		for(k=0;k<NumParams;k++)// randomly select DVs in neighbourhood to perturb 
		{	 
			if (UniformRandom()<Pn) { 
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
		pParamGroup->WriteParams(TestParams); 		
		m_pModel->Execute();
		Ftest=m_pModel->GetObjFuncVal();
      m_CurIter++;

		if (Ftest<=Fbest) // update current (best) solution
		{							 
			Fbest = Ftest;
         pParamGroup->GetSpecialConstraints(Cbest);
         pParamGroup->ConfigureSpecialParams(Fbest, Cbest);

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

	pStatus.pct			= 100.0;
   pStatus.numRuns = m_pModel->GetCounter();
   WriteStatus(&pStatus);
   WriteAlgMetrics(this);

	delete [] BestParams;
	delete [] TestParams;
   delete [] Cbest;
}
/**********************************************************************
		PerturbParam
-----------------------------------------------------------------------
	Generates a neighboring decision variable value for a single
	decision variable value being perturbed by the DDS optimization algorithm.
	New DV value respects the upper and lower DV bounds and must be at least
   +/-1 (since the parameters are assumed to be discrete-valued).
	 
  Returns  new decision variable value (within specified min and max)
**********************************************************************/
double DiscreteDDSAlgorithm::PerturbParam(const double &x_best, //current best decision variable (DV) value
															    ParameterABC * pParam)     
{

	double x_new;

	double x_max = pParam->GetUprBnd();
	double x_min = pParam->GetLwrBnd();

	x_new=x_best+GaussRandom()*m_r_val*(x_max-x_min);

   //check that a discrete change has been made.
   if((int)x_new == (int)x_best)
   {
      if(x_new<x_best) x_new = x_new-1.00;  //decrement
      else x_new = x_new+1.00; //increment
      m_nCorr++;
   }/* end if() */

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
void DiscreteDDSAlgorithm::Calibrate(void)
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
void DiscreteDDSAlgorithm::WriteMetrics(FILE * pFile)
{
   fprintf(pFile, "\nAlgorithm Metrics\n");
   fprintf(pFile, "Algorithm                 : Discretely Dynamically-Dimensioned Search Algorithm (DDDS)\n");
   fprintf(pFile, "Desired Convergence Val   : N/A\n");
   fprintf(pFile, "Actual Convergence Val    : N/A\n");
   fprintf(pFile, "Max Generations           : %d\n", m_MaxIter);
   fprintf(pFile, "Actual Generations        : %d\n", m_MaxIter);
   fprintf(pFile, "Peterbation Value         : %lf\n", m_r_val);
   fprintf(pFile, "Num. Discrete Corrections : %d\n", m_nCorr);
   m_pModel->WriteMetrics(pFile);
   fprintf(pFile, "Algorithm successfully converged on a solution, however more runs may be needed\n");
}/* end WriteMetrics() */

/******************************************************************************
DDS_Program()
Calibrate the model using DDS.
******************************************************************************/
void DiscreteDDS_Program(int argC, StringType argV[])
{
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model();

   NEW_PRINT("DiscreteDDS", 1);
   DiscreteDDSAlgorithm * DDDS = new DiscreteDDSAlgorithm(model);

   MEM_CHECK(DDDS);
   
   if(model->GetObjFuncId() == OBJ_FUNC_WSSE) { DDDS->Calibrate(); }
   else { DDDS->Optimize(); }

   delete DDDS;
   delete model;
} /* end DiscreteDDS_Program() */

