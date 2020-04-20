/******************************************************************************
File     : PADDS.cpp
Author   : L. Shawn Matott
Copyright: 2015, L. Shawn Matott

The PADDS (Pareto Archived Dynamically Dimensioned Search) algorithm is a 
multi-objective version of the DDS algorithm. It has been ported from a C++ 
implementation of Mohammadamin Jahanpour (mjahanpo@uwaterloo.ca).

Version History
05-07-15    lsm   created file
******************************************************************************/

/*------------------------------------------------------------------- 
The copyright, disclaimer, and citation given below affects the 
following list of subroutines:

   bool covers(double* cub, double * regLow);
   bool partCovers(double* cub, double * regUp);
   int containsBoundary(double* cub, double * regLow, int split);
   double getMeasure(double * regLow, double * regUp);
   int isPile(double* cub, double * regLow, double * regUp);
   double getMedian(double * bounds, int size);
   double computeTrellis(double * regLow, double * regUp, double * trellis);
   void stream(double * regionLow, double * regionUp, double ** points, int npoints, int split, double cover);

Copyright (c) 2006 Nicola Beume  This program is free software: you 
can redistribute it and/or modify it under the terms of the GNU 
General Public License as published by the Free Software Foundation, 
either version 2 of the License, or (at your option) any later version. 
This program is distributed in the hope that it will be useful, but 
WITHOUT ANY WARRANTY; without even the implied warranty of 
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
Public License for more details. You should have received a copy of the 
GNU General Public License along with this program.

--------------------------------------------------------------------- 

This program calculates the dominated hypervolume or S-metric of a set 
of d-dimensional points (d>=3). Please refer to the following publication 
for a description of the algorithm: 
   Nicola Beume and Guenter Rudolph. Faster S-Metric Calculation by 
   Considering Dominated Hypervolume as Klee's Measure Problem. 
   In: B. Kovalerchuk (ed.): Proceedings of the Second IASTED Conference 
   on Computational Intelligence (CI 2006), pp. 231-236. 
   ACTA Press: Anaheim, 2006. 

   Extended version published as: Technical Report of the Collaborative 
   Research Centre 531 'Computational Intelligence', CI-216/06, 
   ISSN 1433-3325. University of Dortmund, July 2006. 
-------------------------------------------------------------------*/
#include <string.h>

#include "PADDS.h"
#include "Model.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "ObjectiveFunction.h"

#include "Utility.h"
#include "WriteUtility.h"
#include "Exception.h"

/******************************************************************************
CTOR

Registers the algorithm pointer and creates instances of member variables.
******************************************************************************/
PADDS::PADDS(ModelABC * pModel)
{
   RegisterAlgPtr(this);
   m_pModel = pModel;
   m_pNonDom = NULL;
   m_pDom = NULL;
   m_NumNonDom = 0;
   m_NumDom = 0;
   m_maxiter = 0;
   m_CurIter = 0;
   m_num_dec = 0;
   m_num_objs = 0;
   m_Select_metric = 3; //make "exact" the default
   m_fraction1 = 0.0;
   m_dominance_flag = 0;
   m_seed = 0;
   m_dim = 0;
   m_dimension = 0;
   m_dSqrtDataNumber = 0;
   m_volume = 0;
   m_pInit = NULL;

   IncCtorCount();
}/* end CTOR() */

/******************************************************************************
Destroy()

Free up memory used by PADDS and it's member variables.
******************************************************************************/
void PADDS::Destroy(void)
{
   ArchiveStruct * pCur, * pDel;

   for(pCur = m_pNonDom; pCur != NULL;)
   {
     pDel = pCur;
     pCur = pCur->pNext;
     delete [] pDel->F;
     delete [] pDel->X;
     delete pDel; 
   }

   for(pCur = m_pDom; pCur != NULL;)
   {
     pDel = pCur;
     pCur = pCur->pNext;
     delete [] pDel->F;
     delete [] pDel->X;
     delete pDel; 
   }

   for(int i = 0; i < m_NumInit; i++)
   {
      delete [] m_pInit[i];
   }
   delete [] m_pInit;

   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
WriteMetrics()

Write out algorithm metrics and setup.
******************************************************************************/
void PADDS::WriteMetrics(FILE * pFile) 
{
   char select_str[DEF_STR_SZ];
   //0: Random
   //1: Crowding distance
   //2: Hypervolume Contribution (ESTIMATE)
   //3: Hypervolume Contribution (EXACT)
   if(m_Select_metric == 0) strcpy(select_str, "random");
   else if(m_Select_metric == 1) strcpy(select_str, "crowding distance");
   else if(m_Select_metric == 2) strcpy(select_str, "estimated hypervolume contribution");
   else if(m_Select_metric == 3) strcpy(select_str, "exact hypervolume contribution");
   else strcpy(select_str, "unknown");

   fprintf(pFile, "\nAlgorithm Metrics\n");
   fprintf(pFile, "Algorithm               : PADDS - Pareto Archived Dynamically Dimensioned Search\n");
   fprintf(pFile, "Max Iterations          : %d\n", m_maxiter);
   fprintf(pFile, "Actual Iterations       : %d\n", m_CurIter);
   fprintf(pFile, "Num Decision Variables  : %d\n", m_num_dec);
   fprintf(pFile, "Num Objectives          : %d\n", m_num_objs);
   fprintf(pFile, "Random Seed             : %d\n", m_seed);
   fprintf(pFile, "Perturbation Value      : %lf\n", m_fraction1);
   fprintf(pFile, "Non-Dominated Solutions : %d\n", m_NumNonDom);  
   fprintf(pFile, "Dominated Solutions     : %d\n", m_NumDom);     
   fprintf(pFile, "Selection Metric        : %s\n", select_str);

   m_pModel->WriteMetrics(pFile);
}/* end WriteMetrics() */

/******************************************************************************
InitFromFile()

Read configuration information from the given filename.
******************************************************************************/
void PADDS::InitFromFile(IroncladString pFileName)
{
   FILE * pFile;
   char * line;
   char tmp[DEF_STR_SZ];
   char begin_token[DEF_STR_SZ];
   char end_token[DEF_STR_SZ];
   
   m_maxiter = 50;
   m_fraction1 = 0.2;
   m_Select_metric = 3; //make "exact" the default

   //read in PADDS configuration
   pFile = fopen(pFileName, "r");
   if(pFile == NULL) 
   {
      //couldn't open file, use defaults and log the error.
      LogError(ERR_FILE_IO, "Couldn't open PADDS config. file. Using Defaults");      
      return;
   }/* end if() */   

  //accept multiple section headings
  if(CheckToken(pFile, "BeginPADDSAlg", pFileName) == true)
  {
     strcpy(begin_token, "BeginPADDSAlg");
     strcpy(end_token, "EndPADDSAlg");
  }
  else
  {
     rewind(pFile);
     if(CheckToken(pFile, "BeginPADDS", pFileName) == true)
     {
        strcpy(begin_token, "BeginPADDS");
        strcpy(end_token, "EndPADDS");
     }
     else //default
     {
        strcpy(begin_token, "BeginPADDS");
        strcpy(end_token, "EndPADDS");
     }
   }
   rewind(pFile);

   //make sure correct tokens are present
   if(CheckToken(pFile, begin_token, pFileName) == true)
   {
      FindToken(pFile, end_token, pFileName);
      rewind(pFile);

      FindToken(pFile, begin_token, pFileName);
      line = GetNxtDataLine(pFile, pFileName);
      while(strstr(line, end_token) == NULL)
      {         
         if (strstr(line, "PerturbationValue") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_fraction1);
         }
         else if(strstr(line, "MaxIterations") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_maxiter);
         }
         else if(strstr(line, "SelectionMetric") != NULL)
         {
            sscanf(line, "SelectionMetric %s", tmp);

            //0: Random
            //1: Crowding distance
            //2: Hypervolume Contribution (ESTIMATE)
            //3: Hypervolume Contribution (EXACT)
            MyStrLwr(tmp);
            if(strcmp(tmp, "random") == 0) m_Select_metric = 0;
            else if(strcmp(tmp, "crowdingdistance") == 0) m_Select_metric = 1; 
            else if(strcmp(tmp, "estimatedhypervolumecontribution") == 0) m_Select_metric = 2; 
            else if(strcmp(tmp, "exacthypervolumecontribution") == 0) m_Select_metric = 3; 
            else m_Select_metric = 3;
         }
         line = GetNxtDataLine(pFile, pFileName);
      } /* end while() */
   }/* end if() */   

   /* read in a list of initial guesses */
   m_NumInit = 0;
   int num = 0;
   rewind(pFile);
   if(CheckToken(pFile, "BeginInitParams", pFileName) == true)
   {
      FindToken(pFile, "EndInitParams", pFileName);
      rewind(pFile);

      //allocate space for the parameter list
      num = m_pModel->GetParamGroupPtr()->GetNumParams();

      //count the number of entries
      FindToken(pFile, "BeginInitParams", pFileName);
      line = GetNxtDataLine(pFile, pFileName);
      while(strstr(line, "EndInitParams") == NULL)
      {
         m_NumInit++;
         line = GetNxtDataLine(pFile, pFileName);
      }/* end while() */

      //allocate space for entries
      if(m_NumInit > 0)
      {
         NEW_PRINT("double *", m_NumInit);
         m_pInit = new double * [m_NumInit];
         MEM_CHECK(m_pInit);
         for(int i = 0; i < m_NumInit; i++)
         { 
            NEW_PRINT("double", num);
            m_pInit[i] = new double[num];
            MEM_CHECK(m_pInit[i]);
         }
      }/* end if() */

      //read in entries
      rewind(pFile);
      FindToken(pFile, "BeginInitParams", pFileName);
      line = GetNxtDataLine(pFile, pFileName);
      int i = 0;
      int j = 0;
      char * pTok = NULL;
      while(strstr(line, "EndInitParams") == NULL)
      {
         pTok = line;
         //extract values, one-by-one, making any necessary conversions
         for(int k = 0; k < num; k++)
         {
            j = ExtractString(pTok, tmp);
            j = ValidateExtraction(j, k, num, "PADDS::InitFromFile()");
            pTok += j;            
            m_pInit[i][k] = m_pModel->GetParamGroupPtr()->GetParamPtr(k)->ConvertInVal(atof(tmp));
         }/* end for() */                  
         i++;
         line = GetNxtDataLine(pFile, pFileName);
      }/* end while() */
   }/* end if() */

   fclose(pFile);
} /* end InitFromFile() */

/******************************************************************************
Calibrate()

Solve the Least-Squares minimization problem using PADDS.
******************************************************************************/
void PADDS::Calibrate(void)
{    
   Optimize();
} /* end Calibrate() */

/******************************************************************************
Optimize()

Search for the pareto front using PADDS.
******************************************************************************/
void PADDS::Optimize(void)
{
   bool bBanner = false;
   StatusStruct pStatus;
   ParameterGroup * pGroup = m_pModel->GetParamGroupPtr();
   FILE * pPnFile;

   InitFromFile(GetInFileName());

   m_num_dec = pGroup->GetNumParams();
   m_num_objs = m_pModel->GetObjFuncPtr()->CalcMultiObjFunc(NULL, -1);
   m_seed = GetRandomSeed();

   WriteSetup(m_pModel, "PADDS - Pareto Archived Dynamically Dimensioned Search");
   //write banner
   WriteBanner(m_pModel, "gen   ", "trials remaining");
	
   //number of initial solutions
	int its = (int)(0.005 * m_maxiter); 
   if(its < 5) its =  5;

   ArchiveStruct * stest, * sbest;

	double * S_min, * S_max; //low and high bounds
   S_min = new double[m_num_dec];
   S_max = new double[m_num_dec];

	for (int i = 0; i < m_num_dec; i++)
   { 
		S_min[i] = pGroup->GetParamPtr(i)->GetLwrBnd();
		S_max[i] = pGroup->GetParamPtr(i)->GetUprBnd();
	}

   //per B. Tolson, user-sepecified initial guesses count towards the budget
   WriteInnerEval(WRITE_USR, m_NumInit, '.');
	for (int i = 0; i < m_NumInit; i++)
   {
      if(IsQuit() == true){ break;}

		stest = new ArchiveStruct;
      stest->nF = m_num_objs;
      stest->nX = m_num_dec;
      stest->F = new double[m_num_objs];
      stest->X = new double[m_num_dec];
      stest->Z = -999.999;
      stest->pNext = NULL;

		for (int j = 0; j < m_num_dec; j++)
      {
			stest->X[j] = m_pInit[i][j];
		}		

      WriteInnerEval(i+1, m_maxiter, '.');
      bBanner = false;
		F(stest);

      int result = UpdateArchive(stest->X, stest->nX, stest->F, stest->nF);

      if(result == ARCHIVE_NON_DOM)
      {
         WriteInnerEval(WRITE_ENDED, 0, '.');
         WriteMultiObjRecord(m_pModel, i+1, m_pNonDom, (double)(m_maxiter-i-1));
         if(i < (m_NumInit - 1))
            WriteInnerEval(WRITE_USR, m_NumInit, '.');
         bBanner = true;
      }
	}/* end user guesses */
   if(bBanner == false)
      WriteInnerEval(WRITE_ENDED, 0, '.');

	//creating initials
   WriteInnerEval(WRITE_SMP, its, '.');
	for (int i = 1; i < its + 1; i++)
   {
      if(IsQuit() == true){ break;}

		stest = new ArchiveStruct;
      stest->nF = m_num_objs;
      stest->nX = m_num_dec;
      stest->F = new double[m_num_objs];
      stest->X = new double[m_num_dec];
      stest->Z = -999.999;
      stest->pNext = NULL;

		for (int j = 0; j < m_num_dec; j++)
      {
			stest->X[j] = S_min[j] + (S_max[j] - S_min[j])*UniformRandom();
		}		

      WriteInnerEval(i+m_NumInit, m_maxiter, '.');
      bBanner = false;
		F(stest);

      int result = UpdateArchive(stest->X, stest->nX, stest->F, stest->nF);

      if(result == ARCHIVE_NON_DOM)
      {
         WriteInnerEval(WRITE_ENDED, 0, '.');
         WriteMultiObjRecord(m_pModel, i+m_NumInit, m_pNonDom, (double)(m_maxiter-i-m_NumInit));
         if(i < its)
            WriteInnerEval(WRITE_SMP, its, '.');
         bBanner = true;
      }
	}/* end initials */
   if(bBanner ==  false)
      WriteInnerEval(WRITE_ENDED, 0, '.');

	//finished creating initials
	int iLeft = m_maxiter - its - m_NumInit;

	//Calculating Selection Metric Z:
	Calc_Z(m_pNonDom);

	//MAIN LOOP

   pPnFile = fopen("OstPADDSPn.txt", "w");
   fprintf(pPnFile, "EVAL  Pn\n");
   fclose(pPnFile);

   WriteInnerEval(WRITE_DDS, iLeft, '.');
 	for (int i = 1; i < iLeft + 1; i++)
   {
		if (m_dominance_flag == -1)
      {
			sbest = SelectFrom(m_pNonDom);
		}
		else
      {
			sbest = m_pNonDom; //Archive[Archive.size() - 1]; 
		}

		// %% DDS
		double Pn = 1 - log10((double)i) / log10((double)iLeft);

      pPnFile = fopen("OstPADDSPn.txt", "a");
      fprintf(pPnFile, "%04d  %E\n", i, Pn);
      fclose(pPnFile);

		int dvn_count = 0;

		stest = new ArchiveStruct;
      stest->nF = m_num_objs;
      stest->nX = m_num_dec;
      stest->F = new double[m_num_objs];
      stest->X = new double[m_num_dec];
      for(int j = 0; j < m_num_dec; j++)
      {
         stest->X[j] = sbest->X[j];
      }
      for(int j = 0; j < m_num_objs; j++)
      {
         stest->F[j] = sbest->F[j];
      }
      stest->Z = sbest->Z;
      stest->pNext = NULL;

		for (int j = 0; j < m_num_dec; j++)
      {
			if (UniformRandom() < Pn)
         {
				dvn_count += 1;
				double new_value = neigh_value_continuous(sbest->X[j], S_min[j], S_max[j], m_fraction1);
				stest->X[j] = new_value;// change relevant dec var value in stest
			}
		}
		if (dvn_count == 0)
      {			
			int dec_var = (int)(ceil(m_num_dec * UniformRandom() - 1.00));
			double new_value = neigh_value_continuous(sbest->X[dec_var], S_min[dec_var], S_max[dec_var], m_fraction1);
			stest->X[dec_var] = new_value;
		}

      WriteInnerEval(i+its+m_NumInit, m_maxiter, '.');
      bBanner = false;
      F(stest);

		//checking to see if x_curr dominates x_new
		bool sbest_dominates_stest = false;

		if(dominion_status(stest, sbest) == 2)
      {
			sbest_dominates_stest = true;
			m_dominance_flag = -1;
         UpdateArchive(stest->X, stest->nX, stest->F, stest->nF);
		}
		else
      {
			//check to see if it is a duplicate, if yes, flag = 0, do not Update_Archive, and do not add
			bool stest_is_duplicate = true;
			for (ArchiveStruct * pCur = m_pNonDom; pCur != NULL; pCur = pCur->pNext)
         {
            for(int j = 0; j < pCur->nF; j++)
            {
				   if (stest->F[j] != pCur->F[j])
               {
					   stest_is_duplicate = false;
					   break;
				   }
            }
            if(stest_is_duplicate == false)
            {
               break;
            }
			}/* end for() */

         if(stest_is_duplicate == true)
         {
            m_dominance_flag = 0;
            DestroyArchive(stest);
         }

			if (stest_is_duplicate == false)
         {
            int result = UpdateArchive(stest->X, stest->nX, stest->F, stest->nF);
            if(result == ARCHIVE_NON_DOM)
            {
               WriteInnerEval(WRITE_ENDED, 0, '.');
               WriteMultiObjRecord(m_pModel, (i+its+m_NumInit), m_pNonDom, (double)(m_maxiter-its-i-m_NumInit));
               if((m_maxiter-its-i-m_NumInit) > 0)
               {
                  WriteInnerEval(WRITE_DDS, 0, '.');
               }
               bBanner = true;
            }
         }

			stest_is_duplicate = false;

			if(m_dominance_flag != -1)
         {
				Calc_Z(m_pNonDom);
			}
		}/* end else(dominion status) */

      pStatus.curIter = (i+its+m_NumInit);
      pStatus.maxIter = m_maxiter;
      pStatus.pct = ((float)100.00*(float)(i+its+m_NumInit))/(float)m_maxiter;
      pStatus.numRuns = (i+its+m_NumInit);
      WriteStatus(&pStatus);
	}/* end for loop */
	
   if(bBanner == false)
   {
      WriteInnerEval(WRITE_ENDED, 0, '.');
      WriteMultiObjRecord(m_pModel, m_maxiter, m_pNonDom, 0.00);
   }

   WriteMultiObjOptimal(m_pModel, m_pNonDom, m_pDom);

   pStatus.pct = 100.00;
   pStatus.numRuns = m_pModel->GetCounter();
   WriteStatus(&pStatus);

   //write algorithm metrics
   WriteAlgMetrics(this);
}/* end main() */

/******************************************************************************
Calc_Z()
******************************************************************************/
void PADDS::Calc_Z(ArchiveStruct * archive)
{
   ArchiveStruct * pA;
   ArchiveStruct ** pSorted;
   int i, archive_size;

   //compute size of archive
   archive_size = 0;
   for(pA = archive; pA != NULL; pA = pA->pNext)
   {   
      archive_size++;
   }
   //make space for sorted archive
   pSorted = new ArchiveStruct *[archive_size];
   //initialize array
   i = 0;
   for(pA = archive; pA != NULL; pA = pA->pNext)
   {
      pSorted[i] = pA;
      i++;   
   }   

	switch (m_Select_metric)
	{
	   case 0://RND
	   {
         for(pA = archive; pA != NULL; pA = pA->pNext)
         {   
	         pA->Z = 1.0;
         }
         break;
   	}
	   case 1://CD
	   {
         for(pA = archive; pA != NULL; pA = pA->pNext)
         {
	         pA->Z = 0.0;
         }
			for (int obj = 0; obj < archive->nF; obj++)
         {
				m_dim = obj;

				SortArchive(pSorted, archive_size, obj);

				for (int i = 1; i < archive_size - 1; i++)
            {
               double F0 = pSorted[0]->F[obj];
               double F1 = pSorted[i - 1]->F[obj];
               double F2 = pSorted[i + 1]->F[obj];
               double F4 = pSorted[archive_size - 1]->F[obj];
					pSorted[i]->Z += abs(F1 - F2) / abs(F0 - F4);
				}

            if(archive_size > 1)
				  pSorted[0]->Z = pSorted[1]->Z;
            if(archive_size > 2)
				  pSorted[archive_size - 1]->Z = pSorted[archive_size - 2]->Z;
			}/* end for(each objective function) */
		   break;
   	}/* end case() */
	   case 2://HVC_ESTIMATE
	   {
			for (int i = 0; i < archive_size; i++)
         {
			   pSorted[i]->Z = 0;
			}

			//define boundries
			double * f_low_bound = new double[m_num_objs];
			double * f_high_bound = new double[m_num_objs];

         //compute boundaries
			for (int i = 0; i < m_num_objs; i++)
         {
            m_dim = i;
				SortArchive(pSorted, archive_size, i);

				f_low_bound[i] = pSorted[0]->F[i];
				f_high_bound[i] = pSorted[archive_size - 1]->F[i];
			}

			//prepare monte carlo archive
			int dots_num = 100;
			ArchiveStruct ** mc_points = new ArchiveStruct*[dots_num];

			for (int i = 0; i < dots_num; i++)
         {
            //allocate sample
            ArchiveStruct * dot = new ArchiveStruct;
            dot->F = new double[m_num_objs];
            dot->nF = m_num_objs;
            dot->X = new double[m_num_dec];
            dot->nX = m_num_dec;
            dot->pNext = NULL;
            dot->Z = -999.999;
            
				for (int j = 0; j < m_num_objs; j++)
            {               					   
				   dot->F[j] = f_low_bound[j] + (f_high_bound[j] - f_low_bound[j])*UniformRandom();
				}

				mc_points[i] = dot;
			}/* end for() */

			for (int i = 0; i < dots_num; i++)
         {
		      int jj = archive_size; 
            bool any_good = false;

				for (int j = 0; j < archive_size; j++)
            {
			      if (dominion_status(pSorted[j], mc_points[i]) == 1)
               {
				      jj = j;
						any_good = true;
						break;
					}
				}

				if (any_good == true)
            {
			      for (int k = jj + 1; k < archive_size; k++)
               {
						if (dominion_status(pSorted[k], mc_points[i]) == 1)
                  {
							//this dot means nothing, goto the next dot.
							any_good = false;
							break;
						}
					}
				}

				//this dot is only dominated by archive[jj]
				if (any_good == true)
            {
					pSorted[jj]->Z += 1.00;
					any_good = false;
				}
			}/* end for() */

			//normalize z
         double best_z = 0.00;

			for (int i = 0; i < archive_size; i++)
         {
				pSorted[i]->Z = (pSorted[i]->Z / (double)dots_num);
				if (pSorted[i]->Z > best_z)
            {
					best_z = pSorted[i]->Z;
				}
         }/* end for() */

			for (int i = 0; i < archive_size; i++)
         {
				if (pSorted[i]->Z == 0)
            {
					pSorted[i]->Z = 0.5*best_z;
				}
			}

         delete [] f_low_bound;
		   delete [] f_high_bound;
         for(int i = 0; i < dots_num; i++)
         {
            DestroyArchive(mc_points[i]);
         }
         delete [] mc_points;
         break;
	   }/* end case() */
   	case 3://HVC_EXACT
	   {
			int dataNumber = archive_size;
			int dimension = m_num_objs;

			double * refPoint = new double[dimension];

   	   for(int i = 0; i < dimension; i++)
         {
				m_dim = i;
				SortArchive(pSorted, archive_size, i);
				refPoint[i] = 1.00001*(pSorted[archive_size - 1]->F[i]);
			}/* end for() */

			double ** pointsInitial = new double *[dataNumber];
			for (int i = 0; i < dataNumber; i++)
         {
            pointsInitial[i] = new double[dimension];

				for (int j = 0; j < m_num_objs; j++)
            {
			      pointsInitial[i][j] = pSorted[i]->F[j];
				}
			}

			double HyperVolume = HV(dataNumber, dimension, refPoint, pointsInitial);

			//destroy pointsInitial
			for (int j = 0; j < dataNumber; j++)
         {
			   delete[] pointsInitial[j];
			}
         delete[] pointsInitial;

			double best_z = 0;

			for(int i = 0; i < dataNumber; i++) 
         {
	         int * included_points = new int[dataNumber];
				for(int j = 0; j < dataNumber; j++)
            {
               included_points[j] = j;
               //leave out ith point
               if(j >= i) included_points[j]++;
				}

				double ** pointsInitial_sub = new double * [dataNumber - 1];
				for (int j = 0; j < dataNumber - 1; j++)
            {
				   pointsInitial_sub[j] = new double[dimension];

					for (int k = 0; k < m_num_objs; k++)
               {
					   pointsInitial_sub[j][k] = pSorted[included_points[j]]->F[k];
					}
				}/* end for() */

				pSorted[i]->Z = HyperVolume - HV(dataNumber - 1, dimension, refPoint, pointsInitial_sub);

            if(pSorted[i]->Z > best_z)
            {
               best_z = pSorted[i]->Z; 
            }

				//destroy pointsInitial_sub
				for (int j = 0; j < dataNumber - 1; j++)
            {
				   delete[] pointsInitial_sub[j];
				}
            delete[] pointsInitial_sub;
            delete [] included_points;
			}/* end for() */

			// taking care of the edges
			for (int i = 0; i < dimension; i++)
         {
				m_dim = i;
				SortArchive(pSorted, archive_size, i);

				pSorted[0]->Z = best_z;
				pSorted[archive_size - 1]->Z = best_z;
			}

         delete [] refPoint;
         break;
	   }/* end case() */
	}/* end switch() */

	delete [] pSorted;
}/* end Calc_Z() */

/******************************************************************************
SortPoints()

Sort the array based on the which index.
******************************************************************************/
void PADDS::SortPoints(double ** X, int size, int which)
{
   double F1, F2;
   double * pTmp;
   for(int i = 0; i < size; i++)
   {
      for(int j = (i+1); j < size; j++)
      {
         F1 = X[i][which];
         F2 = X[j][which];
         if(F2 < F1)
         {
            pTmp = X[i];
            X[i] = X[j];
            X[j] = pTmp;
         }
      }
   }/* end for() */
}/* end SortPoints() */

/******************************************************************************
bool_vec_to_ulong()

Convert bitwise array of bools into equivalent integer.

to_ulong() simply calculates bit[0]*2^0 + bit[1]*2^1 + bit[2]*2^2
******************************************************************************/
int PADDS::bool_vec_to_ulong(bool * pB, int size)
{
   int sum = 0;
   for(int i = 0; i < size; i++)
   {
      sum += (int)(pB[i]*pow((double)2, i));
   }/* end for() */
   return sum;
}/* end bool_vec_to_ulong() */

/******************************************************************************
ulong_to_bool_vec()

Convert integer into equivalent bitwise array of bools.
******************************************************************************/
void PADDS::ulong_to_bool_vec(int val, bool * pB, int size)
{
   int mask = 1;
   int tmp;
   for(int i = 0; i < size; i++)
   {
      tmp = (mask & val);
      if(tmp == 0)
        pB[i] = false;
      else
        pB[i] = true;
      mask = mask << 1;
   }/* end for() */
}/* end ulong_to_bool_vec() */

/******************************************************************************
SortArchive()

Sort the archive by objective function or by Z if whichObj == -1.
******************************************************************************/
void PADDS::SortArchive(ArchiveStruct ** pArch, int size, int whichObj)
{
   double F1, F2;
   ArchiveStruct * pTmp;
   for(int i = 0; i < size; i++)
   {
      for(int j = (i+1); j < size; j++)
      {
         if(whichObj == -1)
         {
            F1 = pArch[i]->Z;
            F2 = pArch[j]->Z;
         }
         else
         {
            F1 = pArch[i]->F[whichObj];
            F2 = pArch[j]->F[whichObj];
         }
         if(F2 < F1)
         {
            pTmp = pArch[i];
            pArch[i] = pArch[j];
            pArch[j] = pTmp;
         }
      }
   }/* end for() */
}/* end SortArchive() */

/******************************************************************************
DestroyArchive()

Free up memory of archive.
******************************************************************************/
void PADDS::DestroyArchive(ArchiveStruct * pArch)
{
   ArchiveStruct * pDel;
   for(ArchiveStruct * pCur = pArch; pCur != NULL;)
   {
     pDel = pCur;
     pCur = pCur->pNext;
     delete [] pDel->F;
     delete [] pDel->X;
     delete pDel; 
   }
}/* end DestroyArchive() */

/******************************************************************************
UpdateArchive()

Update the dominated and non-dominated archives with latest sample.
******************************************************************************/
int PADDS::UpdateArchive(double * pX, int nX, double * pF, int nF)
{
   int i;
   double Fcur, Ftst;
   ArchiveStruct * pArch, * pCur, * pPrev, * pNxt;
   bool bDominates, bIsDominated, bMarkForInsertion;
   pArch =  new ArchiveStruct;
   pArch->F = pF;
   pArch->X = pX;
   pArch->nX = nX;
   pArch->nF = nF;
   pArch->pNext = NULL;
   pArch->Z = -999.999;
   
   //first entry is always non-dominated
   if((m_NumDom == 0) && (m_NumNonDom == 0))
   {
      m_pDom = NULL;
      m_pNonDom = pArch;
      m_NumNonDom++;
      return ARCHIVE_NON_DOM;
   }

   //assume solution is non-dominated until we discover otherwise
   bMarkForInsertion = true;

   //compare against current list of non-dominated solutions
   ArchiveStruct * pDummy;
   for(pCur = m_pNonDom; pCur != NULL;)
   {
      //save next item since pCur->pNext may be changed during processing
      pDummy = pCur->pNext;

      //does new solution (Ftst) dominate the existing solution (Fcur)?
      bDominates = true;
      for(i = 0; i < pArch->nF; i++)
      {
         Fcur = pCur->F[i];
         Ftst = pArch->F[i];
         if(Fcur < Ftst)
         {
            bDominates = false;
            break;
         }/* end if() */
      }/* end for() */

      //is new solution (Ftst) dominated by an existing solution (Fcur)?
      if(bDominates == false)
      {
         bIsDominated = true;
         for(i = 0; i < pArch->nF; i++)
         {
            Fcur = pCur->F[i];
            Ftst = pArch->F[i];
            if(Ftst < Fcur)
            {
               bIsDominated = false;
               break;
            }/* end if() */
         }/* end for() */
      }/* end if() */

      /* -----------------------------------------------------------------------------
      Existing solution is dominated. Remove it from list of non-dominated solutions
      and mark new solution for insertion into the non-dominated list.
      ----------------------------------------------------------------------------- */
      if(bDominates == true)
      {
         //solution to be removed is at head of list.
         if(pCur == m_pNonDom)
         {
            m_pNonDom = pCur->pNext;
            m_NumNonDom--;
         }/* end if() */
         else //somewhere in middle or at end.
         {
            for(pPrev = m_pNonDom; pPrev->pNext != pCur;)
            {
               pPrev = pPrev->pNext;
            }
            pNxt = pCur->pNext;
            pPrev->pNext = pNxt; //removes item from list
            m_NumNonDom--;
         }/* end else() */

         //insert at head of dominated list
         pCur->pNext = NULL;
         pNxt = m_pDom;
         m_pDom = pCur;
         m_pDom->pNext = pNxt;
         m_NumDom++;         
      }/* end if() */
      /* -----------------------------------------------------------------------------
      New solution is dominated. Make note so that it is not inserted into the list
      of non-dominated solutions.
      ----------------------------------------------------------------------------- */
      else if(bIsDominated == true)
      {
         bMarkForInsertion = false;   
      }/* end if() */

      //advance to next item
      pCur = pDummy;
   }/* end for() */

   //insert new solution into list of non-dominated solutions?
   if(bMarkForInsertion == true)
   {
      //insert at head of non-dominated list
      pNxt = m_pNonDom;
      m_pNonDom = pArch;
      m_pNonDom->pNext = pNxt;
      m_NumNonDom++;               
      return ARCHIVE_NON_DOM;
   }/* end if() */
   else
   {
      //insert at head of dominated list
      pNxt = m_pDom;
      m_pDom = pArch;
      m_pDom->pNext = pNxt;
      m_NumDom++;               
      return ARCHIVE_DOM;
   }/* end else() */
}/* end UpdateArchive() */

/******************************************************************************
dominion_status()

Determine whether solution x1 dominates solution x2.
******************************************************************************/
int PADDS::dominion_status(ArchiveStruct * x1, ArchiveStruct * x2)
{
	int ds = 1;

	for (int i = 0; i < m_num_objs; i++)
   {
		if (x1->F[i] > x2->F[i])
      {
			goto try_2;
		}
	}
	return 1;

try_2:

	ds = 2;
	for (int i = 0; i < m_num_objs; i++)
   {
		if (x1->F[i] < x2->F[i])
      {
			return 0;
		}
	}

	return 2;
}/* end dominion_status() */

/******************************************************************************
SelectFrom()

Select and entry from the archive.
******************************************************************************/
ArchiveStruct * PADDS::SelectFrom(ArchiveStruct * pArchive)
{
   ArchiveStruct * pA;
   ArchiveStruct ** archive;
   int i, archive_size;

   //compute size of archive
   archive_size = 0;
   for(pA = pArchive; pA != NULL; pA = pA->pNext)
   {   
      archive_size++;
   }
   //make space for sorted archive
   archive = new ArchiveStruct *[archive_size];
   //initialize array
   i = 0;
   for(pA = pArchive; pA != NULL; pA = pA->pNext)
   {
      archive[i] = pA;
      i++;   
   }   

	double * z_cum = new double[archive_size];
	z_cum[0] = 0;
	for (int i = 0; i < archive_size; i++)
   {
		if (i == 0)
      {
			z_cum[i] = archive[i]->Z;
		}
		else
      {
			z_cum[i] = z_cum[i - 1] + archive[i]->Z;
		}
	}/* end for() */

	double t = UniformRandom() * z_cum[archive_size - 1];

	int ii = 0;

	for (int i = 0; i < archive_size; i++)
   {
		if (z_cum[i] >= t)
      {
			ii = i;
			break;
		}
	}

   delete [] z_cum;
   pA = archive[ii];
   delete [] archive;
	return pA;
}/* end SelectFrom() */

/******************************************************************************
neigh_value_continuous()
******************************************************************************/
double PADDS::neigh_value_continuous(double s, double s_min, double s_max, double r)
{
	double s_range = s_max - s_min;

	double snew = s + GaussRandom() * r * s_range;

	double P_Abs_or_Ref = UniformRandom();

	if (snew < s_min) 
   {
		if (P_Abs_or_Ref <= 0.5)
      {
			snew = s_min + (s_min - snew);
		}
		else
      {
			snew = s_min;
		}
		if (snew > s_max)
      { 
         snew = s_min; 
      }
	}
	else if (snew > s_max)
   {
		if (P_Abs_or_Ref <= 0.5)
      {
			snew = s_max - (snew - s_max);
		}
		else
      {
			snew = s_max;
		}

		if (snew < s_min)
      {
			snew = s_max;
		}
	}
	return snew;
}/* end neigh_value_continuous() */

/******************************************************************************
HV()

Hypervolume calculation.
******************************************************************************/
double PADDS::HV(int data_n, int dim_n, double * ref, double ** points)
{
	int i, j;
	int dataNumber = data_n;

   //set global
   m_dimension = dim_n;

	double ** pointsInitial = new double *[dataNumber];
	for (int n = 0; n < dataNumber; n++)
   {
		pointsInitial[n] = new double[dim_n];
		for (int i = 0; i < dim_n; i++)
      {
			pointsInitial[n][i] = points[n][i];
		}
	}
	double* refPoint = new double[dim_n];
	for (i = 0; i < dim_n; i++)
   {
		refPoint[i] = ref[i];
	}

	// initialize volume
	m_volume = 0.0;
	// sqrt of dataNumber
	m_dSqrtDataNumber = sqrt((double)dataNumber);

	// initialize region
	double* regionLow = new double[dim_n - 1];
	double* regionUp = new double[dim_n - 1];
	for (j = 0; j < dim_n - 1; j++) 
   {
		// determine minimal j coordinate
		double min = NEARLY_HUGE;
		for (i = 0; i < dataNumber; i++)
      {
			if (pointsInitial[i][j] < min) 
         {
				min = pointsInitial[i][j];
			}
		}
		regionLow[j] = min;
		regionUp[j] = refPoint[j];
	}

	// sort pointList according to d-th dimension
	SortPoints(pointsInitial, dataNumber, dim_n - 1);

	// call stream initially
	stream(regionLow, regionUp, pointsInitial, dataNumber, 0, refPoint[m_dimension - 1]);

   for (int n = 0; n < dataNumber; n++)
   {
		delete [] pointsInitial[n];
	}
   delete [] pointsInitial;
	delete [] refPoint;
	delete [] regionLow;
	delete [] regionUp;

	// print hypervolume
	return m_volume;
} /* end HV() */

/******************************************************************************
F()

The objective functions.
******************************************************************************/
void PADDS::F(ArchiveStruct * pA)
{
   m_pModel->GetParamGroupPtr()->WriteParams(pA->X);
   m_pModel->Execute(pA->F, pA->nF);
}/* end F() */

/******************************************************************************
covers()
******************************************************************************/
bool PADDS::covers(double* cub, double * regLow) 
{
	int i_hvc;
	for (i_hvc = 0; i_hvc < m_dimension - 1; i_hvc++) 
   {
		if (cub[i_hvc] > regLow[i_hvc]) 
      {
			return false;
		}
	}
	return true;
}/* end covers() */

/******************************************************************************
partCovers()
******************************************************************************/
bool PADDS::partCovers(double* cub, double * regUp)
{
	int i_hvc;
	for (i_hvc = 0; i_hvc < m_dimension - 1; i_hvc++)
	{
		if (cub[i_hvc] >= regUp[i_hvc])
		{
			return false;
		}
	}
	return true;
}/* end partCovers() */

/******************************************************************************
containsBoundary()
******************************************************************************/
int PADDS::containsBoundary(double* cub, double * regLow, int split) 
{
	// condition only checked for split>0
	if (regLow[split] >= cub[split])
   {
		// boundary in dimension split not contained in region, thus
		// boundary is no candidate for the splitting line
		return -1;
	}
	else 
   {
		int j_hvc;
		for (j_hvc = 0; j_hvc < split; j_hvc++) 
      { // check boundaries
			if (regLow[j_hvc] < cub[j_hvc]) 
         {
				// boundary contained in region
				return 1;
			}
		}
	}
	// no boundary contained in region
	return 0;
}/* end containsBoundary() */

/******************************************************************************
getMeasure()
******************************************************************************/
double PADDS::getMeasure(double * regLow, double * regUp) 
{
	double vol_hvc;
	int i_hvc;
	vol_hvc = 1.0;
	for (i_hvc = 0; i_hvc < m_dimension - 1; i_hvc++) 
   {
		vol_hvc *= (regUp[i_hvc] - regLow[i_hvc]);
	}
	return vol_hvc;
}/* end getMeasure() */

/******************************************************************************
isPile()
******************************************************************************/
int PADDS::isPile(double* cub, double * regLow, double * regUp) 
{
	int pile_hvc;
	int k_hvc;

	pile_hvc = m_dimension;
	// check all dimensions of the node
	for (k_hvc = 0; k_hvc < m_dimension - 1; k_hvc++) 
   {
		// k-boundary of the node's region contained in the cuboid? 
		if (cub[k_hvc] > regLow[k_hvc]) 
      {
			if (pile_hvc != m_dimension) 
         {
				// second dimension occured that is not completely covered
				// ==> cuboid is no pile
				return -1;
			}
			pile_hvc = k_hvc;
		}
	}
	// if pile == this.dimension then
	// cuboid completely covers region
	// case is not possible since covering cuboids have been removed before

	// region in only one dimenison not completly covered 
	// ==> cuboid is a pile 
	return pile_hvc;
}/* isPile() */

/******************************************************************************
computeTrellis()
******************************************************************************/
double PADDS::computeTrellis(double * regLow, double * regUp, double * trellis) 
{
	int i_hvc, j_hvc;
	double vol_hvc;
	int numberSummands_hvc;
	double summand_hvc;
	bool * bitvector_hvc = new bool[m_dimension];

	vol_hvc = 0.0;
	summand_hvc = 0.0;
	numberSummands_hvc = 0;

	// calculate number of summands
	bool * nSummands = new bool[m_dimension];
	for (i_hvc = 0; i_hvc < m_dimension - 1; i_hvc++) 
   {
		nSummands[i_hvc] = 1;
	}
	numberSummands_hvc = bool_vec_to_ulong(nSummands, m_dimension);

	double* valueTrellis = new double[m_dimension - 1];
	double* valueRegion = new double[m_dimension - 1];
	for (i_hvc = 0; i_hvc < m_dimension - 1; i_hvc++) 
   {
		valueTrellis[i_hvc] = trellis[i_hvc] - regUp[i_hvc];
	}
	for (i_hvc = 0; i_hvc < m_dimension - 1; i_hvc++) 
   {
		valueRegion[i_hvc] = regUp[i_hvc] - regLow[i_hvc];
	}

	double* dTemp = new double[numberSummands_hvc / 2 + 1];

	// sum
	for (i_hvc = 1; i_hvc <= numberSummands_hvc / 2; i_hvc++) 
   {

		// set bitvector length to fixed value 16
		// TODO Warning: dimension-1 <= 16 is assumed
		ulong_to_bool_vec(i_hvc, bitvector_hvc, m_dimension);

		// construct summand
		// 0: take factor from region
		// 1: take factor from cuboid
		summand_hvc = 1.0;
		for (j_hvc = 0; j_hvc < m_dimension - 2; j_hvc++) 
      {
			if (bitvector_hvc[j_hvc]) 
         {
				summand_hvc *= valueTrellis[j_hvc];
			}
			else 
         {
				summand_hvc *= valueRegion[j_hvc];
			}
		}
		summand_hvc *= valueRegion[m_dimension - 2];

		// determine sign of summand
		vol_hvc -= summand_hvc;
		dTemp[i_hvc] = -summand_hvc;

		// add summand to sum
		// sign = (int) pow((double)-1, (double)counterOnes+1);
		// vol += (sign * summand); 
	}

	ulong_to_bool_vec(i_hvc, bitvector_hvc, m_dimension);

	summand_hvc = 1.0;
	for (j_hvc = 0; j_hvc < m_dimension - 1; j_hvc++) 
   {
		if (bitvector_hvc[j_hvc]) 
      {
			summand_hvc *= valueTrellis[j_hvc];
		}
		else 
      {
			summand_hvc *= valueRegion[j_hvc];
		}
	}
	vol_hvc -= summand_hvc;

	for (i_hvc = 1; i_hvc <= numberSummands_hvc / 2; i_hvc++) 
   {
		summand_hvc = dTemp[i_hvc];
		summand_hvc *= regUp[m_dimension - 2] - trellis[m_dimension - 2];
		summand_hvc /= valueRegion[m_dimension - 2];
		vol_hvc -= summand_hvc;
	}
	//she has already set it{
	delete [] valueTrellis;
	delete [] valueRegion;
   delete [] dTemp;
   delete [] bitvector_hvc;
   delete [] nSummands;
	//she has already set it}

	return vol_hvc;
}/* end trellis() */

/******************************************************************************
getMedian()

return median of the list of boundaries considered as a set
******************************************************************************/
double PADDS::getMedian(double * bounds, int size) 
{
   double * pSort;
   double tmp, median;

	// do not filter duplicates	
	if (size == 1) 
   {
		return bounds[0];
	}
	else if (size == 2) 
   {
		return bounds[1];
	}

   pSort = new double[size];
   for(int i = 0; i < size; i++)
   {
      pSort[i] = bounds[i];
   }
   for(int i = 0; i < size; i++)
   {
      for(int j = i+1; j < size; j++)
      {
         if(pSort[j] < pSort[i])
         {
            tmp = pSort[j];
            pSort[j] = pSort[i];
            pSort[i] = tmp;
         }
      }
   }

   median = pSort[size/2];

   delete [] pSort;

   return median;
}/* end getMedian() */

/******************************************************************************
stream()

recursive calculation of hypervolume
******************************************************************************/
void PADDS::stream(double * regionLow, double * regionUp, double ** points, int npoints, 
            int split, double cover) 
{
	//--- init --------------------------------------------------------------//
	double coverOld_hvc;
	coverOld_hvc = cover;
	unsigned int coverIndex_hvc = 0;
	int c_hvc;

	//--- cover -------------------------------------------------------------//

	// identify first covering cuboid
	double dMeasure = getMeasure(regionLow, regionUp);
	while (cover == coverOld_hvc && coverIndex_hvc < (unsigned int)npoints)
   {
		if (covers(points[coverIndex_hvc], regionLow)) 
      {
			// new cover value
			cover = points[coverIndex_hvc][m_dimension - 1];
			m_volume += dMeasure * (coverOld_hvc - cover);
		}
		else coverIndex_hvc++;
	}/* end while() */

	/* coverIndex shall be the index of the first point in points which
	* is ignored in the remaining process
	*
	* It may occur that that some points in front of coverIndex have the same
	* d-th coordinate as the point at coverIndex. This points must be discarded
	* and therefore the following for-loop checks for this points and reduces
	* coverIndex if necessary.
	*/
	for (c_hvc = coverIndex_hvc; c_hvc > 0; c_hvc--) 
   {
		if (points[c_hvc - 1][m_dimension - 1] == cover) 
      {
			coverIndex_hvc--;
		}
	}/* end for() */

	// abort if points is empty
	if (coverIndex_hvc == 0) 
   {
		return;
	}
	// Note: in the remainder points is only considered to index coverIndex

	//--- allPiles  ---------------------------------------------------------//

	bool allPiles = true;
	unsigned int iii;

	int* piles_hvc = new int[coverIndex_hvc];
	for (iii = 0; iii < coverIndex_hvc; iii++) 
   {
		piles_hvc[iii] = isPile(points[iii], regionLow, regionUp);
		if (piles_hvc[iii] == -1) 
      {
			allPiles = false;

			//she has already set it{
			delete[] piles_hvc;
			//she has already set it{

			break;
		}
	}/* end for() */

	/*
	* tre llis[i] contains the values of the minimal i-coordinate of
	* the i-piles.
	* If there is no i-pile the default value is the upper bpund of the region.
	* The 1-dimensional KMP of the i-piles is: reg[1][i] - tre llis[i]
	*
	*/
	if (allPiles) 
   { // sweep

		// initialize trellis with region's upper bound
		double* trellis = new double[m_dimension - 1];
		for (c_hvc = 0; c_hvc < m_dimension - 1; c_hvc++) 
      {
			trellis[c_hvc] = regionUp[c_hvc];
		}

		double current = 0.0;
		double next = 0.0;
		iii = 0;
		do 
      { // while(next != coverNew)
			current = points[iii][m_dimension - 1];
			do 
         { // while(next == current)
				if (points[iii][piles_hvc[iii]] < trellis[piles_hvc[iii]]) 
            {
					trellis[piles_hvc[iii]] = points[iii][piles_hvc[iii]];
				}
				iii++; // index of next point
				if (iii < coverIndex_hvc) 
            {
					next = points[iii][m_dimension - 1];
				}
				else 
            {
					next = cover;
				}

			} while (next == current);
			m_volume += computeTrellis(regionLow, regionUp, trellis) * (next - current);
		} while (next != cover);
	}
	//--- split -------------------------------------------------------------//
	// inner node of partition tree
	else
   {
		double bound = -1.0;
      int boundaries_size = 0;
      int noBoundaries_size = 0;
		double * boundaries = new double[coverIndex_hvc];
		double * noBoundaries = new double[coverIndex_hvc];

		do 
      {
			for (iii = 0; iii < coverIndex_hvc; iii++) 
         {
				int contained = containsBoundary(points[iii], regionLow, split);
				if (contained == 1) 
            {
					boundaries[boundaries_size] = points[iii][split];
               boundaries_size++;
				}
				else if (contained == 0) 
            {
					noBoundaries[noBoundaries_size] = points[iii][split];
               noBoundaries_size++;
				}
			}

			if (boundaries_size >  0) 
         {
				bound = getMedian(boundaries, boundaries_size);
			}
			else if (noBoundaries_size > m_dSqrtDataNumber) 
         {
				bound = getMedian(noBoundaries, noBoundaries_size);
			}
			else 
         {
				split++;
			}
		} while (bound == -1.0);

		//destroy boundaries vector
		delete [] boundaries;
		delete [] noBoundaries;

		double dLast;
      int pointsChild_size = 0;
		double ** pointsChild = new double *[coverIndex_hvc];

		// left child
		// reduce maxPoint
		dLast = regionUp[split];
		regionUp[split] = bound;
		for (iii = 0; iii < coverIndex_hvc; iii++) 
      {
			if (partCovers(points[iii], regionUp)) 
         {
				pointsChild[pointsChild_size] = points[iii];
            pointsChild_size++;
			}
		}
		if (pointsChild_size > 0) 
      {
			stream(regionLow, regionUp, pointsChild, pointsChild_size, split, cover);
		}

		// right child
		// increase minPoint
		pointsChild_size = 0;

		regionUp[split] = dLast;
		dLast = regionLow[split];
		regionLow[split] = bound;
		for (iii = 0; iii < coverIndex_hvc; iii++) 
      {
			if(partCovers(points[iii], regionUp)) 
         {
				pointsChild[pointsChild_size] = points[iii];
            pointsChild_size++;
			}
		}
		if (pointsChild_size > 0) 
      {
			stream(regionLow, regionUp, pointsChild, pointsChild_size, split, cover);
		}
		regionLow[split] = dLast;

      delete [] pointsChild;
	}// end inner node
} /* end stream() */

/******************************************************************************
PADDS_Program()

Calibrate or optimize the model using PADDS.
******************************************************************************/
void PADDS_Program(int argC, StringType argV[])
{
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model();

   NEW_PRINT("PADDS", 1);
   PADDS * TheAlg = new PADDS(model);
   MEM_CHECK(TheAlg);

   if(model->GetObjFuncId() == OBJ_FUNC_WSSE) { TheAlg->Calibrate(); }
   else { TheAlg->Optimize(); }

   delete TheAlg;
   delete model;
} /* end PADDS_Program() */
