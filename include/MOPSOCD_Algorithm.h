/******************************************************************************
File     : MOPSOCD.h

Authors  : L. Shawn Matott, C.R. Raquel and P.C. Naval, Jr.

Copyrights: 2005 C.R. Raquel and P.C. Naval, Jr., University of the Philippines-Diliman
               --- original implementation as standalone C program
               --- MOPSO-CD Version 0.5b Feb. 17, 2006
            2017, L. Shawn Matott 
               --- adaptation into OSTRICH toolkit

Multi-Objective Particle Swarm Optimization with Crowding Distance (MOPSOCD).

This is an implementation of MOPSO-CD,a multiobjective particle swarm 
optimization algorithm using crowding distance.

For details please see:
 C.R. Raquel and P.C. Naval, Jr., An Effective Use of Crowding         
    Distance in Multiobjective Particle Swarm Optimization.            
    In Proc. of Genetic and Evolutionary Computation Conference        
    (GECCO 2005), Washington DC, June 2005.
                                                                       
E-mail address       : cvmig@engg.upd.edu.ph
Version              : mopsocd05b
Last updated         : Fri Feb 17, 2006
                                                                      
Random Generator Source code has been taken from Random Library found 
at http://www.swin.edu.au/astronomy/pbourke/software/random/          
                                                                      
Permission to use MOPSO-CD codes is hereby granted for academic and   
research purposes ONLY. Commercial usage of these codes is prohibited 
without prior knowledge of the authors.  In no way will the authors   
be held responsible for any possible faulty operation of              
software/hardware arising from the use of these codes.                

Version History
12-28-17    lsm   added copyright information and initial comments.
******************************************************************************/

#ifndef MOPSOCD_H
#define MOPSOCD_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "MOPSOCD_RandomLib.h" 

#define MOPSOCD_OPTIMIZATION_TYPE (0) /* optimization type, 0 for min, 1 for max */
#define MOPSOCD_VERBOSE (0)           /* verbosity level 0,1 */
#define MOPSOCD_PRINT_EVERY (10)      /* how frequently should output be generated */

// parent class
#include "AlgorithmABC.h"

// forward decs
class ModelABC;

/******************************************************************************
class MOPSOCD
Multi-Objective Particle Swarm Optimization with Crowding Distance
******************************************************************************/
class MOPSOCD : public AlgorithmABC
{
public:
   MOPSOCD(ModelABC * pModel);
   ~MOPSOCD(void){ DBG_PRINT("MOPSOCD::DTOR"); Destroy(); }
   void InitFromFile(IroncladString pFileName);
   void Optimize(void);
   void Calibrate(void);
   void Destroy(void);
   void WriteMetrics(FILE * pFile);
   void WarmStart(void){ return; }
   int  GetCurrentIteration(void) { return m_CurIter; }

private:
   int UpdateLists(double * pX, int nX, double * pF, int nF);
   void WriteLicenseInfo(void);

   void initialize_rand(void);
   void initialize_pop(void);
   void initialize_vel(void);
   void evaluate(void);
   void evaluate_parallel(int myrank, int nprocs);
   void broadcast_population(int myrank, int nprocs);
   void gather_results(int myrank, int nprocs);
   void store_pbests(void);
   void insert_nondom(void);
   void delete_particle(unsigned int k);
   unsigned int check_constraints(double * consVar);
   void crowding(void);
   void qsortFitness(unsigned int f, unsigned int begin, unsigned int lastPart);
   void qsortCrowd(unsigned int begin, unsigned int lastPart);
   void compute_distance(unsigned int f);
   void mutate(unsigned int t);
   void get_ranges(double *minvalue, double *maxvalue);
   void maintain_particles(void);
   void compute_velocity(void);
   void update_archive(void);
   unsigned int check_nondom(unsigned int i);
   void update_pbests(void);
   void save_results(char * archiveName);

   ModelABC * m_pModel;
   ArchiveStruct * m_pNonDomList; //complete list non-dominated solutions
   ArchiveStruct * m_pDomList; //complete list of dominated solutions
   int m_NonDomListSize;
   int m_DomListSize;
   int m_PopSize;     /* number of particles in the population */
   int m_MaxGen;      /* maximum number of generations         */
   int m_CurIter;
   int m_ArchiveSize; /* capacity of archive                   */
   int m_NumFun;      /* number of objective functions         */
   int m_NumVar;      /* number of variables                   */ 
   double m_PI;

   double ** m_ArchiveVar; /* [m_ArchiveSize][m_NumVar]; variable values of particles in the archive  */
   double ** m_ArchiveFit; /* [m_ArchiveSize][m_NumFun]; fitness values of particles in the archive   */
   double ** m_PopVar;     /* [m_PopSize][m_NumVar]; variable values of particles in the population   */
   double ** m_PopFit;     /* [m_PopSize][m_NumFun]; fitness values of particles in the population    */
   double ** m_PbestsVar;  /* [m_PopSize][m_NumVar]; personal bests of particles in the population    */
   double ** m_PbestsFit;  /* [m_popsize][m_NumFun]; personal bests of particles in the population    */
   double ** m_Velocity;   /* [m_popsize][m_NumVar]; velocity of particles in the population          */
   double *  m_CrowdDist;  /* [m_ArchiveSize];		  crowding distance values of particles in archive */
   double m_ProbMut;       /* probability of mutation                          */
   int m_NumNonDom;        /* number of non-dominated solutions in the archive */
}; /* end class MOPSOCD */

extern "C" {
   void MOPSOCD_Program(int argC, StringType argV[]);
}

#endif /* MOPSOCD_H */


