/******************************************************************************
File     : StatsClass.h
Author   : L. Shawn Matott
Date     : 06/12/2003
Version  : V01
Copyright: 2003, L. Shawn Matott

The StatsClass is used to compute statistical measures, following a successful
calibration. Many of the statistics require the Jacobian matrix, so this class
also provides routines for evaluation of the Jacobian, the transpose of the
Jacobian, and also the Normal matrix (JQJ).

The following statistical measures are avaiable:

   variance/standard deviation   
   covariance/standard error
   correlation coefficients
   Beale's Nonlinearity Measure
   Linssen's Nonlinearity Measure
   Cook's D Observation Influence Measure
   DFBETAS Observation Influence Measure
   Confidence Intervals
   Parameter Sensitivities
   Normal Probability Plot (R2N)

Version History
06-12-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
11-19-03    lsm   Parameter output uses WRITE_TX_BNR definition.
11-25-03    lsm   Added sensitivity measures, corrected FSTAT in Beale/Linssen
03-24-04    lsm   Added GetLinearity() to support PSO-LevMar hybrid
                  Added CalcOptimalStepSize() option
01-11-05    lsm   Added metrics
03-21-05    lsm   Added support for parameter-specific relative increments
01-01-07    lsm   StatsClass now uses abstract model base class (ModelABC). 
                  The Jacobian can now be calculated in parallel. Added support
                  for the detection of temporarily insensitive parameters and
                  observations; when detected such insensitive terms are "held",
                  i.e. removed from consideration. Added support for additional 
                  finte difference increments. Added Durbin-Watson, Runs Test 
                  and MMRI statistics.
******************************************************************************/
#ifndef STATS_CLASS_H
#define STATS_CLASS_H

#include "MyHeaderInc.h"

// forward decs
class ModelBackup;
class ResponseVarGroup;
class ModelABC;

typedef struct RUNS_STRUCT
{
   int pos, neg, runs;
   int clwr, cupr;
   bool bSuccess;
}RunsStruct;

typedef struct AUTORUN_STRUCT
{
   int sur, def, n1;
   double r1, var, vpx, med, clwr, cupr;
}AutorunStruct;

typedef struct MMRI_STRUCT
{
   double AIC, AICc, AICu, BIC, HQ;
   bool bSuccess;
}MMRI_Struct;

/******************************************************************************
class StatsClass

******************************************************************************/
class StatsClass
{
   public:
      int GetNumHeldObs(void){ return m_NumHeldObs;}
      int GetNumHeldParams(void){ return m_NumHeldParams;}
      Unchangeable2DArray CalcJacobian(bool bOkToHoldParams, bool bOkToHoldObs, double * pBestSavedF);
      Unchangeable2DArray CalcJacobian(double * pBestSavedF);
      Unchangeable2DArray GetJacobT(void);
      Unchangeable2DArray GetJacobUW(void);
      Unmoveable1DArray GetMinJac(void){ return m_pMinJac;}
      void AdjustJacobian(void);
      Unchangeable2DArray CalcNormal(void);      
      void AdjustVector(double * vec, bool obs);
      void CalcStats(void);
      void WriteStats(FILE * pFile);
      void PrintStats(void);
      void WriteMetrics(FILE * pFile);      

      StatsClass(ModelABC * pModel);
      ~StatsClass(void){ DBG_PRINT("StatsClass::DTOR"); Destroy(); }
      void Destroy(void);

      Unchangeable1DArray CalcResiduals(void);
      void AdjustResiduals(void);
      void WriteResiduals(int step, char * prefix);

   private:
      void BcastMinJac(void);
      void BcastJacobian(void);
      void EvalJacSerial(double * pBestSavedF);
      void EvalJacParallel(void);
      void EvalJacSuperMUSE(void);
      void InitFromFile(IroncladString pStatsFileName); 
      void CalcBealeAndLinssen(void);
      void CalcCooksD(void);
      void CalcDFBETAS(void);
      void CalcHatAndChangeMatrices(void);
      void CalcNormProbPlot(void);
      void CalcBestBoxCox(void);
      void CalcRawRy(void);
	   void CalcWeightedRy(void);
      void CalcCI(void);
      void CalcMMRI(bool bInv);
      void CalcSensitivities(void);
      void CalcPredictions(bool bStats, double ** pV, int np);
      double CalcOptimalStepSize(int idx, double * params);
      double AdjustObjFunc(double val);

      //Configuration variables (to be input by user)
      FiniteDiffType m_DiffType; //finite diff. method used to calc. derivatives
      double * m_pDiffInc; //finite diff. increment array
      FiniteDiffIncType m_DiffIncType; //type of FD increment
      double   m_MinInc; //smallest allowable increment
           
      //these flags are set to true if user wishes the given stat to be computed.      
      bool m_bNoStats; ////default = false
      bool m_StdDevFlag; //default = true
      bool m_StdErrFlag; //default = true
      bool m_CorrCoefFlag; //default = true
      bool m_NormPlotFlag; //default = false
      bool m_BealeFlag; //default = false
      bool m_LinssenFlag; //default = false
      bool m_CooksFlag; //default = false
      bool m_DfbetasFlag; //default = false       
      bool m_MatricesFlag; //default = false, if true write J, JQJ, Inv(JQJ), and residuals
      bool m_CIflag;      //confidence interval on each parameter
      bool m_SensFlag;    //parameter sensitivities
      bool m_RunsTestFlag; //default = false
      bool m_MMRI_Flag; //default = false
      bool m_bInv; //no default, indicates whehter or not normal matrix could be inverted
      bool m_bDOF; //indicates whether sufficient dof for AICc, AICu metrics
      bool m_bOkToHoldParams;
      bool m_bOkToHoldObs;
      bool m_BestBoxCoxFlag; //if true, determine optimal lambda value for BoxCox transformation
      bool m_AutorunFunctionFlag; // a third autocorrelation test, breaks ties between Runs and Durbin-Watson test
      bool m_bWriteIterationResiduals; //if true, write residuals file at each iteration

      //need three model backups for FD computations
      ModelBackup * m_pMidBkup;
      ModelBackup * m_pLowBkup;
      ModelBackup * m_pHiBkup;
      
      double * m_pOrdResid; //ordered (after weighting) observation residuals
      double * m_pExpResid; //expected value of the std. normal ordered residuals      
      double * m_pResid;    //observation residuals
      double * m_pMinJac;   //best (min) configuration discovered during a given Jacobian calculation
      double ** m_pJacob;   //Jacobian matrix (J)
      double ** m_pJacobUW;  //Unweighted Jacobian matrix (J)
      double ** m_pJacobT;  //Transpose of Jacobian matrix (J)
      double ** m_pJacPred; //partial derivatives of predictions
      double ** m_pPbyO1;   //Temp. matrix, N = num. params, m = num. obs.
      double ** m_pNormal;  //normal matrix (JQJ)
      double ** m_pInvNormal;  //Inverse of normal normal matrix (JQJ)^-1
      double ** m_pHat;     //Hat matrix, (Q^1/2)J[(JQJ)^-1](Q^1/2)J
      double ** m_pChange;  //Change matrix, [(JQJ)^-1](Q^1/2)J
      double ** m_pCovar;   //covariance matrix: var*(JQJ)^-1
      double *  m_pCooksD;  //Cook's D for each observation
      double ** m_pDFBETAS; //DFBETAS matrix (n by p)
      double ** m_pScaledSens; //scaled sensitivity matrix
      double *  m_pCompScaledSens; //composite scaled sensitivity matrix
      double ** m_pPctScaledSens; //one percent scaled sensitivity matrix

      double    m_Variance; //variance: phi/(num. obs. - num. params.)

      double    m_BealeStat;
      double    m_LinssenStat;
      double    m_BestBoxCoxVal;
      double    m_NonLinThresh; //Beale threshold for non-linearity
      double    m_EffLinThresh; //Beale threshold for effecively linear
      double    m_CooksInfluThresh;  //Cook's D 50% influence threshold
      int       m_NumInfluCooks;
      double    m_CooksAvgLvg;  //Cook's D average leverage
      int       m_NumInfluLvg;
      double    m_DfbetaInfluThresh; //DFBETAS influence threshold
      int       m_NumInfluDfbeta;
      double    m_OrdCorrCoeff; //correlation of residual to normalal dist.
      double    m_WeightedRy; //correlation of weighted measured and weighted simulated values
      double	 m_RawRy; //correlation of measured and simulated values (non-weighted)
      double  * m_pCIlwr;     //lower confidence intervals
      double  * m_pCIupr;     //upper confidence intervals
      double    m_CIpct;      //confidence interval percentage

      /*The joint confidence ellipsoid percentage that will create a region 
      with volume equal to the volume of the CI block specified by the 
      individual parameter confidence intervals, IF all parameters were 
      uncorrelated. */
      double    m_EllipsePct; 
   
      RunsStruct m_Runs;
      AutorunStruct m_AR;
      MMRI_Struct m_MMRI;
      
      bool m_bAdjustedJac; //true if Jacobian has been adjusted
      int m_NumParams;  //number of parameters
      int m_NumHeldParams; //number of insensitive parameters
      bool * m_bHoldParam;
      int m_NumObs;     //number of observations
      bool * m_bHoldObs;
      int m_NumHeldObs; //number of insensitive observations

      ModelABC * m_pModel;

      //Finite difference matrices
      double ** m_ParaMat; //3x3 matrix used in parabolic solution
      double ** m_ParaInv; //3x3 inverse of paraMat

      //arrays used in parallel Jacobian calculation 
      double * m_pBuf;
      double * m_pMyBuf;

      //storage for Jacobian vars, used in SuperMUSE evaluation
      FiniteDiffType * m_pDType;
      double * m_pDx;
      double * m_pMid;
      double * m_pHi;
      double * m_pLow;

      double m_Phi;

      //predictions
      ResponseVarGroup * m_pPredictions;
      double * m_Pred;
      double * m_PredSD;
      double * m_PredCIlwr;
      double * m_PredCIupr;

      //metrics
      int m_DiffCount;
      int m_StepCount;
      int m_StatsCount;
}; /* end class StatsClass */

extern "C" {
void STATS_Program(int argc, StringType argv[]);
void Jacobian_Program(int argc, StringType argv[]);
void EVAL_Program(int argc, StringType argv[]);
int ResumeEvaluations(ModelABC * pModel, int id, int nprocs, double * pList);
}
#endif /* STATS_CLASS_H */

