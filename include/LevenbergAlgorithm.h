/******************************************************************************
File     : LevenbergAlgorithm.h
Author   : L. Shawn Matott
Copyright: 2003, L. Shawn Matott

The Levenberg-Marquardt algorithm is a hybrid numerical optimization method 
that initially uses the Steepest-Descent technique. However, since it is 
known that the Steepest-Descent algorithm converges very slowly near the 
optimum point, it is desirable to smoothly transition to a polynomial 
approximation method near the optimum. This implementation of the Levenberg 
algorithm is based upon the description/solution provided in the WinPEST 
user's manual, pages 9-42.

Version History
03-20-03    lsm   added copyright information and initial comments.
08-26-03    lsm   created version history field and updated comments.
08-17-04    lsm   RAM fragmentation fixes, metrics collection and reporting
                  Added support for user-requested program termination
01-11-05    lsm   Added algorithm metrics
01-01-07    lsm   Algorithm now uses abstract model base class (ModelABC).
******************************************************************************/
#ifndef LEVENBERG_ALGORITHM_H
#define LEVENBERG_ALGORITHM_H

#include "MyHeaderInc.h"

//parent class
#include "AlgorithmABC.h"

//forward decs
class ModelABC;
class ModelBackup;
class StatsClass;

/******************************************************************************
class LevenbergAlgorithm

******************************************************************************/
class LevenbergAlgorithm : public AlgorithmABC
{
   public:
      LevenbergAlgorithm(ModelABC * pModel, bool bMulti);
      ~LevenbergAlgorithm(void){ DBG_PRINT("LevenbergAlgorithm::DTOR"); Destroy(); }
      void Destroy(void);            
      void Optimize(void);
      void Calibrate(void);
      void WriteMetrics(FILE * pFile);
      void WarmStart(void);
      int  GetCurrentIteration(void) { return m_CurIter; }

  private:
      void InitFromFile(IroncladString pLevFileName); 
      void CalibrateGML(void);
      void AdjustLambda(void);
      double TryLambda(double lambda);
      void CalcJacobian(void);
      void CalcNormal(void);      
      void CalcScale(void);      
      void CalcAlpha(double lambda);
      void CalcUpgrade(void);
      void CalcGamma(void);
      void CalcBeta(void);
      void AdjModelParams(void);      

      //GML-MS routines
      void   GetRndParamSet(MyPoint *Point);
      double GetMinDist(MyPoint * Point, ParameterList * List);
      void CopyPoint(MyPoint * from, MyPoint * to, int np);
      void InsertParamSet(void);

      //Configuration variables (to be input by user)
      double m_Lambda;    //Marquardt lambda
      double m_LamSF;     //Marquardt lamdba scale factor
      double m_Converge;  //obj.func. convergence value
      double m_RatioConv;  //lambda convergence val. associated with m_PhiRatio
      double m_RelRedConv; //lambda convergence val. associated with m_PhiRelRed
      int m_MaxLambdas; //maximum number of lambdas to try      
      int m_MaxIter;    //maximum number of iterations  
      int m_CurIter;    
      double m_MoveLimit; //max. param. adj. (fraction of overall range)
      bool m_bMS; //true if using GML-MS
      int m_NumMS; //number of multi-starts (for GML-MS)
      ParameterList * m_pList;
  
      double m_Alpha;     //Marquardt parameter
      double m_Beta;      //SF for the optimal length of the upgrade vector      
      double m_Phi;       //objective function value
      double m_BestSavedPhi;
      double m_PhiRatio;  //ratio of old and cur obj. function vals
      double m_PhiRelRed; //relative reduction in old and cur obj. func. vals      
      
      //need four model backups for lamdba trials
      //m_pInitBkup = original model (before lambda trials)
      //m_pNonBkup  = model after non-ajusted lambda
      //m_pDecBkup  = model after decreased lambda
      //m_pIncBkup  = model after increased lambda
      ModelBackup * m_pInitBkup;
      ModelBackup * m_pNonBkup;
      ModelBackup * m_pDecBkup;
      ModelBackup * m_pIncBkup;
                  
      Unchangeable1DArray m_pResid; //vector of residuals (r)
      double * m_pUpgrade;  //upgrade vector (u)
      double * m_pTmpVec;   //temp. upgrade/parameter vector
      double * m_pGamma;    //gamma = Ju
      double ** m_pPbyP1;   //Temp. matrix, N = num. params
      double ** m_pPbyP2;   //Temp. matrix, N = num. params
      double ** m_pPbyO1;   //Temp. matrix, N = num. params, m = num. obs.
      double ** m_pPbyO2;   //Temp. matrix, N = num. params, m = num. obs.
      Unchangeable2DArray m_pJacob;   //Jacobian matrix (J)
      Unchangeable2DArray m_pJacobT;  //Transpose of Jacobian matrix (J)
      Unchangeable2DArray m_pJacobUW;  //Transpose of Jacobian matrix (J)
      double ** m_pScale;   //scaling matrix (S)      
      Unchangeable2DArray m_pNormal;  //normal matrix (JQJ)

      int m_NumParams;  //number of parameters
      int m_NumObs;     //number of observations

      StatsClass * m_pStats; //statistics class

      ModelABC * m_pModel;

      //alg. metrics
      int m_NumIters;
      int m_NumEvals;
      int m_NumUprViols;
      int m_NumLwrViols;
      int m_NumMoveViols; //move limit violations
}; /* end class LevenbergAlgorithm */

extern "C" {
void LEV_Program(int argc, StringType argv[]);
void GMLMS_Program(int argc, StringType argv[]);
}

#endif /* LEVENBERG_ALGORITHM_H */

