/******************************************************************************
File     : ModelBackup.h
Author   : L. Shawn Matott and Vijaykumar Raghavan
Copyright: 2003, L. Shawn Matott and Vijaykumar Raghavan

The ModelBackup class stores and restores snapshots of the Parameter and 
Observation Groups of the Model, along with the objective function value. This 
class is convenient for algorithms such as simulated annealing, which must make 
several trial moves from the same starting point before settling on the best 
move.  Additionally, the backup class is useful for finite difference 
computations, which must perterb model parameters without adversely affecting 
the overall optimization process.

Version History
03-09-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
08-13-04    lsm   distinguished between semi and full restore functionality.
01-01-07    lsm   ModelBackup now uses an abstract model base class (ModelABC).
******************************************************************************/
#ifndef MODEL_BACKUP_H
#define MODEL_BACKUP_H

#include "MyHeaderInc.h"

// forward decs
class ModelABC;
class ResponseVarGroup;

/******************************************************************************
class ModelBackup
******************************************************************************/
class ModelBackup
{
   public :
      ModelBackup(ModelABC * pModel);
      ~ModelBackup(void){ DBG_PRINT("ModelBackup::DTOR"); Destroy(); }
      void Destroy(void);

      void Store(void);
      void SemiRestore(void);
      void FullRestore(void);
      void SetResponseVarGroup(ResponseVarGroup * pRV);

      double GetParam(int i){ return m_pParams[i];} //parameters
      double GetObs(int i, bool bTransformed, bool bWeighted); //observations
      double GetPred(int i){ return m_pPred[i];} //predictions

   private:
      void StoreParamVals(void);
      void StoreObsVals(void);
      void RestoreParamVals(void);
      void RestoreObsVals(void);
      void StorePredictedVals(void);
      void RestorePredictedVals(void);

      ModelABC * m_pModel;
      double * m_pParams;
      int m_NumParams;
      double * m_pObs;
      int m_NumObs;
      ResponseVarGroup * m_pRV;
      double * m_pPred;
      int m_NumPred;
      double m_ObjFuncVal;
}; /* end class ModelBackup */

#endif /* MODEL_BACKUP_H */


