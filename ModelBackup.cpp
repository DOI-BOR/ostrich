/******************************************************************************
File     : ModelBackup.cpp
Author   : L. Shawn Matott and Vijaykumar Raghavan
Copyright: 2003, L. Shawn Matott and Vijaykumar Raghavan

The ModelBackup class stores and restores snapshots of the Parameter and 
Observation Groups of the Model. This class is convenient for algorithms such
as simulated annealing, which must make several trial moves from the same 
starting point before settling on the best move.  Additionally, the backup 
class is useful for finite difference computations, which must perterb 
model parameters without adversely affecting the overall optimization process.

Version History
03-09-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
08-13-04    lsm   distinguished between semi and full restore functionality.
01-01-07    lsm   ModelBackup now uses an abstract model base class (ModelABC).
******************************************************************************/
#include "ModelBackup.h"
#include "ModelABC.h"
#include "ResponseVarGroup.h"
#include "RespVarABC.h"
#include "ObservationGroup.h"
#include "Observation.h"
#include "ParameterGroup.h"
#include "ObjectiveFunction.h"

#include "Exception.h"

/******************************************************************************
CTOR

Size the parameter and observation storage arrays using the Observation and 
Parameter Groups of the Model and the response variable group pRV.
******************************************************************************/
ModelBackup::ModelBackup(ModelABC * pModel)
{
   ParameterGroup * pParamGroup;
   ObservationGroup * pObsGroup;

   m_pModel = pModel;
   
   pObsGroup = m_pModel->GetObsGroupPtr();     
   m_pObs = NULL;
   if(pObsGroup != NULL)
   { 
      m_NumObs = pObsGroup->GetNumObs();

      NEW_PRINT("double", m_NumObs);
      m_pObs  = new double[m_NumObs];
      MEM_CHECK(m_pObs);
   }/* end if() */

   pParamGroup = m_pModel->GetParamGroupPtr();
   m_NumParams = pParamGroup->GetNumParams();

   NEW_PRINT("double", m_NumParams);
   m_pParams = new double[m_NumParams];
   MEM_CHECK(m_pParams);

   m_NumPred = 0;
   m_pRV = NULL;
   m_pPred = NULL;
 
   IncCtorCount();
} /* end CTOR */

/******************************************************************************
SetResponseVarGroup()

Size the predictions storage arrays using the response variable group pRV. Also
save the pointer to the response variable group for future reference.
******************************************************************************/
void ModelBackup::SetResponseVarGroup(ResponseVarGroup * pRV)
{
   delete [] m_pPred;

   m_pRV = pRV;
   if(m_pRV != NULL)
   {
      m_NumPred = m_pRV->GetNumRespVars();
      NEW_PRINT("double", m_NumPred);
      m_pPred  = new double[m_NumPred];
      MEM_CHECK(m_pPred);
   }
}/* SetResponseVarGroup() */

/******************************************************************************
Destroy()

Free up the observation and parameter storage arrays.
******************************************************************************/
void ModelBackup::Destroy(void)
{  
   delete [] m_pParams;
   delete [] m_pObs;
   delete [] m_pPred;

   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
StoreParamVals()

Copy parameter group values into the storage array.
******************************************************************************/
void ModelBackup::StoreParamVals(void)
{  
   ParameterGroup * pGroup;

   pGroup = m_pModel->GetParamGroupPtr();
   pGroup->ReadParams(m_pParams);
} /* end StoreParamVals() */

/******************************************************************************
RestoreParamVals()

Copy stored parameters into the model parameter group.
******************************************************************************/
void ModelBackup::RestoreParamVals(void)
{
   ParameterGroup * pGroup;

   pGroup = m_pModel->GetParamGroupPtr();
   pGroup->WriteParams(m_pParams);
}/* end RestoreParamVals() */

/******************************************************************************
StoreObsVals()

Copy model computed observation group values into the storage array.
******************************************************************************/
void ModelBackup::StoreObsVals(void)
{ 
   int i;
 
   for(i = 0; i < m_NumObs; i++)
   {
      m_pObs[i] = m_pModel->GetObsGroupPtr()->GetObsPtr(i)->GetComputedVal(false, false);
   }
}/* end StoreObsVals() */

/******************************************************************************
RestoreObsVals()

Copy stored model computed observations into the model observation group.
******************************************************************************/
void ModelBackup::RestoreObsVals(void)
{
   ObservationGroup * pGroup;
   Observation * pObs;
   int i;
   double val;

   pGroup = m_pModel->GetObsGroupPtr();

   for(i = 0; i < m_NumObs; i++)
   {
      pObs = pGroup->GetObsPtr(i);
      val = m_pObs[i];
      pObs->SetComputedVal(val);
   }/* end for() */
}/* end RestoreObsVals() */

/******************************************************************************
StorePredictedVals()

Copy model computed predictions into the storage array.
******************************************************************************/
void ModelBackup::StorePredictedVals(void)
{
   if(m_NumPred <= 0) return;

   m_pRV->ExtractVals();

   int i;
  
   for(i = 0; i < m_NumPred; i++)
   {
      m_pPred[i] = m_pRV->GetRespVarPtr(i)->GetCurrentVal();
   }
}/* end StoreObsVals() */

/******************************************************************************
RestorePredictedVals()

Copy stored model computed predictions into the response variable group.
******************************************************************************/
void ModelBackup::RestorePredictedVals(void)
{
   if(m_NumPred <= 0) return;

   int i;

   for(i = 0; i < m_NumPred; i++)
   {
      m_pRV->GetRespVarPtr(i)->SetCurrentVal(m_pPred[i]);
   }/* end for() */
}/* end RestorePredictedVals() */

/******************************************************************************
Store()

Copy model parameter group and computed observation group values into the 
storage arrays.
******************************************************************************/
void ModelBackup::Store(void)
{ 
   StoreParamVals();
   if(m_pObs != NULL){ StoreObsVals();}
   StorePredictedVals();
   m_ObjFuncVal = m_pModel->GetObjFuncVal();
}/* end Store() */

/******************************************************************************
SemiRestore()

Copy stored model computed observations into the model observation group and 
copy stored paramters into the parameter group, also restore the objective 
function value. 

NOTE: It's a semi-restore becasue the model is not rerun using the restored 
parameters (therefore, tied parameters, response variables, and constraints 
will not be properly restored).
******************************************************************************/
void ModelBackup::SemiRestore(void)
{ 
   RestoreParamVals();
   if(m_pObs != NULL){ RestoreObsVals();}
   RestorePredictedVals();
   m_pModel->SetObjFuncVal(m_ObjFuncVal);
}/* end SemiRestore() */

/******************************************************************************
FullRestore()

Copy stored paramters into the parameter group and rerun the model. This will
enusure that all parts of the model are consistent with the restored set of 
parameters.
******************************************************************************/
void ModelBackup::FullRestore(void)
{ 
   RestoreParamVals();
   m_pModel->Execute();
   if(m_pRV != NULL){ m_pRV->ExtractVals();}
}/* end FullRestore() */

/******************************************************************************
FullRestore()

Copy stored paramters into the parameter group and rerun the model. This will
enusure that all parts of the model are consistent with the restored set of 
parameters.
******************************************************************************/
double ModelBackup::GetObs(int i, bool bTransformed, bool bWeighted)
{ 
   double y, w;
   y = m_pObs[i];
   w = GetObsWeight(m_pModel->GetObsGroupPtr()->GetObsPtr(i));

   if(bTransformed == true) //transformed implies also weighted
   {
     y = BoxCox(y*w);
   }
   else if(bWeighted == true)
   {
     y = (y*w);
   }

   return (y);
} /* end GetObs() */


