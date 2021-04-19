/******************************************************************************

Abstract Base Class for algorithms

******************************************************************************/
#ifndef MODEL_ABC_H
#define MODEL_ABC_H

#include "MyHeaderInc.h"
#include <vector>

//forward decs
class ObservationGroup;
class ObjectiveFunction;
class ParameterGroup;

/******************************************************************************
class AlgorthmABC2

Abstract base class for algorithm (main [complex] model and surrogate models)
******************************************************************************/
class AlgorithmABC2 {
   public:
      virtual ~AlgorithmABC2(void){ DBG_PRINT("AlgorthmABC2::DTOR"); }
      virtual void Destroy(void)=0;

      // Original members
      virtual void Destroy(void);
      virtual double** CreateInitialSample(int sampleSize);
      virtual double** CreateSample(int sampleSize);
      virtual void Optimize(void);
      virtual void Calibrate(void);
      virtual void WriteMetrics(FILE* pFile);
      virtual void WarmStart(void);
      virtual int  GetCurrentIteration(void);

      // Model members
      //virtual ObservationGroup * GetObsGroupPtr(void) = 0;
      //virtual ParameterGroup *  GetParamGroupPtr(void) = 0;
      //virtual ObjectiveFunction * GetObjFuncPtr(void) = 0;
      //virtual double GetObjFuncVal(void) = 0;
      //virtual void SetObjFuncVal(double curVal) = 0;
      //virtual int GetCounter(void) = 0;
      //virtual ObjFuncType GetObjFuncId(void) = 0;
      //virtual UnchangeableString GetObjFuncStr(void) = 0;
      //virtual UnchangeableString GetModelStr(void) = 0;
      virtual void SaveBest(int id) = 0;
      virtual void WriteMetrics(FILE * pFile);
      //virtual int GetNumDigitsOfPrecision(void) = 0;
      //virtual TelescopeType GetTelescopingStrategy(void) = 0;
      virtual void PerformParameterCorrections(void) = 0;
      //virtual bool CheckWarmStart(void) = 0;


private:
    // Communication functions
    virtual void ConfigureWorkerDirectory(int workerRank);
    virtual void ConfigureWorkerExtraFiles(int workerRank);
    virtual void ConfigureWorkerFilePairs(int workerRank);
    virtual void ConfigureWorkerObservations(int workerRank);
    virtual void ConfigureWorkerParameterGroups(int workerRank);
};

#endif 
