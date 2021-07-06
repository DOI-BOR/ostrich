/******************************************************************************

The Model class encapsulates the interaction of the Ostrich optimization tools
with the externally executed groundwater modeling program. The class divides
model components into three groups: the parameter group, the observation group
and the objective function group. In addition to being able to execute the
groundwater model, the Model class provides Ostrich algorithms with access to
these groups.

******************************************************************************/
// TODO: Doc string

#ifndef ALGORITHM_H
#define ALGORITHM_H

#include "MyHeaderInc.h"
#include "MyTypes.h"
#include <mpi.h>
#include <filesystem>
#include "WriteUtility2.h"
#include <cstring>
#include <vector>


//forward decs
class ObservationGroup;
class ParameterGroup;
class ObjectiveFunction;
class DecisionModule;
class SuperMUSE;
class FilePair;
class FileList;
class DatabaseABC;
class SurrogateParameterGroup;
class ParameterCorrection;



/******************************************************************************
class Algorthm

******************************************************************************/
class Algorithm {
public:
    Algorithm(void);
    ~Algorithm(void) { DBG_PRINT("AlgorithmABC2::DTOR"); Destroy(); }
    void Destroy(void);

    //retrieve member variables
    ObservationGroup* GetObsGroupPtr(void);
    ParameterGroup* GetParamGroupPtr(void);
    ObjectiveFunction* GetObjFuncPtr(void);
    //double GetObjFuncVal(void) { return m_CurObjFuncVal; }
    //void SetObjFuncVal(double curVal) { m_CurObjFuncVal = curVal; }
    int GetCounter(void);
    void SetCounter(int count);
    ObjFuncType GetObjFuncId(void) { return m_ObjFuncId; }
    UnchangeableString GetObjFuncStr(void);
    UnchangeableString GetModelStr(void) { return m_ExecCmd; }
    void PerformParameterCorrections(void);
    //misc. member functions     
    void   CheckGlobalSensitivity(void);
    
    //void   WriteMetrics(FILE* pFile);
    //void   Bookkeep(bool bFinal);
    //int GetNumDigitsOfPrecision(void) { return m_Precision; }
    //bool CheckWarmStart(void) { return m_bWarmStart; }
    
    //FilePair* GetFilePairs(void) { return m_FileList; }
    //TelescopeType GetTelescopingStrategy(void) { return m_Telescope; }

    double m_BestObjective = INFINITY;                          // Best objective
    std::vector<double> m_BestAlternative;                  // Best alternative
    
    void ConfigureWorkers(void);
    void ManageSingleObjectiveIterations(std::vector<std::vector<double>> parameters, int startingIndex, int numberOfParameters, std::vector<double>& objectives);
    void TerminateWorkers();

private:
    // Set the default variables for the class
    DatabaseABC* m_DbaseList = NULL;
    std::vector<std::vector<std::string>> fileListPairs;
    int m_Counter = 0;
    int m_NumCacheHits = 0;
    int m_Precision = 6;
    StringType  m_ExecCmd = NULL;
    char m_DirPrefix[DEF_STR_SZ];
    char pDirName[DEF_STR_SZ];                                          // Stem of the worker path without the rank
    FileList* m_pFileCleanupList = NULL;
    bool m_InternalModel = false;
    bool m_bCheckGlobalSens = false;
    bool m_bUseSurrogates = false;
    
    // Model preservation options
    bool m_bPreserveModelOutput = false;                                // Preserve the output from each of the model runs
    std::filesystem::path  m_PreserveOutputCmd;                         // Custom command to preserve the output from each model run
    bool m_bPreserveModelBest = false;                                  // Preserve the best solution from the analysis
    std::filesystem::path  m_PreserveBestCmd;                           // Custom command to preserve best solution from the analysis

    // Model solution options
    bool m_bWarmStart = false;
    bool m_bSolveOnPrimary = false;
    bool m_bCaching = false;
    
    bool m_bDiskless = false;
    bool m_bMultiObjProblem = false;
    double* m_CurMultiObjF = NULL;

    // Solution functions
    double** CreateInitialSample(int sampleSize);
    double** CreateSample(int sampleSize);
    void Optimize(void);
    void Calibrate(void);
    void WriteMetrics(FILE* pFile);
    void WarmStart(void);
    int  GetCurrentIteration(void);   
    
    void AddDatabase(DatabaseABC* pDbase);

    // MPI communication functions
    void ConfigureWorkerDirectory(int workerRank);
    void ConfigureWorkerSolveCommand(int workerRank);
    void ConfigureWorkerArchiveCommand(int workerRank);
    void ConfigureWorkerExtraFiles(int workerRank);
    void ConfigureWorkerFilePairs(int workerRank);
    void ConfigureWorkerObservations(int workerRank);
    void ConfigureWorkerParameterGroups(int workerRank);

    void SendWorkerContinue(int workerRank, bool workerContinue);
    void SendWorkerPreserveBest(int workerRank, bool preserveModel);
    void SendWorkerParameters(int workerRank, int alternativeIndex, std::vector<double> parameters);

    

protected:
    // Set the default groups for the class
    ObjFuncType m_ObjFuncId = OBJ_FUNC_WSSE;                            // Objective function type
    ObservationGroup* m_pObsGroup = NULL;
    ObjectiveFunction* m_pObjFunc = NULL;
    ParameterGroup* m_pParamGroup = NULL;
    DecisionModule* m_pDecision = NULL;
    ParameterCorrection* m_pParameterCorrection = NULL;
    TelescopeType m_Telescope;                                          // Telescoping bounds strategy

//protected: //can be called by DecisionModule
//    double StdExecute(double viol);
//    friend class DecisionModule;

//protected: //can be called by SuperMUSE class
//    double GatherTask(char* pDir);
//    friend class SuperMUSE;
}; /* end class Model */

/******************************************************************************
class SurrogateModel

******************************************************************************/
/*class SurrogateModel : public AlgorithmABC2 {
public:
    SurrogateModel(UnmoveableString pFileName, AlgorithmABC2* pComplex, char* pType);
    ~SurrogateModel(void) { DBG_PRINT("SurrogateModel::DTOR"); Destroy(); }
    void Destroy(void);

    //retrieve member variables
    ObservationGroup* GetObsGroupPtr(void);
    ParameterGroup* GetParamGroupPtr(void) { return NULL; }
    ObjectiveFunction* GetObjFuncPtr(void);
    double GetObjFuncVal(void) { return m_CurObjFuncVal; }
    void SetObjFuncVal(double curVal) { m_CurObjFuncVal = curVal; }
    int                 GetCounter(void);
    ObjFuncType GetObjFuncId(void) { return m_ObjFuncId; }
    UnchangeableString GetObjFuncStr(void);
    UnchangeableString GetModelStr(void) { return m_ExecCmd; }

    SurrogateParameterGroup* GetSurrogateParamGroupPtr(void);
    void Bookkeep(bool bFinal) { return; }
    int GetNumDigitsOfPrecision(void) { return 6; }
    bool CheckWarmStart(void) { return false; }

    //misc. member functions     
    double Execute(void);
    void Execute(double* pF, int nObj) { return; }
    void Write(double objFuncVal);
    void WriteMetrics(FILE* pFile);
    void SaveBest(int id) { return; }
    TelescopeType GetTelescopingStrategy(void) { return TSCOPE_NONE; }
    void PerformParameterCorrections(void) { return; }

private:
    ObjFuncType         m_ObjFuncId;
    ObservationGroup* m_pObsGroup;
    ObjectiveFunction* m_pObjFunc;
    SurrogateParameterGroup* m_pParamGroup;

    FilePair* m_FileList;
    int m_Counter;
    StringType  m_ExecCmd;
    StringType  m_pTypeStr;
    double m_CurObjFuncVal;

    void SetCmdToExecModel(IroncladString cmd);
    void AddFilePair(FilePair* pFilePair);
}; /* end class SurrogateModel */
#endif /* MODEL_H */

