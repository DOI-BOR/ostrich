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

// Include C classes
#include <mpi.h>
#include <math.h>

// Include C++ classes
#include <string>
#include <iostream>
#include <chrono>
#include <thread>
#include <filesystem>
#include <cstring>
#include <vector>
#include <algorithm>

// Include custom classes
#include "MyHeaderInc.h"
#include "MyTypes.h"
#include "WriteUtility2.h"
#include "ModelWorker.h"
#include "ObservationGroup.h"
#include "Observation.h"
#include "ResponseVarGroup.h"
#include "SurrogateParameterGroup.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "TiedParamABC.h"
#include "GeomParamABC.h"
#include "FilePair.h"
#include "FileList.h"
#include "AccessConverter.h"
#include "NetCDFConverter.h"
#include "PumpAndTreat.h"
#include "DecisionModule.h"
#include "SuperMUSE.h"
#include "ParameterCorrection.h"
#include "GenConstrainedOpt.h"

#include "IsoParse.h"
#include "BoxCoxModel.h"
#include "Utility.h"
#include "Exception.h"
#include "SuperMuseUtility.h"


class Algorithm {
public:
    // Define constructors and destructors
    Algorithm(void);
    ~Algorithm(void) { DBG_PRINT("AlgorithmABC2::DTOR"); Destroy(); }
    void Destroy(void);

    // Solution variables
    bool m_bWarmStart = false;                                                      // Start the solution from a previously terminated run
    StringType  m_ExecCmd = NULL;                                                   // Command used to solve the model
    bool m_bSolveOnPrimary = false;                                                 // Solve on the primary worker in addition to the secondary workers

    // Objective variables
    double m_BestObjective = INFINITY;                                              // Best objective
    std::vector<double> m_BestAlternative;                                          // Best alternative

    // Defiine objective functions
    ObjectiveFunction* GetObjFuncPtr(void);
    double GetObjectiveFunctionValue(void) { return m_BestObjective; }              // Gets the current best objective function value
    void SetObjectiveFunctionValue(double curVal) { m_BestObjective = curVal; }     // Sets the current best objective function value
    ObjFuncType GetObjectiveFunctionType(void) { return m_ObjFuncId; }              // Gets the current object function type
    UnchangeableString GetObjectiveFunctionString(void);                            // Gets the string describing the current objective function
    
    void GetBestSingleObjective(std::vector<double> objectives, double& bestObjective, int& bestIndex); //  Gets the best objective function and index from a vector

    // Define public functions
    ObservationGroup* GetObsGroupPtr(void);                                         // Function to get the observation group pointer
    ParameterGroup* GetParamGroupPtr(void);                                         // Function to get the parameter group pointer
    void PerformParameterCorrections(void);                                         // Function to correct any parameters
    //void   CheckGlobalSensitivity(void);                                          // Function to check the gobal sensitivity
    int GetNumDigitsOfPrecision(void) { return m_Precision; }                       // Function to get the number of model precision digits
    
    // Expose functions to inheriting subclasses
    void ConfigureWorkers(void);                                                    // Configure the workders for the solution
    void ManageSingleObjectiveIterations(std::vector<std::vector<double>> parameters, int startingIndex, int numberOfParameters, 
                                         std::vector<double>& objectives);          // Solve using a single objective function
    void TerminateWorkers();                                                        // Terminate all workers

private:
    // Working directory information
    char pDirName[DEF_STR_SZ];                                                      // Working directory for the analysis
    char m_DirPrefix[DEF_STR_SZ];                                                   // Stem of the worker path without the rank

    // File information
    std::vector<std::vector<std::string>> fileListPairs;                            // File template/output pairs
    FileList* m_pFileCleanupList = NULL;                                            // Files to cleanup after a run

    // Caching setup information
    bool m_bCaching = false;                                                        // Flag to indicate if caching should be enabled
    int m_NumCacheHits = 0;                                                         // Counter for the number of cache hits
    std::vector<std::vector<double>> m_CacheMembers;                                // Vector to keep track of the previously calculated alternatives
    std::vector<double> m_CacheObjectives;                                          // Objectives associated with each chace member

    // Solution information
    int m_Precision = 6;                                                            // Precision that should be used unless otherwise specified
    int m_NumSolves = 0;                                                            // Number of times the model has been solved
    ModelWorker m_primaryWorker;                                                    // ModelWorker slot if using solve on primary
    
    // Sensitivity information
    bool m_bCheckGlobalSens = false;                                                // Flag to indicate if sensitivity information should be calculated

    // Surrogate model information
    bool m_bUseSurrogates = false;                                                  // Flag to indicate if surrogate information should be used
    
    // Model preservation options
    bool m_bPreserveModelOutput = false;                                            // Preserve the output from each of the model runs
    std::filesystem::path  m_PreserveOutputCmd;                                     // Custom command to preserve the output from each model run
    
    bool m_bPreserveModelBest = false;                                              // Preserve the best solution from the analysis
    std::filesystem::path  m_PreserveBestCmd;                                       // Custom command to preserve best solution from the analysis

    // Model solution options
    bool m_bDiskless = false;                                                       // Flag  to indicate if diskless solution should be done
    
    // Objective information    
    bool m_bMultiObjProblem = false;                                                // Flag to indicate if multi-objective optimization should be used
    double* m_CurMultiObjF = NULL;

    // Solution functions
    std::vector<std::vector<double>> CreateInitialSample(int sampleSize);           // Function to create an initial sample
    std::vector<std::vector<double>> CreateSample(int sampleSize);                  // Function to create a sample for subsquent iterations
    void Optimize(void);                                                            // Function to call the algorithm in optimization mode
    void Calibrate(void);                                                           // Function to call the algorithm in LS fitting mode                                              
    void WarmStart(void);                                                           // Function to warm start from a previous solution
    void AddDatabase(DatabaseABC* pDbase);                                          // Function to output to a database
    void ManagePreserveBest(double& solutionObjective, double alternativeObjective, MPI_Status mpiStatus); // Function to preserve model solves

    // MPI communication functions
    void ConfigureWorkerDirectory(int workerRank);
    void ConfigureWorkerSolveCommand(int workerRank);
    void ConfigureWorkerArchiveCommand(int workerRank, bool bMPI);
    void ConfigureWorkerExtraFiles(int workerRank, bool bMPI);
    void ConfigureWorkerFilePairs(int workerRank, bool bMPI);
    void ConfigureWorkerObservations(int workerRank, bool bMPI);
    void ConfigureWorkerParameterGroups(int workerRank, bool mMPI);

    void SendWorkerContinue(int workerRank, bool workerContinue);
    void SendWorkerPreserveBest(int workerRank, bool preserveModel);
    void SendWorkerParameters(int workerRank, int alternativeIndex, std::vector<double> parameters);

    void ReceiveWorkerPreserveBest(void);    

protected:
    // Set the default groups for the class
    ObjFuncType m_ObjFuncId = OBJ_FUNC_WSSE;                                        // Objective function type
    ObservationGroup* m_pObsGroup = NULL;                                           // Observation group
    ObjectiveFunction* m_pObjFunc = NULL;                                           // Objective function
    ParameterGroup* m_pParamGroup = NULL;                                           // Parameter group
    DecisionModule* m_pDecision = NULL;                                             // Decision module
    ParameterCorrection* m_pParameterCorrection = NULL;                             // Parameter correction
    //TelescopeType m_Telescope;                                                    // Telescoping bounds strategy


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
    double GetObjectiveFunctionValue(void) { return m_CurObjFuncVal; }
    void SetObjectiveFunctionValue(double curVal) { m_CurObjFuncVal = curVal; }
    int                 GetCounter(void);
    ObjFuncType GetObjectiveFunctionType(void) { return m_ObjFuncId; }
    UnchangeableString GetObjectiveFunctionString(void);
    UnchangeableString GetModelStr(void) { return m_ExecCmd; }

    SurrogateParameterGroup* GetSurrogateParamGroupPtr(void);
    int GetNumDigitsOfPrecision(void) { return 6; }
    bool CheckWarmStart(void) { return false; }

    //misc. member functions     
    double Execute(void);
    void Execute(double* pF, int nObj) { return; }
    void Write(double objFuncVal);
    void SaveBest(int id) { return; }
    TelescopeType GetTelescopingStrategy(void) { return TSCOPE_NONE; }
    void PerformParameterCorrections(void) { return; }

private:
    ObjFuncType         m_ObjFuncId;
    ObservationGroup* m_pObsGroup;
    ObjectiveFunction* m_pObjFunc;
    SurrogateParameterGroup* m_pParamGroup;

    FilePair* m_FileList;
    int m_NumSolves;
    StringType  m_ExecCmd;
    StringType  m_pTypeStr;
    double m_CurObjFuncVal;

    void SetCmdToExecModel(IroncladString cmd);
    void AddFilePair(FilePair* pFilePair);
}; /* end class SurrogateModel */
#endif 

