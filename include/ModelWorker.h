// TODO: Doc string

#ifndef MODEL_WORKER_H
#define MODEL_WORKER_H

#include "MyHeaderInc.h"
#include "MyTypes.h"
#include <mpi.h>
#include <string>
#include <vector>
#include <filesystem>
#include <algorithm>

#include <Observation.h>
#include <FilePair.h>
#include <Utility.h>
#include <ParameterGroup.h>
#include "ParameterCorrection.h"
#include <ParameterABC.h>
#include "Exception.h"
#include "ValueExtractor.h"
#include "ObservationGroup.h"
#include "ObjectiveFunction.h"
#include "WriteUtility2.h"


enum OBJECTIVE_TYPE { single, multi };

/******************************************************************************
class ModelWorker

******************************************************************************/
class ModelWorker {
public:
	ModelWorker(void);
	void Destroy(void);
	double ExecuteSingle(void);
	void ExecuteMulti(double* pF, int nObj);
	double DisklessExecute(void);
	void PreserveModel(int rank, int trial, int counter, IroncladString ofcat);
	void   Write(double objFuncVal);

	//bool CheckCache(double* val);

	void Work(void);


private:
	IroncladString GetObjFuncCategory(double* pF, int nObj);
	void SetCmdToExecModel(IroncladString cmd);

	// Preservation variables
	bool preserveModelOutput = false;
	std::string preserveCommand;

	// Configuration variables
	std::string workerDirectory;
	std::vector<std::string> fileCleanupList;
	std::vector<std::vector<std::string>> filePairs;
	ObservationGroup *observationGroup;
	ParameterGroup *paramGroup;
	
	// Solution variables
	int solveCounter = 0;
	bool disklessModel = false;
	bool internalModel = false;
	int objectiveType = single;
	std::string solveCommand;


	// Workflow functions
	void SetupFromPrimary(void);
	void SetupWork(void);
	void CommenceWork(void);
	void TerminateWork(void);

	void ReadObservations(void);
	double StdExecute(double viol);


	// MPI communication functions
	int rank;

	std::string ReceiveString(int tag_number);
	int ReceiveInteger(int tag_number);
	double ReceiveDouble(int tag_number);
	char ReceiveChar(int tag_number);

	void SendDouble(int tag_number, double value);
	void SendInt(int tag_number, int value);

	// Setup functions 
	void ReceiveWorkerDirectory(void);
	void ReceiveWorkerSolveCommand(void);
	void ReceiveWorkerExtraFiles(void);
	void ReceiveWorkerFilePairs(void);
	void ReceiveWorkerObservations(void);
	void ReceiveWorkerParameters(void);

	int RequestParameters(void);
	bool RequestContinue(void);


};

#endif 