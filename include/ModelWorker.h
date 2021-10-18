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
#include <iostream>
#include <chrono>
#include <thread>
#include <cctype>
#include <locale>

#include <Observation.h>
#include <FilePair.h>
#include <Utility.h>
#include "ParameterGroup.h"
#include "TiedParamABC.h"
#include <ParameterGroupWorker.h>
#include "ParameterCorrection.h"
#include <ParameterWorker.h>
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
	ModelWorker();
	ModelWorker(bool bMPI);
	void Destroy(void);
	double ExecuteSingle(void);
	//void ExecuteMulti(double* pF, int nObj);		// Todo: reimplement this
	//double DisklessExecute(void);					// Todo: reimplement this
	void PreserveModel(bool preserveBest);
	void Write(double objFuncVal);

	// Solution functions
	void SetupWorker(void);
	void WorkMPI(void);
	double StdExecute(double viol);

	// Primary worker setup functions
	void SetWorkerDirectory(std::string workerDirectory);
	void SetWorkerSolveCommand(std::string solveCommand);
	void SetWorkerPreserveArchiveCommand(int preserveStatus, std::string preserveCommand);
	void SetWorkerPreserveBestCommand(int preserveBest, std::string bestCommand);
	void SetWorkerExtraFiles(std::vector<std::string> extraFiles);
	void SetWorkerFilePairs(std::vector<std::vector<std::string>> filePairs);
	void SetWorkerObservations(ObservationGroup* m_pObsGroup);
	void SetWorkerParameters(ParameterGroup* m_pParamGroup);
	int SetStandardParameters(std::vector<double> inputParameters);


private:
	IroncladString GetObjFuncCategory(double* pF, int nObj);

	// Preservation variables
	bool preserveModelBest = false;
	std::string preserveBestCommand;
	bool preserveModelOutput = false;
	std::string preserveOutputCommand;

	// Configuration variables
	std::string m_workerDirectory;
	std::vector<std::string> m_fileCleanupList;
	std::vector<std::vector<std::string>> m_filePairs;
	ObservationGroup *observationGroup;
	ParameterGroupWorker *paramGroup;
	
	// Solution variables
	int solveCounter = 0;
	bool disklessModel = false;
	bool internalModel = false;
	int objectiveType = single;
	std::string solveCommand;

	// Solution functions
	void ReadObservations(void);

	// MPI workflow functions
	void SetupMPI(void);
	void CommenceMPIWork(void);
	void TerminateMPIWork(void);

	// MPI communication functions
	int rank;

	std::string ReceiveString(int tag_number);
	int ReceiveInteger(int tag_number);
	double ReceiveDouble(int tag_number);
	char ReceiveChar(int tag_number);

	void SendDouble(int tag_number, double value);
	void SendInteger(int tag_number, int value);

	// MPI setup functions 
	void ReceiveWorkerDirectory(void);
	void ReceiveWorkerSolveCommand(void);
	void ReciveWorkerArchiveCommand(void);
	void ReceiveWorkerExtraFiles(void);
	void ReceiveWorkerFilePairs(void);
	void ReceiveWorkerObservations(void);
	void ReceiveWorkerParameters(void);

	int RequestParameters(void);
	bool RequestContinue(void);

};

#endif 