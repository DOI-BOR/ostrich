// TODO: Doc string

#ifndef MODEL_WORKER_H
#define MODEL_WORKER_H

#include "MyHeaderInc.h"
#include <mpi.h>
#include <string>
#include <vector>
#include <filesystem>

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
#include "WriteUtility.h"



enum MPI_TAGS {
	tag_directory, tag_textLength, tag_textFile, tag_fileLength, tag_filePairs, tag_obsLengthNum, tag_obsLengthGroup,
	tag_obsName, tag_obsValue, tag_obsWeight, tag_obsFile, tag_obsKeyword, tag_obsLine, tag_obsColumn,
	tag_obsToken, tag_obsAugmented, tag_obsGroup, tag_paramTotalNum, tag_paramTotalReal,
	tag_paramRealName, tag_paramRealInit, tag_paramRealLower, tag_paramRealUpper, tag_paramRealIn,
	tag_paramRealOst, tag_paramRealFmt, tag_paramTotalInt, tag_paramInitName, tag_paramIntInit,
	tag_paramIntLower, tag_paramIntUpper, tag_data, tag_continue
};

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



private:
	IroncladString GetObjFuncCategory(double* pF, int nObj);
	void SetCmdToExecModel(IroncladString cmd);

	// Preservation variables
	bool preserveModelOutput = false;
	std::string preserveCommand = NULL;

	// Configuration variables
	std::string workerDirectory = NULL;
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
	void Work(void);
	void SetupWork(void);
	void CommenceWork(void);
	void TerminateWork(void);
	
	void ReadObservations(void);


	// MPI communication functions
	int rank;
	void ReceiveWorkerDirectory(void);
	void ReceiveWorkerExtraFiles(void);
	void ReceiveWorkerFilePairs(void);
	void ReceiveWorkerObservations(void);
	void ReceiveWorkerParameters(void);

	int RequestParameters(void);
	bool RequestContinue(void);

	void SendDouble(int tag_number, double value);




protected: //can be called by DecisionModule
    double StdExecute(double viol);

};

#endif 