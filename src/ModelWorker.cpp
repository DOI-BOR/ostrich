// TODO: doc string

#include "ModelWorker.h"
#include "Exception.h"

#include <iostream>


/*
------------------------------------------------------------------------------------------------------------------------------
 Constructor for the model worker class
------------------------------------------------------------------------------------------------------------------------------
*/
ModelWorker::ModelWorker(void) {

    // Get the MPI information
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Get the information through MPI from the primary worker
    SetupFromPrimary();  
}

void ModelWorker::Work(void) {

    // Setup the worker directory 
    SetupWork();

    // Enter the solution loop
    CommenceWork();

    // Cleanup the worker

    // Terminate the worker
    TerminateWork();
}


void ModelWorker::Destroy(void) {
    

    IncDtorCount();
};

void ModelWorker::SetupFromPrimary() {  

    // Get the secondary worker directory tag from the primary worker
    ReceiveWorkerDirectory();

    // Get the solve command
    ReceiveWorkerSolveCommand();

    // Get the extra files from the directory
    ReceiveWorkerExtraFiles();

    // Get the file pairs for the templates and working files
    ReceiveWorkerFilePairs();

    // Get the observations for comparison
    ReceiveWorkerObservations();

    // Get the parameters for the analysis
    ReceiveWorkerParameters();
}


void ModelWorker::ReceiveWorkerDirectory(void) {
    // TODO: Doc string

    // Request the directory from the primary worker
    workerDirectory = ReceiveString(tag_directory);

    // Remove the last character. This is necessary for the numbers to append correctly.
    workerDirectory.pop_back();

    // Append the worker number onto the tag
    workerDirectory += std::to_string(rank);

}

void ModelWorker::ReceiveWorkerSolveCommand(void) {
    // TODO: Doc string

    // Request the directory from the primary worker
    solveCommand = ReceiveString(tag_solve);

    // Remove the last character.
    solveCommand.pop_back();

}


void  ModelWorker::ReceiveWorkerExtraFiles(void) {
    // TODO: doc string

    // Get the number of entries that will be coming
    int numberOfFiles = ReceiveInteger(tag_textLength);

    // Loop and request each file from the primary worker
    for (int entryFile = 0; entryFile < numberOfFiles; entryFile++) {
        // Receive and append into the holding vector
        fileCleanupList.push_back(ReceiveString(tag_textFile));
    }
}


 void ModelWorker::ReceiveWorkerFilePairs(void) {
    // TODO: doc string

    // Get the number of entries that will be coming
    int numberOfFiles = ReceiveInteger(tag_fileLength);

    // Create a vector to hold the files
    std::vector<std::string> templateFiles;
    std::vector<std::string> destinationFiles;

    // Loop and request each file from the primary worker. This will first be given as the template followed by the 
    // destination
    for (int entryFile = 0; entryFile < numberOfFiles; entryFile++) {
        for (int entryPair = 0; entryPair < 2; entryPair++) {
            // Receive the string and append into the holding vector
            if (entryPair == 0) {
                templateFiles.push_back(ReceiveString(tag_filePairs));
            }
            else {
                destinationFiles.push_back(ReceiveString(tag_filePairs));
            }
        }
    }

    // Reformat into a list pair
    for (int entryPair = 0; entryPair < templateFiles.size(); entryPair++) {
        std::vector<std::string> temp{ templateFiles[entryPair], destinationFiles[entryPair] };
        filePairs.push_back(temp);
    }

}


void ModelWorker::ReceiveWorkerObservations(void) {
    // TODO: Doc string


    // Get the number of observatiosn being transferred from the primary worker
    int numberOfObservations = ReceiveInteger(tag_obsLengthNum);
    int numberOfGroups = ReceiveInteger(tag_obsLengthGroup);

    // Create an array of observations
    Observation** observationArray = new Observation*[numberOfObservations];

    // Create a vector containing the unique files
    std::vector<std::string> uniqueFiles;

    // Loop and request each observation from the primary worker. This will loop through all obsevations, transmitting all the 
    // information for each observation before moving onto the next
    for (int entryFile = 0; entryFile < numberOfObservations; entryFile++) {

        // Receive the observation name
        std::string tempName = ReceiveString(tag_obsName);
        IroncladString obsName = &tempName[0];

        // Receive the observation value
        double obsValue = ReceiveDouble(tag_obsValue);

        // Receive the observation weight
        double obsWeight = ReceiveDouble(tag_obsWeight);

        // Receive the observation file
        std::string tempFile = ReceiveString(tag_obsFile);
        IroncladString obsFile = &tempFile[0];

        // Check if the observation file is already in use
        if (std::find(uniqueFiles.begin(), uniqueFiles.end(), tempFile) == uniqueFiles.end()) {
            // File is not in the list of files. Add to the vector.
            uniqueFiles.push_back(tempFile);
        }

        // Receive the observation keyword
        std::string tempKeyword = ReceiveString(tag_obsKeyword);
        IroncladString obsKeyword = &tempKeyword[0];

        // Recieve the observation line
        int obsLine = ReceiveInteger(tag_obsLine);

        // Receive the observation column
        int obsColumn = ReceiveInteger(tag_obsColumn);

        // Receive the observation token
        char obsToken = ReceiveChar(tag_obsToken);

        // Receive the augmented boolean
        bool obsAugmented = false;
        int obsAugmentedValue = ReceiveInteger(tag_obsAugmented);

        if (obsAugmentedValue == 1) {
            obsAugmented = true;
        }

        // Receive the observation group
        std::string tempGroup = ReceiveString(tag_obsGroup);
        IroncladString obsGroup = &tempGroup[0];

        // Create an observation object from the transferred data
        observationArray[entryFile] = new Observation(obsName, obsValue, obsWeight, obsFile, obsKeyword, obsLine, obsColumn, obsToken,
                                                      obsAugmented, obsGroup);

    }

    // Construct the value extractors
    // Create the value extractors for the observation group
    std::string temp = uniqueFiles[0];
    IroncladString temp2 = &temp[0];

    ValueExtractor *observationExtractors = new ValueExtractor(temp2, false, 1e6);
    for (int entryExtractor = 1; entryExtractor < uniqueFiles.size(); entryExtractor++) {
        // TODO: Update thse with nondefault values
        std::string tempExtractor = uniqueFiles[entryExtractor];
        observationExtractors->Insert(&tempExtractor[0]);
    }

    for (int i = 0; i < uniqueFiles.size(); i++) {
        observationExtractors->GetNext();
    }

    // Construct the group
    observationGroup = new ObservationGroup(observationArray, observationExtractors, numberOfObservations, numberOfGroups);

}


void ModelWorker::ReceiveWorkerParameters(void) {
    // TODO: Doc string

    // Get the total number of paramters
    int numberOfTotalParameters = ReceiveInteger(tag_paramTotalNum);

    // Create a counter for the number of excluded parameters
    int numberOfExcluded = 0;

    // Create a parameter group object to store everything in
    ParameterABC** m_pList = new ParameterABC *[numberOfTotalParameters];
    char** m_ParamNameList = new char* [numberOfTotalParameters];
    ParameterABC** m_pExcl = new ParameterABC *[numberOfTotalParameters];
    int positionCounter = 0;

    // Real parameters
    // Get the number of real parameters
    int numberOfRealParameters = ReceiveInteger(tag_paramTotalReal);

    // Request the values from the primary worker
    for (int entryReal = 0; entryReal < numberOfRealParameters; entryReal++) {

        // Receive the values from the primary worker
        std::string tempName = ReceiveString(tag_paramRealName);
        IroncladString paramName = &tempName[0];

        double paramInitial = ReceiveDouble(tag_paramRealInit);
        double paramLowerBound = ReceiveDouble(tag_paramRealLower);
        double paramUpperBound = ReceiveDouble(tag_paramRealUpper);
        
        std::string tempIn = ReceiveString(tag_paramRealIn);
        IroncladString paramTxIn = &tempIn[0];

        std::string tempOst = ReceiveString(tag_paramRealOst);
        IroncladString paramTxOst = &tempOst[0];

        std::string tempOut = ReceiveString(tag_paramRealOut);
        IroncladString paramTxOut = &tempOut[0];

        std::string tempFmt = ReceiveString(tag_paramRealFmt);
        IroncladString paramFmt = &tempFmt[0];

        m_pList[positionCounter] = new RealParam(paramName, paramInitial, paramLowerBound, paramUpperBound, paramTxIn, paramTxOst, paramTxOut, paramFmt);
        
        m_ParamNameList[positionCounter] = new char[strlen(paramName) + 1];
        strcpy(m_ParamNameList[positionCounter], paramName);
        
        m_pExcl[positionCounter] = NULL; // No values are assumed in the exclude list

        // Increment the position counter
        positionCounter++;
    }

    // Integer parameters
    // Get the number of integer parameters
    int numberOfIntegerParameters = ReceiveInteger(tag_paramTotalInt);

    // Request the values from the primary worker
    for (int entryInt = 0; entryInt < numberOfIntegerParameters; entryInt++) {

        // Receive the values tag_paramInitName the primary worker
        std::string tempIn = ReceiveString(tag_paramInitName);
        IroncladString paramName = &tempIn[0];

        int paramInitial = ReceiveInteger(tag_paramIntInit);
        int paramLowerBound = ReceiveInteger(tag_paramIntLower);
        int paramUpperBound = ReceiveInteger(tag_paramIntUpper);

        m_pList[positionCounter] = new IntParam(paramName, paramInitial, paramLowerBound, paramUpperBound);

        m_ParamNameList[positionCounter] = new char[strlen(paramName) + 1];
        strcpy(m_ParamNameList[positionCounter], paramName);

        m_pExcl[positionCounter] = NULL; // No values are assumed in the exclude list

        // Increment the position counter
        positionCounter++;
    }

    // Special parameters
    SpecialParam **m_pSpecial = new SpecialParam * [0];

    // Tied parameters
    TiedParamABC **m_pTied = new TiedParamABC * [0];

    // Geom parameters
    GeomParamABC **m_pGeom = new GeomParamABC * [0];

    // Set the arrays into a group object
    // Initialize a parameter group
    paramGroup = new ParameterGroup(false);
    
    // Set the parameter values into the group
    paramGroup->SetGroupValues(m_pList, m_pExcl, m_pTied, m_pGeom, m_pSpecial, m_ParamNameList, numberOfTotalParameters,  0,  0, 0,  numberOfExcluded);
}


std::string  ModelWorker::ReceiveString(int tag_number) {
    // TODO: Doc string

    // Initialize the variables required for the probe
    MPI_Status status;
    int characterLength;

    // Conduct the probe on the request
    MPI_Probe(0, tag_number, MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, MPI_CHAR, &characterLength);

    // Allocate the character array dynamically
    char* directoryFromPrimary = new char[characterLength];

    // Request the directory from the primary worker
    MPI_Recv(directoryFromPrimary, characterLength, MPI_CHAR, 0, tag_number, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Convert to a string and return to the calling function
    std::string convertedArray = std::string(directoryFromPrimary, characterLength);
    delete directoryFromPrimary;

    // Return the array to the calling function
    return convertedArray;
}


int  ModelWorker::ReceiveInteger(int tag_number) {
    // TODO: Doc string

    // Initialize the variable to fill
    int value;

    // Request the value
    MPI_Recv(&value, 1, MPI_INT, 0, tag_number, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    return value;
}


double  ModelWorker::ReceiveDouble(int tag_number) {
    // TODO: Doc string

    // Initialize the variable to fill
    double value;

    // Request the value
    MPI_Recv(&value, 1, MPI_DOUBLE, 0, tag_number, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    return value;
}


char  ModelWorker::ReceiveChar(int tag_number) {
    // TODO: Doc string

    // Initialize the variable to fill
    char value;

    // Request the value
    MPI_Recv(&value, 1, MPI_CHAR, 0, tag_number, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    return value;
}


void  ModelWorker::SendDouble(int tag_number, double value) {
    // TODO: Doc string

    // Send the value
    MPI_Send(&value, 1, MPI_DOUBLE, 0, tag_number, MPI_COMM_WORLD);

}

void  ModelWorker::SendInt(int tag_number, int value) {
    // TODO: Doc string

    // Send the value
    MPI_Send(&value, 1, MPI_INT, 0, tag_number, MPI_COMM_WORLD);

}


void ModelWorker::SetupWork(void) {

    // Create the working directory
    if (!std::filesystem::exists(workerDirectory)) {
        // Directory does not already exist. Create it.
        std::filesystem::create_directories(workerDirectory);
    }

    // Copy all the files to it from the source directory
    for (int entryFile = 0; entryFile < fileCleanupList.size(); entryFile++) {

        // GEt the current working directory
        std::filesystem::path workerDirectoryPath = workerDirectory;
        
        // Remove the last character which causes issues.
        // TODO: remove this during the read of the control file
        std::string temp = fileCleanupList[entryFile];
        temp.pop_back();

        // Create the filespath
        std::filesystem::path filePath = workerDirectoryPath /= temp;
        std::filesystem::path testPath = filePath.remove_filename();

        // If the directory structure does not exist, create it
        if (!std::filesystem::exists(testPath)) {
            std::filesystem::create_directories(testPath);
        }
        
        // Create the source path
        std::filesystem::path sourcePath = fileCleanupList[entryFile];

        // Copy from the source path to the worker path
        std::error_code code;
        std::filesystem::copy(sourcePath, filePath, code);

        // Handle any errors in a user friendly mannter.
        if (code.value() == 2){
            std::cout << "Extra file " << fileCleanupList[entryFile] << " is not found." << std::endl;
            LogError(ERR_FILE_IO, "Using default algorithm setup.");
        }
    }

    // Perform any checks on the model setup  

}


void  ModelWorker::CommenceWork() {
    // TODO: doc string

    // Define control variable to break the solution loop
    bool continueWork = true;

    // Enter the secondar working processing loop
    while (continueWork) {

        // Determine if the worker should continue. The primary must send work to secondary if it is told to continue.
        continueWork = RequestContinue();

        if (continueWork) {
            // Request parameter values from the primary worker
            int alternativeIndex = RequestParameters();

            // Solve the model
            if (objectiveType == single) {
                // Call the solve function
                double objective = ExecuteSingle();

                // Stack the alternative and the objective into a single array to transimit 
                double returnArray[2]{ (double)alternativeIndex, objective };
                
                // Send to the primary worker
                MPI_Send(returnArray, 2, MPI_DOUBLE, 0, tag_data, MPI_COMM_WORLD);

            }
            else {
                // Calculate the objective function
                // TODO: Add multiobject evaluation
            }

            // Perform any inter-call cleanup
        }
        
    }
}

void  ModelWorker::TerminateWork(void) {
    // Todo: update when additional operations are required at work termination

}

int ModelWorker::RequestParameters(void) {

    // Create holder for the parameters. This is one larger than the number of parameters to allow the alternative index to be transferred
    // with the alternatives
    int numberOfParameters = paramGroup->GetNumParams();
    double *parametersTemp = new double[numberOfParameters + 1];

    // Receive the parameters
    MPI_Recv(parametersTemp, numberOfParameters+1, MPI_DOUBLE, 0, tag_data, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Separate out the alternative index from the parameters
    int alternativeIndex = (int)parametersTemp[0];

    double *parameters = new double[numberOfParameters];
    for (int entryParameter = 0; entryParameter < numberOfParameters; entryParameter++) {
        parameters[entryParameter] = parametersTemp[entryParameter + 1];
    }

    // Set the parameters into the the parameter group
    paramGroup->WriteParams(parameters);

    // Cleanup the memory space
    delete parametersTemp;

    // Return the alternative index to the calling function
    return alternativeIndex;
}

bool ModelWorker::RequestContinue(void) {
    

    // Send that it's reddy to receive work
    SendInt(tag_continue, 1);

    // Request the continue flag back
    int continueValue = ReceiveInteger(tag_continue);

    // Convert the int to a bool
    if (continueValue == 0) {
        return false;
    }
    else {
        return true;
    }

}


/******************************************************************************
 PreserveModel()
******************************************************************************/
void ModelWorker::PreserveModel(int rank, int trial, int counter, IroncladString ofcat)
{
    char tmp[DEF_STR_SZ];

    if (preserveModelOutput == false) {
        return;
    }

    /* use built in preservation option --- tries to save everything */
    if (preserveCommand.empty()) {
        // Get the current path to the worker
        std::filesystem::path directory = std::filesystem::current_path();

        // Move up one directory to project directory
        directory = directory.parent_path();

        // Construct the archive folder
        directory /= std::string("archive");
        if (!std::filesystem::exists(directory)) {
            std::filesystem::create_directory(directory);
        };

        // Construct the worker folder
        std::filesystem::path workerDirectoryPath = workerDirectory;
        directory /= workerDirectoryPath.parent_path().filename();
        directory += std::to_string(rank);

        if (!std::filesystem::exists(directory)) {
            std::filesystem::create_directory(directory);
        }

        // Construct the run folder 
        directory /= "run_" + std::to_string(counter - 1); // This adjusts by one to align with the model output files
        std::filesystem::create_directory(directory);

        // Move the files from the worker to the archive
        std::filesystem::copy(std::filesystem::current_path(), directory, std::filesystem::copy_options::recursive);
    }
    else{
        /* run the user-supplied preservation command */
        sprintf(tmp, "%s %d %d %d %s", preserveCommand, rank, trial, counter, ofcat);

        #ifdef _WIN32 //windows version
            strcat(tmp, " > OstPreserveModelOut.txt");
        #else //Linux (bash, dash, and csh)
            // '>&' redircts both output and error
            strcat(tmp, " > OstPreserveModelOut.txt 2>&1");
        #endif

        system(tmp);    
    }
}/* end PreserveModel() */


/*****************************************************************************
Execute()
   Executes the model (or surrogate) and returns the objective function value.
******************************************************************************/
double ModelWorker::ExecuteSingle(void) {
    
    // Define the objective value
    double objective;
    objective = StdExecute(0.00);

    // TODO: reenable other execute methods
    /*if (disklessModel == true && internalModel == true) {
        objective = DisklessExecute();
    }
    else {
        objective = StdExecute(0.00);
    }
    else {
        objective = Execute();
    }*/

    //if desired update log of residuals
    WriteIterationResiduals();

    // Return to the calling function
    return objective;
}

/*****************************************************************************
Execute()
   Executes the model returns a vector of objective function values.
******************************************************************************/
/*void ModelWorker::ExecuteMulti(double* pF, int nObj) {
    IroncladString dirName = GetExeDirName();
    FilePair* pCur;
    FilePipe* pPipe;
    int rank;

    //initialize costs
    for (int i = 0; i < nObj; i++)
    {
        pF[i] = NEARLY_HUGE;
    }

    //exit early if the user has requested program termination
    if (IsQuit() == true)
    {
        return;
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //inc. number of times model has been executed
    m_Counter++;

    //make substitution of parameters into model input file
    pCur = m_FileList;
    while (pCur != NULL)
    {
        pPipe = pCur->GetPipe();
        m_pParamGroup->SubIntoFile(pPipe);
        pCur = pCur->GetNext();
    }

    //cd to model subdirectory, if needed
    if (dirName[0] != '.') { MY_CHDIR(dirName); }

    //make substitution of parameters into model input databases
    if (m_DbaseList != NULL)
    {
        m_pParamGroup->SubIntoDbase(m_DbaseList);
    } 

    //invoke system command to execute the model   
    system(m_ExecCmd);

    //extract computed reponses from model output database(s)
    if (m_DbaseList != NULL)
    {
        DatabaseABC* pCur;

        //first clean up ASCII files
        for (pCur = m_DbaseList; pCur != NULL; pCur = pCur->GetNext())
        {
            pCur->DeleteASCIIFile();
        }

        //convert the responses to ASCII file
        for (pCur = m_DbaseList; pCur != NULL; pCur = pCur->GetNext())
        {
            pCur->ReadResponse();
        }
    }

    //extract computed observations from model output file(s)
    if (m_pObsGroup != NULL) { m_pObsGroup->ExtractVals(); }

    //compute obj. func.
    m_pObjFunc->CalcMultiObjFunc(pF, nObj);

    //categorize the obj. func.
    if (dirName[0] != '.') { MY_CHDIR(".."); }
    IroncladString pCatStr = GetObjFuncCategory(pF, nObj);
    if (dirName[0] != '.') { MY_CHDIR(dirName); }

    //preserve model output, if desired
    PreserveModel(rank, GetTrialNumber(), m_Counter, pCatStr);

    //cd out of model subdirectory, if needed
    if (dirName[0] != '.') { MY_CHDIR(".."); }

    //ouput results (first objective, others should be tagged as augmented responses)
    double val = pF[0];
    m_CurObjFuncVal = val;

    //store copy of latest result, needed for printing purposes
    if (m_CurMultiObjF == NULL) m_CurMultiObjF = new double[nObj];
    for (int i = 0; i < nObj; i++)
    {
        m_CurMultiObjF[i] = pF[i];
    }

    Write(val);
} */

/*****************************************************************************
StdExecute()
   Executes the standard (complex) model and returns the objective function
   value.
******************************************************************************/
double ModelWorker::StdExecute(double viol) {
    //IroncladString dirName = &workerDirectory[0];

    bool isGoodTopo;

    //inc. number of times model has been executed
    solveCounter++;

    /*
    adjust geometries to conform to topology rules, if the topology cannot be 'fixed' then false will be returned....
    */
    isGoodTopo = paramGroup->FixGeometry();
    if (isGoodTopo == false) {
        LogError(ERR_MODL_EXE, "Could not correct model topology");
    }

    // Make substitution of parameters into model input files
    for (int entryPair = 0; entryPair < filePairs.size(); entryPair++) {
        // Setup the input file and path
        std::string inTemp = filePairs[entryPair][0];
        IroncladString inFile = &inTemp[0];

        // Setup the output file and path
        std::string outFileString = filePairs[entryPair][1];
        outFileString.pop_back();
        std::filesystem::path outTemp = workerDirectory;
        outTemp = outTemp /= outFileString;

        outFileString = outTemp.string();
        IroncladString outFile = &outFileString[0];

        // Define a pair and pipe to allow substitution
        FilePair *filePair = new FilePair(inFile, outFile);
        FilePipe *filePile = filePair->GetPipe();

        // Substitute the parameters into the file
        paramGroup->SubIntoFile(filePile);

        // Delete the file pair 
        delete filePair;

    } 

    // Move to model subdirectory, if needed
    std::filesystem::path parentPath = std::filesystem::current_path();
    std::filesystem::current_path(workerDirectory);


    //make substitution of parameters into model input databases
    /*if (m_DbaseList != NULL)
    {
        m_pParamGroup->SubIntoDbase(m_DbaseList);
    } */

    /* -----------------------------------------------------------
    If caching is enabled, attempt to read model evaluation
    from OstModel0.txt file.
    ----------------------------------------------------------- */
    /*if (m_bCaching == true)
    {
        bool bCached = CheckCache(&val);
        if (bCached == true) //found previous model result
        {
            m_NumCacheHits++;
            Write(val);
            m_CurObjFuncVal = val;
            return (val);
        }
    }*/

    // Solve the model 
    /*if (internalModel == true) {
        if (strcmp(m_ExecCmd, "Isotherm()") == 0) { Isotherm(m_bDiskless); }
        else if (strcmp(m_ExecCmd, "Orear()") == 0) { Orear(); }
        else if (strcmp(m_ExecCmd, "McCammon()") == 0) { McCammon(m_bDiskless); }
        else if (strcmp(m_ExecCmd, "Kinniburgh()") == 0) { Kinniburgh(m_bDiskless); }
        else if (strcmp(m_ExecCmd, "AdvancedKinniburgh()") == 0) { AdvancedKinniburgh(); }
        else if (strcmp(m_ExecCmd, "BoxCox()") == 0) { BoxCoxModel(); }
        else { LogError(ERR_BAD_ARGS, "Unknown internal model"); ExitProgram(1); }
    }
    else
    {
        //invoke system command to execute the model   
        system(m_ExecCmd);
    }*/
    
    // Solve the model
    system(&solveCommand[0]);

    //extract computed reponses from model output database(s)
    /*if (m_DbaseList != NULL) {
        DatabaseABC* pCur;

        //first clean up ASCII files
        for (pCur = m_DbaseList; pCur != NULL; pCur = pCur->GetNext())
        {
            pCur->DeleteASCIIFile();
        }

        //convert the responses to ASCII file
        for (pCur = m_DbaseList; pCur != NULL; pCur = pCur->GetNext())
        {
            pCur->ReadResponse();
        }
    }*/

    //extract computed observations from model output file(s)
    observationGroup->ExtractVals(); 

    //compute obj. func.
    // TODO: add other objective types
    // TODO: add box cox variables
    WSSE objFunc = WSSE(observationGroup, false, 0);
    double objValue = objFunc.CalcObjFunc();

    //add in penalty for violation of parameter bounds
    objValue += viol * MyMax(1.00, objValue);

    //categorize the obj. func.
    //IroncladString pCatStr = GetObjFuncCategory(&objValue, 1);

    //preserve model output, if desired
    if (preserveModelOutput) {
        // TODOL Reenable preservation with a user command 
        //PreserveModel(rank, GetTrialNumber(), solveCounter, NULL);
    }

    // Output results
    Write(objValue);

    // Change back to the parent path for the next iterator
    std::filesystem::current_path(parentPath);

    // Return to the calling function
    return objValue;
} 



/*****************************************************************************
DisklessExecute()
   Executes an inernal model without using I/O.
******************************************************************************/
/*double ModelWorker::DisklessExecute(void) {
    //inc. number of times model has been executed
    m_Counter++;

    if (strcmp(m_ExecCmd, "Isotherm()") == 0) { DisklessIsotherm(m_pParamGroup, m_pObsGroup); }
    //else if(strcmp(m_ExecCmd, "Orear()") == 0){ DisklessOrear();}
    else if (strcmp(m_ExecCmd, "McCammon()") == 0) { DisklessMcCammon(m_pParamGroup, m_pObsGroup); }
    else if (strcmp(m_ExecCmd, "Kinniburgh()") == 0) { DisklessKinniburgh(m_pParamGroup, m_pObsGroup); }
    //else if(strcmp(m_ExecCmd, "AdvancedKinniburgh()") == 0){ DisklessAdvancedKinniburgh();}
    //else if(strcmp(m_ExecCmd, "BoxCox()") == 0){ DisklessBoxCoxModel();}
    else { return 0.00; }

    //compute obj. func.
    double val = m_pObjFunc->CalcObjFunc();

    m_CurObjFuncVal = val;

    return (val);
} */

/*****************************************************************************
GetObjFuncCategory()
   Determine the obj. func. category argument to pass to the PreserveModel
   script.
******************************************************************************/
/*IroncladString ModelWorker::GetObjFuncCategory(double* pF, int nObj)
{
    switch (GetProgramType())
    {
        
    case(GA_PROGRAM):
    case(BGA_PROGRAM):
    case(SA_PROGRAM):
    case(CSA_PROGRAM):
    case(VSA_PROGRAM):
    case(PSO_PROGRAM):
    case(PSO_LEV_PROGRAM):
    case(LEV_PROGRAM):
    case(POWL_PROGRAM):
    case(BIS_PROGRAM):
    case(STEEP_PROGRAM):
    case(FLRV_PROGRAM):
    case(DDS_PROGRAM):
    case(GMLMS_PROGRAM):
    case(SCEUA_PROGRAM):
    case(DDDS_PROGRAM):
    case(SMP_PROGRAM):
    case(PDDS_PROGRAM):
    case(APPSO_PROGRAM):
    case(BEERS_PROGRAM):
    {
        if ((m_Counter <= 1) || (pF[0] < GetBestObjFunc(m_pParamGroup->GetNumParams())))
        {
            return(ObjFuncBest);
        }
        return ObjFuncOther;
    }

    // rejection samplers
    case(RJSMP_PROGRAM):
    case(METRO_PROGRAM):
    {
        return ObjFuncOther;
    }

    // behaviorial samplers
    case(GLUE_PROGRAM):
    case(DDSAU_PROGRAM):
    {
        if (pF[0] < GetObjFuncThreshold())
        {
            return(ObjFuncBehavioral);
        }
        return ObjFuncNonBehavioral;
    }

    // multi-objective optimizers
    case(SMOOTH_PROGRAM):
    case(PADDS_PROGRAM):
    case(PARA_PADDS_PROGRAM):
    {
        if ((m_Counter <= 1) || IsNonDominated(pF, nObj))
        {
            return ObjFuncNonDominated;
        }
        return ObjFuncDominated;
    }

    // miscellaneous
    case(SET_INFILE):
    case(STATS_PROGRAM):
    case(UTIL_PROGRAM):
    case(GRID_PROGRAM):
    case(EVAL_PROGRAM):
    case(JACOBIAN_PROGRAM):
    case(HESSIAN_PROGRAM):
    case(GRADIENT_PROGRAM):
    case(QUIT_PROGRAM):
    default:
    {
        return ObjFuncOther;
    }
    }
    return ObjFuncOther;
}*/

/******************************************************************************
SetCmdToExecModel()
   Sets the syntax which is used to execute the model
*******************************************************************************/
/*void ModelWorker::SetCmdToExecModel(IroncladString cmd)
{
    int len;
    len = (int)strlen(cmd) + 1;
    NEW_PRINT("char", len);
    m_ExecCmd = new char[len];
    MEM_CHECK(m_ExecCmd);

    strcpy(m_ExecCmd, cmd);
} */

/******************************************************************************
Write()
   Store parameter and objective function value to model output file.
******************************************************************************/
void ModelWorker::Write(double objFuncVal)
{
    //ResponseVarGroup* pRespVarGroup;
    FILE* pFile;
    char name[DEF_STR_SZ];

    /*pRespVarGroup = NULL;
    if (m_pObjFunc != NULL)
    {
        pRespVarGroup = (ResponseVarGroup*)(m_pObjFunc->GetResponseVarGroup());
    }*/

    // Create the name of the log file for the secondary worker
    sprintf(name, "OstModel%d.txt", rank);

    // Log the solve
    pFile = fopen(name, "a+");
    if (pFile == NULL) {
        LogError(ERR_FILE_IO, "Write(): Couldn't open OstModel.txt file");
        ExitProgram(1);
    }
    
    // Setup the header on the first iteration
    if (solveCounter == 1) {
        if (objectiveType == single) {
            fprintf(pFile, "Run   obj.function   ");
        }
        else {
            fprintf(pFile, "Run   ");
        }

        paramGroup->Write(pFile, WRITE_BNR);
        fprintf(pFile, "\n");
        fclose(pFile);
    }

    
    pFile = fopen(name, "a+");
    fprintf(pFile, "%-4d  ", solveCounter);

    if (objectiveType == single) {
        WritePreciseNumber2(pFile, objFuncVal);
        fprintf(pFile, "  ");
    }

    //if (observationGroup != NULL) observationGroup->Write(pFile, WRITE_SCI, m_CurMultiObjF);
    //if (pRespVarGroup != NULL) pRespVarGroup->Write(pFile, WRITE_SCI);
    paramGroup->Write(pFile, WRITE_SCI);
    fprintf(pFile, "\n");
    fclose(pFile);
} /* end Write() */



/******************************************************************************
CheckCache()
   Attempt to read previous result stored in OstModel0.txt file. If successful
   the corresponding obj. function value will be stored in val argument. Other-
   wise the function returns false and val is unchanged.
******************************************************************************/
/*bool ModelWorker::CheckCache(double* val)
{
    int i, j, max_line_size;
    char* line;
    char valStr[DEF_STR_SZ];
    char* pTmp;
    double objFunc, paramVal;
    bool bFound;
    FILE* pIn;

    ParameterABC* pParam;
    ParameterGroup* pGroup = GetParamGroupPtr();
    if (pGroup == NULL) return false;
    int np = pGroup->GetNumParams();

    max_line_size = GetMaxLineSizeInFile((char*)"OstModel0.txt");
    line = new char[max_line_size + 1];
    line[0] = NULLSTR;
    pIn = fopen("OstModel0.txt", "r");
    if (pIn == NULL)
    {
        delete[] line;
        return false;
    }
    while (!feof(pIn))
    {
        fgets(line, max_line_size, pIn);
        pTmp = line;
        MyTrim(pTmp);
        j = ExtractColString(pTmp, valStr, ' ');
        pTmp += j;
        MyTrim(pTmp);
        j = ExtractColString(pTmp, valStr, ' ');
        pTmp += j;
        MyTrim(pTmp);
        objFunc = atof(valStr);

        bFound = true;
        for (i = 0; i < np; i++)
        {
            pParam = pGroup->GetParamPtr(i);
            j = ExtractColString(pTmp, valStr, ' ');
            paramVal = pParam->ConvertInVal(atof(valStr));
            if (fabs(paramVal - pParam->GetEstVal()) > 1E-10)
            {
                bFound = false;
                break;
            }
            pTmp += j;
            MyTrim(pTmp);
        }
        if (bFound == true)
        {
            fclose(pIn);
            *val = objFunc;
            delete[] line;
            return true;
        }
    }
    fclose(pIn);
    delete[] line;
    return false;
}*/