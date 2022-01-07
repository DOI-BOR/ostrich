// TODO: doc string

#include "ModelWorker.h"
#include "Exception.h"

/**************************************************************************************************************************************************************
 ModelWorker()
 
 Default Constructor for the model worker class
**************************************************************************************************************************************************************/
ModelWorker::ModelWorker() {

}

/**************************************************************************************************************************************************************
ModelWorker()
 
Constructor for the model worker class
**************************************************************************************************************************************************************/
ModelWorker::ModelWorker(bool bMPI) {

    // Get the MPI information
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Get the information through MPI from the primary worker
    if (bMPI) {
        SetupMPI();
    }
    
}

/**************************************************************************************************************************************************************
Destroy()

Destructor for the model worker class
**************************************************************************************************************************************************************/
void ModelWorker::Destroy(void) {

    IncDtorCount();
};

/**************************************************************************************************************************************************************
WorkMPI

Runs the sequence of operations for a secondary worker called in MPI mode
**************************************************************************************************************************************************************/
void ModelWorker::WorkMPI(void) {

    // Setup the worker directory 
    SetupWorker();

    // Enter the solution loop
    CommenceMPIWork();

    // Cleanup the worker

    // Terminate the worker
    TerminateMPIWork();
}

/**************************************************************************************************************************************************************
SetupMPI

Receives configuration information from the primary worker through the MPI interface in the correct order.
**************************************************************************************************************************************************************/
void ModelWorker::SetupMPI() {  

    // Get the secondary worker directory tag from the primary worker
    ReceiveWorkerDirectory();

    // Get the solve command
    ReceiveWorkerSolveCommand();

    // Get the archive commands
    ReciveWorkerArchiveCommand();

    // Get the extra files from the directory
    ReceiveWorkerExtraFiles();

    // Get the file pairs for the templates and working files
    ReceiveWorkerFilePairs();

    // Get the observations for comparison
    ReceiveWorkerObservations();

    // Get the parameters for the analysis
    ReceiveWorkerParameters();
}

/**************************************************************************************************************************************************************
ReceiveWorkerDirectory

Receives the worker directory from the primary worker using the MPI inerface
**************************************************************************************************************************************************************/
void ModelWorker::ReceiveWorkerDirectory(void) {

    // Request the directory from the primary worker
    m_workerDirectory = ReceiveString(tag_directory);

    // Append the worker number onto the tag
    m_workerDirectory += std::to_string(rank);

}

/**************************************************************************************************************************************************************
ReceiveWorkerSolveCommand

Receives the worker solve command from the primary worker using the MPI inerface
**************************************************************************************************************************************************************/
void ModelWorker::ReceiveWorkerSolveCommand(void) {

    // Request the directory from the primary worker
    solveCommand = ReceiveString(tag_solve);

}

/**************************************************************************************************************************************************************
ReciveWorkerArchiveCommand

Receives the worker archive command from the primary worker using the MPI inerface
**************************************************************************************************************************************************************/
void ModelWorker::ReciveWorkerArchiveCommand(void) {
    
    // Request the preservation status from the primary worker
    int preserveStatus = ReceiveInteger(tag_archive);

    if (preserveStatus == 0) {
        // Model will not be archived. Set the status to false
        preserveModelOutput = false;

    } else if (preserveStatus == 1) {
        // Model will be archived using standard method
        // Toggle the preserve status to true
        preserveModelOutput = true;

    } else if (preserveStatus == 2) {
        // Model will be archived using a custom command.
        // Toggle the preserve status to true
        preserveModelOutput = true;

        // Request the preservation command from the worker
        preserveOutputCommand = ReceiveString(tag_archive);

    }

    // Request if the best model will be preserved
    preserveStatus = ReceiveInteger(tag_archive);

    if (preserveStatus == 0) {
        // Model will not be archived. Set the status to false
        preserveModelBest = false;

    }
    else if (preserveStatus == 1) {
        // Model will be archived using standard method
        // Toggle the preserve status to true
        preserveModelBest = true;

    }
    else if (preserveStatus == 2) {
        // Model will be archived using a custom command.
        // Toggle the preserve status to true
        preserveModelBest = true;

        // Request the preservation command from the worker
        preserveBestCommand = ReceiveString(tag_archive);

    }
}

/**************************************************************************************************************************************************************
ReceiveWorkerExtraFiles

Receives the worker extra files from the primary worker using the MPI inerface
**************************************************************************************************************************************************************/
void  ModelWorker::ReceiveWorkerExtraFiles(void) {

    // Get the number of entries that will be coming
    int numberOfFiles = ReceiveInteger(tag_textLength);

    // Loop and request each file from the primary worker
    for (int entryFile = 0; entryFile < numberOfFiles; entryFile++) {
        // Receive and append into the holding vector
        m_fileCleanupList.push_back(ReceiveString(tag_textFile));
    }
}

/**************************************************************************************************************************************************************
ReceiveWorkerFilePairs

Receives the worker file pairs from the primary worker using the MPI inerface
**************************************************************************************************************************************************************/
 void ModelWorker::ReceiveWorkerFilePairs(void) {

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
        m_filePairs.push_back(temp);
    }
}

 /**************************************************************************************************************************************************************
 ReceiveWorkerObservations

 Receives the worker observations from the primary worker using the MPI inerface
 **************************************************************************************************************************************************************/
void ModelWorker::ReceiveWorkerObservations(void) {

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

/**************************************************************************************************************************************************************
ReceiveWorkerParameters

Receives the worker parameters from the primary worker using the MPI inerface
**************************************************************************************************************************************************************/
void ModelWorker::ReceiveWorkerParameters(void) {

    // Get the total number of paramters
    int numberOfTotalParameters = ReceiveInteger(tag_paramTotalNum);

    // Create a counter for the number of excluded parameters
    int numberOfExcluded = 0;

    // Create a parameter group object to store everything in
    ParameterWorkerABC** m_pList = new ParameterWorkerABC *[numberOfTotalParameters];
    char** m_ParamNameList = new char* [numberOfTotalParameters];
    ParameterWorkerABC** m_pExcl = new ParameterWorkerABC *[numberOfTotalParameters];
    int positionCounter = 0;

    // Real parameters
    // Get the number of real parameters
    int numberOfRealParameters = ReceiveInteger(tag_paramTotalReal);

    // Request the values from the primary worker
    for (int entryReal = 0; entryReal < numberOfRealParameters; entryReal++) {

        // Receive the name from the primary worker
        std::string tempName = ReceiveString(tag_paramRealName);
        IroncladString paramName = &tempName[0];

        // Recive the value
        double paramInitial = ReceiveDouble(tag_paramRealInit);

        // Receive the format
        std::string tempFmt = ReceiveString(tag_paramRealFmt);
        IroncladString paramFmt = &tempFmt[0];

        // Create the parameter
        m_pList[positionCounter] = new RealParamWorker(paramName, paramInitial, paramFmt);
        
        // Construct the name list
        m_ParamNameList[positionCounter] = new char[strlen(paramName) + 1];
        strcpy(m_ParamNameList[positionCounter], paramName);
        
        // Update the exclusion list
        m_pExcl[positionCounter] = NULL; // No values are assumed in the exclude list

        // Increment the position counter
        positionCounter++;
    }

    // Integer parameters
    // Get the number of integer parameters
    int numberOfIntegerParameters = ReceiveInteger(tag_paramTotalInt);

    // Request the values from the primary worker
    for (int entryInt = 0; entryInt < numberOfIntegerParameters; entryInt++) {

        // Receive the name from the primary worker
        std::string tempIn = ReceiveString(tag_paramInitName);
        IroncladString paramName = &tempIn[0];

        // Receive the value
        int paramInitial = ReceiveInteger(tag_paramIntInit);

        // Create the parameter
        m_pList[positionCounter] = new IntParamWorker(paramName, paramInitial);

        // Update the name list
        m_ParamNameList[positionCounter] = new char[strlen(paramName) + 1];
        strcpy(m_ParamNameList[positionCounter], paramName);

        // Update the exclusion list
        m_pExcl[positionCounter] = NULL; // No values are assumed in the exclude list

        // Increment the position counter
        positionCounter++;
    }
     
    // Tied parameters
    // Get the number of tied parameters
    int numberOfTiedParameters = ReceiveInteger(tag_paramTotalReal);

    // Request the values from the primary worker
    for (int entryTied = 0; entryTied < numberOfTiedParameters; entryTied++) {

        // Receive the values from the primary worker
        std::string tempName = ReceiveString(tag_paramRealName);
        IroncladString paramName = &tempName[0];

        // Recive the value
        double paramInitial = ReceiveDouble(tag_paramRealInit);

        // Recive the format
        std::string tempFmt = ReceiveString(tag_paramRealFmt);
        IroncladString paramFmt = &tempFmt[0];

        // Create a real parameter instead of the tied parameter
        m_pList[positionCounter] = new RealParamWorker(paramName, paramInitial, paramFmt);

        // Update the name list
        m_ParamNameList[positionCounter] = new char[strlen(paramName) + 1];
        strcpy(m_ParamNameList[positionCounter], paramName);

        // Update the exclusion list
        m_pExcl[positionCounter] = NULL; // No values are assumed in the exclude list

        // Increment the position counter
        positionCounter++;
    }

    // Geometric parameters




    // Set the arrays into a group object
    // Initialize a parameter group
    paramGroup = new ParameterGroupWorker();
    
    // Set the parameter values into the group
    paramGroup->SetGroupValues(m_pList, m_pExcl, m_ParamNameList, numberOfTotalParameters, numberOfExcluded);

    /*
    --------------------------------------------------------------------------------------------------------------------------
    Check template files against parameters, each parameter should appear in at least one template file or at least one 
    database entry.
    --------------------------------------------------------------------------------------------------------------------------
    */
    paramGroup->CheckTemplateFiles(m_filePairs, m_workerDirectory);

    /*
    --------------------------------------------------------------------------------------------------------------------------
    Check parameters for uniqueness, each parameter should be unique and should not be a substring of another parameter.
    --------------------------------------------------------------------------------------------------------------------------
    */
    paramGroup->CheckMnemonics();

}

/**************************************************************************************************************************************************************
ReceiveString

Receives a string from the primary worker using the MPI inerface
**************************************************************************************************************************************************************/
std::string  ModelWorker::ReceiveString(int tag_number) {

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

    // Remove the last character that's the null 
    convertedArray.pop_back();

    // Return the array to the calling function
    return convertedArray;
}

/**************************************************************************************************************************************************************
ReceiveInteger

Receives an integer from the primary worker using the MPI inerface
**************************************************************************************************************************************************************/
int  ModelWorker::ReceiveInteger(int tag_number) {

    // Initialize the variable to fill
    int value;

    // Request the value
    MPI_Recv(&value, 1, MPI_INT, 0, tag_number, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    return value;
}

/**************************************************************************************************************************************************************
ReceiveDouble

Receives a double from the primary worker using the MPI inerface
**************************************************************************************************************************************************************/
double  ModelWorker::ReceiveDouble(int tag_number) {

    // Initialize the variable to fill
    double value;

    // Request the value
    MPI_Recv(&value, 1, MPI_DOUBLE, 0, tag_number, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    return value;
}

/**************************************************************************************************************************************************************
ReceiveChar

Receives a character from the primary worker using the MPI inerface
**************************************************************************************************************************************************************/
char  ModelWorker::ReceiveChar(int tag_number) {

    // Initialize the variable to fill
    char value;

    // Request the value
    MPI_Recv(&value, 1, MPI_CHAR, 0, tag_number, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    return value;
}

/**************************************************************************************************************************************************************
SendDouble

Sends a double from the primary worker using the MPI inerface
**************************************************************************************************************************************************************/
void  ModelWorker::SendDouble(int tag_number, double value) {
    // TODO: Doc string

    // Send the value
    MPI_Ssend(&value, 1, MPI_DOUBLE, 0, tag_number, MPI_COMM_WORLD);

}

/**************************************************************************************************************************************************************
SendInteger

Sends an integer from the primary worker using the MPI inerface
**************************************************************************************************************************************************************/
void  ModelWorker::SendInteger(int tag_number, int value) {

    // Send the value
    MPI_Ssend(&value, 1, MPI_INT, 0, tag_number, MPI_COMM_WORLD);

}

/**************************************************************************************************************************************************************
SetupWorker

Sets up the worker based on the configuration information and preparis it for subsequent use by optimization algorithms.
**************************************************************************************************************************************************************/
void ModelWorker::SetupWorker(void) {

    // Create the working directory
    if (!std::filesystem::exists(m_workerDirectory)) {
        // Directory does not already exist. Create it.
        std::filesystem::create_directories(m_workerDirectory);
    }

    // Create holders to store the additional files added in the analysis
    std::vector<std::filesystem::path> additionalFiles;
    std::vector<int> indices2delete;

    // Copy all the files to it from the source directory
    for (int entryFile = 0; entryFile < m_fileCleanupList.size(); entryFile++) {

        // Get the current working directory
        std::filesystem::path workerDirectoryPath = m_workerDirectory;
	
	    // Create the source path. Remove the last character that creates issues.
	    std::filesystem::path sourcePath = m_fileCleanupList[entryFile];

        // Create an error code object to capture any errors
        std::error_code code;

        // Copy from the source path to the worker path. This handles directories and files differently due to the need to create the folder structure
        if (std::filesystem::is_directory(sourcePath)) {
            // Recursively copy files in the directory
            // Create the recursive iterator
            std::filesystem::recursive_directory_iterator iterator = std::filesystem::recursive_directory_iterator(sourcePath);

            // Iteratate over the files
            for (std::filesystem::path sourcePathRecursive : iterator) {
                
                if (!std::filesystem::is_directory(sourcePathRecursive)) {
                    // Copy the file
                    std::filesystem::path filePath = workerDirectoryPath; 
                    filePath /= sourcePathRecursive;

                    // Create directories to make sure the file can be transferred
                    if (!std::filesystem::is_directory(filePath.remove_filename())) {
                        std::filesystem::create_directories(filePath.remove_filename());
                    }

                    // Copy the file
		            std::filesystem::copy(sourcePathRecursive, filePath, code);

                    // Add to the cleanup list
                    additionalFiles.push_back(sourcePathRecursive.string());
                    indices2delete.push_back(entryFile);

                    // Handle the error code
                    if (code.value() == 2) {
                        std::cout << "Extra file " << sourcePathRecursive << "\t" << filePath << " is not found." << std::endl;
                    }
                }                 
            }

        } else {
            // It is a file. Copy it directly
            std::filesystem::path filePath = workerDirectoryPath /= sourcePath;

	        // Make sure the directory exists to the the files
	        if (!std::filesystem::is_directory(filePath.remove_filename())){
            	std::filesystem::create_directories(filePath.remove_filename());
	        }

	        // Perform the copy
            std::filesystem::copy(sourcePath, filePath, code);

            // Handle an error code
            if (code.value() == 2) {
                std::cout << "Extra file " << std::filesystem::current_path() << "\t" << sourcePath << "\t" << filePath << " is not found." << std::endl;
            }
        }
    }

    // Cleanup the initial list by removing directory names
    // Remove duplicates from the file cleanup index list
    std::sort(indices2delete.begin(), indices2delete.end());
    std::vector<int>::iterator ip = std::unique(indices2delete.begin(), indices2delete.end());

    // Resizing the vector so as to remove the undefined terms
    indices2delete.resize(std::distance(indices2delete.begin(), ip));

    // Starting from the end, delete entries from the file list
    for (int entry = indices2delete.size() - 1; entry >= 0; entry--) {
        // Remove the current entry from the file list
        m_fileCleanupList.erase(m_fileCleanupList.begin() + indices2delete[entry]);
    }

    // Add the list of additional files into the file cleanup list
    for (int entry = 0; entry < additionalFiles.size(); entry++) {
        m_fileCleanupList.push_back(additionalFiles[entry].string());
    }

    // Check that the vector only has unique values
    // Remove duplicates from the file cleanup index list
    std::sort(m_fileCleanupList.begin(), m_fileCleanupList.end());

    indices2delete.clear();
    if (m_fileCleanupList.size() > 1) {
        for (int entry = 1; entry < m_fileCleanupList.size(); entry++) {
            if (m_fileCleanupList[entry] == m_fileCleanupList[entry - 1]) {
                indices2delete.push_back(entry);
            }
        }
    }

    for (int entry = indices2delete.size() - 1; entry >= 0; entry--) {
        // Remove the current entry from the file list
        m_fileCleanupList.erase(m_fileCleanupList.begin() + indices2delete[entry]);
    }

}

/**************************************************************************************************************************************************************
CommenceMPIWork

Solves a model using the MPI workflow to request and send parameters
**************************************************************************************************************************************************************/
void  ModelWorker::CommenceMPIWork() {

    // Set the MPI variables
    MPI_Status mpiStatus;
    MPI_Request mpiRequest;
    int requestFlag = 0;

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
                double objective;
                try {
                    objective = ExecuteSingle();
                } catch (...) {
                    objective = INFINITY;
                }

                // Stack the alternative and the objective into a single array to transimit 
                double returnArray[2]{ (double)alternativeIndex, objective };
                
                // Send to the primary worker
                MPI_Ssend(returnArray, 2, MPI_DOUBLE, 0, tag_data, MPI_COMM_WORLD);

                // Check for best archive
                if (preserveModelBest) {
                    // Request archive instructions from the primary worker
                    int preserveValue = ReceiveInteger(tag_preserveBest);

                    // Preserve the model if requested
                    if (preserveValue == 1) {
                        // Move the files
                        PreserveBestModel(true);

                        // Send the preserve is complete
                        SendInteger(tag_preserveBest, 1);
                    }
                }

            } else {
                // Calculate the objective function
                // TODO: Add multiobject evaluation
            }

            // Archive the model if it is enabled. Disable saving to the best directory.
            if (preserveModelOutput) {
                PreserveArchiveModel();
            }
        }
    }
}

/**************************************************************************************************************************************************************
TerminateMPIWork

Terminates a secondary worker called using MPI
**************************************************************************************************************************************************************/
void  ModelWorker::TerminateMPIWork(void) {
    // Todo: update when additional operations are required at work termination

}

/**************************************************************************************************************************************************************
RequestParameters

Requests parameters from the primary worker when the seconday worker is using MPI
**************************************************************************************************************************************************************/
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
    delete parameters;

    // Return the alternative index to the calling function
    return alternativeIndex;
}

/**************************************************************************************************************************************************************
RequestContinue

Requests whether to continue from the primary worker when the seconday worker is using MPI
**************************************************************************************************************************************************************/
bool ModelWorker::RequestContinue(void) {

    // Setup the request space
    MPI_Message mpiMessage;
    MPI_Status mpiStatus;
    MPI_Request mpiRequest;
    int requestFlag = 0;

    // Send that it's reddy to receive work
    int sendValue = 1;
    MPI_Isend(&sendValue, 1, MPI_INT, 0, tag_continue, MPI_COMM_WORLD, &mpiRequest);
    MPI_Wait(&mpiRequest, &mpiStatus);

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

/**************************************************************************************************************************************************************
SetWorkerDirectory

Sets the worker directory into the worker when using solve on primary
**************************************************************************************************************************************************************/
void ModelWorker::SetWorkerDirectory(std::string workerDirectory) {

    // Set the string as the directory
    m_workerDirectory = workerDirectory;

    // Append the worker number onto the tag
    m_workerDirectory += std::to_string(rank);
    
};

/**************************************************************************************************************************************************************
SetWorkerSolveCommand

Sets the solve command into the worker when using solve on primary
**************************************************************************************************************************************************************/
void ModelWorker::SetWorkerSolveCommand(std::string solveCommandInput) {

    // Set the string as the solve command
    solveCommand = solveCommandInput;

};

/**************************************************************************************************************************************************************
SetWorkerPreserveArchiveCommand

Sets the archive command into the worker when using solve on primary
**************************************************************************************************************************************************************/
void ModelWorker::SetWorkerPreserveArchiveCommand(int preserveStatus, std::string preserveCommand) {

    // Setup the solution preservation
    if (preserveStatus == 0) {
        // Model will not be archived. Set the status to false
        preserveModelOutput = false;

    } else if (preserveStatus== 1) {
        // Model will be archived using standard method
        // Toggle the preserve status to true
        preserveModelOutput = true;

    } else if (preserveStatus == 2) {
        // Model will be archived using a custom command.
        // Toggle the preserve status to true
        preserveModelOutput = true;

        // Set the preservation command
        preserveOutputCommand = preserveCommand;

    }
}

/**************************************************************************************************************************************************************
SetWorkerPreserveBestCommand

Sets the preserve command into the worker when using solve on primary
**************************************************************************************************************************************************************/
void ModelWorker::SetWorkerPreserveBestCommand(int preserveBest, std::string bestCommand) {

    // Setup the best perservation
    if (preserveBest == 0) {
        // Model will not be archived. Set the status to false
        preserveModelBest = false;

    } else if (preserveBest == 1) {
        // Model will be archived using standard method
        // Toggle the preserve status to true
        preserveModelBest = true;

    } else if (preserveBest == 2) {
        // Model will be archived using a custom command.
        // Toggle the preserve status to true
        preserveModelBest = true;

        // Set the pservation command
        preserveBestCommand = bestCommand;

    }
};

/**************************************************************************************************************************************************************
SetWorkerExtraFiles

Sets the extra files into the worker when using solve on primary
**************************************************************************************************************************************************************/
void ModelWorker::SetWorkerExtraFiles(std::vector<std::string> extraFiles) {

    // Set the file cleanup list as the extra file list
    m_fileCleanupList = extraFiles;

};

/**************************************************************************************************************************************************************
SetWorkerFilePairs

Sets the file pairs into the worker when using solve on primary
**************************************************************************************************************************************************************/
void ModelWorker::SetWorkerFilePairs(std::vector<std::vector<std::string>> filePairs) {

    // Set the file pairs as object
    m_filePairs = filePairs;

};

/**************************************************************************************************************************************************************
SetWorkerObservations

Sets the observation group into the worker when using solve on primary
**************************************************************************************************************************************************************/
void ModelWorker::SetWorkerObservations(ObservationGroup* m_pObsGroup) {

    // Construct the group from the value passed in
    observationGroup = m_pObsGroup;

};

/**************************************************************************************************************************************************************
SetWorkerParameters

Sets the parameters initially into the worker when using solve on primary
**************************************************************************************************************************************************************/
void ModelWorker::SetWorkerParameters(ParameterGroup* m_pParamGroup) {

    // Get the total number of paramters
    int numberOfTotalParameters = m_pParamGroup->GetNumParams() + m_pParamGroup->GetNumTiedParams();

    // Create a counter for the number of excluded parameters
    int numberOfExcluded = 0;

    // Create a parameter group object to store everything in
    ParameterWorkerABC** m_pList = new ParameterWorkerABC * [numberOfTotalParameters];
    char** m_ParamNameList = new char* [numberOfTotalParameters];
    ParameterWorkerABC** m_pExcl = new ParameterWorkerABC * [numberOfTotalParameters];
    int positionCounter = 0;

    // Real parameters
    int numberOfRealParameters = m_pParamGroup->GetNumRealParams();
    for (int entryParameter = 0; entryParameter < numberOfRealParameters; entryParameter++) {
        // Get the parameter from the group, casing to the correct type
        RealParam* param = (RealParam*)m_pParamGroup->GetParamPtr(entryParameter);

        // Get the information from the parameter
        std::string name = param->GetName();
        IroncladString paramName = &name[0];

        double sampleValue = param->GetEstimatedValueTransformed();
        sampleValue = param->ConvertOutVal(sampleValue);

        std::string fixFmt = param->GetFixFmt();
        IroncladString paramFmt = &fixFmt[0];

        // Create the parameter
        m_pList[positionCounter] = new RealParamWorker(paramName, sampleValue, paramFmt);

        // Construct the name list
        m_ParamNameList[positionCounter] = new char[strlen(paramName) + 1];
        strcpy(m_ParamNameList[positionCounter], paramName);

        // Update the exclusion list
        m_pExcl[positionCounter] = NULL; // No values are assumed in the exclude list

        // Increment the position counter
        positionCounter++;
    }

    // Integer parameters
    // Get the number of integer parameters
    int numberOfIntegerParameters = m_pParamGroup->GetNumIntParams();

    // Request the values from the primary worker
    for (int entryInt = 0; entryInt < numberOfIntegerParameters; entryInt++) {
        // Get the parameter from the group, casing to the correct type
        IntParam* param = (IntParam*)m_pParamGroup->GetParamPtr(positionCounter);

        // Get the information from the parameter
        std::string name = param->GetName();
        IroncladString paramName = &name[0];

        int sampleValue = param->GetEstimatedValueTransformed();
        sampleValue = param->ConvertOutVal(sampleValue);

        // Create the parameter
        m_pList[positionCounter] = new IntParamWorker(paramName, sampleValue);

        // Update the name list
        m_ParamNameList[positionCounter] = new char[strlen(paramName) + 1];
        strcpy(m_ParamNameList[positionCounter], paramName);

        // Update the exclusion list
        m_pExcl[positionCounter] = NULL; // No values are assumed in the exclude list

        // Increment the position counter
        positionCounter++;
    }

    // Tied parameters
    // Get the number of tied parameters
    int numberOfTiedParameters = m_pParamGroup->GetNumTiedParams();

    // Request the values from the primary worker
    for (int entryTied = 0; entryTied < numberOfTiedParameters; entryTied++) {

        // Get the parameter from the group
        TiedParamABC* param = m_pParamGroup->GetTiedParamPtr(entryTied);

        // Get the information from the parameter
        std::string name = param->GetName();
        IroncladString paramName = &name[0];

        double sampleValue = param->GetEstimatedValueTransformed();
        
        std::string fixFmt = param->GetFixFmt();
        IroncladString paramFmt = &fixFmt[0];

        // Create a real parameter instead of the tied parameter
        m_pList[positionCounter] = new RealParamWorker(paramName, sampleValue, paramFmt);

        // Update the name list
        m_ParamNameList[positionCounter] = new char[strlen(paramName) + 1];
        strcpy(m_ParamNameList[positionCounter], paramName);

        // Update the exclusion list
        m_pExcl[positionCounter] = NULL; // No values are assumed in the exclude list

        // Increment the position counter
        positionCounter++;
    }

    // Geometric parameters




    // Set the arrays into a group object
    // Initialize a parameter group
    paramGroup = new ParameterGroupWorker();

    // Set the parameter values into the group
    paramGroup->SetGroupValues(m_pList, m_pExcl, m_ParamNameList, numberOfTotalParameters, numberOfExcluded);


    /*
    --------------------------------------------------------------------------------------------------------------------------
    Check template files against parameters, each parameter should appear in at least one template file or at least one database entry.
    --------------------------------------------------------------------------------------------------------------------------
    */
    paramGroup->CheckTemplateFiles(m_filePairs, m_workerDirectory);

    /*
    --------------------------------------------------------------------------------------------------------------------------
    Check parameters for uniqueness, each parameter should be unique and should not be a substring of another parameter.
    --------------------------------------------------------------------------------------------------------------------------
    */
    paramGroup->CheckMnemonics();

};

/**************************************************************************************************************************************************************
SetStandardParameters

Sets the parameters prior to solve into the worker when using solve on primary
**************************************************************************************************************************************************************/
void ModelWorker::SetStandardParameters(std::vector<double> inputParameters) {

    // Set the parameters into the the parameter group
    paramGroup->WriteParams(&inputParameters[0]);

};


/**************************************************************************************************************************************************************
PreserveBestModel()

Preserves the best model solution by moving it from the working directory.
**************************************************************************************************************************************************************/
void ModelWorker::PreserveBestModel(bool preserveBest) {

    // Model preservation is disabled. Take no action.
    if (preserveModelOutput == false && preserveBest == false) {
        return;
    }


    // Preserve the best model output
    if (preserveBestCommand.empty()) {
        // Get the current path to the worker
        std::filesystem::path directory = std::filesystem::current_path();

        // Construct the archive folder
        directory /= std::string("archive");
        if (!std::filesystem::exists(directory)) {
            std::filesystem::create_directories(directory);
        };

        // Construct the worker folder
        std::filesystem::path workerDirectoryPath = m_workerDirectory;

        // Append the best into the archive path
        directory /= "best";

        // Delete the current contents of the best directory
        std::filesystem::remove_all(directory);

        // Create the directory if it doesn't exist
        if (!std::filesystem::exists(directory)) {
            std::filesystem::create_directories(directory);
        };

        // Move the run
        MoveModel(m_workerDirectory, directory);

    } else {
        // Run user commands to save specific model files
        // Construct the command
        std::string preserveCommandActive = preserveBestCommand;
        preserveCommandActive.append(' ' + std::to_string(rank));
        preserveCommandActive.append(' ' + std::to_string(solveCounter));
        std::cout << preserveCommandActive << std::endl;

        // Send it to the system to execute
        system(preserveCommandActive.c_str());

    }
}


/**************************************************************************************************************************************************************
PreserveArchiveModel()

Preserves the model solution by moving it from the working directory.
**************************************************************************************************************************************************************/
void ModelWorker::PreserveArchiveModel(){

    // Preserve all model output
    if (preserveOutputCommand.empty()) {
        // Get the current path to the worker
        std::filesystem::path directory = std::filesystem::current_path();

        // Construct the archive folder
        directory /= std::string("archive");
        if (!std::filesystem::exists(directory)) {
            std::filesystem::create_directories(directory);
        };

        // Construct the worker folder
        directory /= m_workerDirectory;

        if (!std::filesystem::exists(directory)) {
            std::filesystem::create_directory(directory);
        }

        // Construct the run folder 
        directory /= "run_" + std::to_string(solveCounter - 1); // This adjusts by one to align with the model output files
        std::filesystem::create_directory(directory);

        // Move the run
        MoveModel(m_workerDirectory, directory);

    } else {
        // Preserve to the archive folder for the run
        // Construct the command
        std::string preserveCommandActive = preserveOutputCommand;
        preserveCommandActive.append(' ' + std::to_string(rank));
        preserveCommandActive.append(' ' + std::to_string(solveCounter));
        std::cout << preserveCommandActive << std::endl;

        // Send it to the system to execute
        system(preserveCommandActive.c_str());
    }
}

/**************************************************************************************************************************************************************
MoveModel()

Moves the solution output to a target directory
**************************************************************************************************************************************************************/
void ModelWorker::MoveModel(std::filesystem::path workerDirectoryPath, std::filesystem::path archiveDestinationPath) {

    // Move the files from the worker to the archive
        // Create the iterator over the files in the working directory
    std::filesystem::recursive_directory_iterator iterator = std::filesystem::recursive_directory_iterator(workerDirectoryPath);

    // Iterate over the directory contents
    for (std::filesystem::path path : iterator) {
        // Confirm that the path is to a file
        if (!std::filesystem::is_directory(path)) {
            // Check that the file is not an input file
            bool archiveFlag = true;

            for (int entry = 0; entry < m_fileCleanupList.size(); entry++) {
                // Construct the path to the file in the worker directory
                std::filesystem::path temp = m_workerDirectory;
                temp /= m_fileCleanupList[entry];

                // Compare to the active path to the cleanup list value
                if (std::filesystem::equivalent(temp, path)) {
                    archiveFlag = false;
                    break;
                }
            }

            // If the file is not in the cleanup list, archive it
            if (archiveFlag) {
                // 
                std::filesystem::path copyPath = archiveDestinationPath;
                for (auto it = path.begin(); it != path.end(); ++it) {
                    if (*it != workerDirectoryPath) {
                        copyPath /= *it;
                    }
                }

                // Copy the file
                std::filesystem::create_directories(copyPath.parent_path());
                std::filesystem::copy(path, copyPath);
            }
        }
    }
}


/**************************************************************************************************************************************************************
ExecuteSingle()

Executes the model (or surrogate) and returns the objective function value.
**************************************************************************************************************************************************************/
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

/**************************************************************************************************************************************************************
Execute()

Executes the model returns a vector of objective function values.
**************************************************************************************************************************************************************/
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
    m_NumSolves++;

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
    PreserveModel(rank, GetTrialNumber(), m_NumSolves, pCatStr);

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

/**************************************************************************************************************************************************************
StdExecute()

Executes the standard model and returns the objective function value.
**************************************************************************************************************************************************************/
double ModelWorker::StdExecute(double viol) {
    //IroncladString dirName = &workerDirectory[0];

    bool isGoodTopo;

    //inc. number of times model has been executed
    solveCounter++;

    // Make substitution of parameters into model input files
    for (int entryPair = 0; entryPair < m_filePairs.size(); entryPair++) {
        // Setup the input file and path
        std::string inTemp = m_filePairs[entryPair][0];
        IroncladString inFile = &inTemp[0];

        // Setup the output file and path
        std::string outFileString = m_filePairs[entryPair][1];
        std::filesystem::path outTemp = m_workerDirectory;
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
    std::filesystem::current_path(m_workerDirectory);


    //make substitution of parameters into model input databases
    /*if (m_DbaseList != NULL)
    {
        m_pParamGroup->SubIntoDbase(m_DbaseList);
    } */

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



/**************************************************************************************************************************************************************
DisklessExecute()
   Executes an inernal model without using I/O.
**************************************************************************************************************************************************************/
/*double ModelWorker::DisklessExecute(void) {
    //inc. number of times model has been executed
    m_NumSolves++;

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


/**************************************************************************************************************************************************************
Write()

Store parameter and objective function value to model output file.
**************************************************************************************************************************************************************/
void ModelWorker::Write(double objFuncVal) {
    //ResponseVarGroup* pRespVarGroup;
    FILE* pFile;

    /*pRespVarGroup = NULL;
    if (m_pObjFunc != NULL)
    {
        pRespVarGroup = (ResponseVarGroup*)(m_pObjFunc->GetResponseVarGroup());
    }*/

    // Create the name of the log file for the secondary worker
    std::filesystem::path parentPath = m_workerDirectory;
    parentPath /= "..";
    parentPath /= "..";
    parentPath /= "OstModel" + std::to_string(rank) + ".txt";
    parentPath = parentPath.lexically_normal();

    // Log the solve
    pFile = fopen(parentPath.string().c_str(), "a+");
    if (pFile == NULL) {
        std::cout << "File write error" << std::endl;
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
    }

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
} 

