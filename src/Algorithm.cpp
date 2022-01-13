/**************************************************************************************************************************************************************
The algorithm class is the base class from which all other algorithms inheirit. It provides core functionality shared across all other algorithms.

**************************************************************************************************************************************************************/


#include "Algorithm.h"


#define JOB_SUCCEDDED (0)
#define JOB_FAILED    (1)
#define JOB_TIMED_OUT (2)

//objective function categories passed to PreserveModel output script
IroncladString ObjFuncBest = "best";
IroncladString ObjFuncBehavioral = "behavioral";
IroncladString ObjFuncNonBehavioral = "non-behavioral";
IroncladString ObjFuncDominated = "dominated";
IroncladString ObjFuncNonDominated = "non-dominated";
IroncladString ObjFuncOther = "other";


/**************************************************************************************************************************************************************
 default CTOR
**************************************************************************************************************************************************************/
Algorithm::Algorithm(void) {
    FilePair* pFilePair;
    char* line;
    char tmp1[DEF_STR_SZ];
    char tmp2[DEF_STR_SZ];
    char tmp3[DEF_STR_SZ];
    FILE* pInFile;
    int id, i, j;
    bool quoteWrap; //if true, then executable must be wrapped in quotes


    // Get the current directory and the input filename
    IroncladString inFileName = GetInFileName();
    UnmoveableString pDirName = GetExeDirName();

    // Open the input file
    pInFile = fopen(inFileName, "r");
    
    // Check that the input file is available
    if (pInFile == NULL) {
        FileOpenFailure("Algorithm::CTOR", inFileName);
    }

    /* 
    --------------------------------------------------------------------------------------------------------------------------
    Verify the structure of the input file
    --------------------------------------------------------------------------------------------------------------------------
    */
    // Check the input file for critical entries which have no reasonable defaults
    FindToken(pInFile, "BeginFilePairs", inFileName);
    FindToken(pInFile, "EndFilePairs", inFileName);
    rewind(pInFile);
    FindToken(pInFile, "ModelExecutable", inFileName);

    /*
    --------------------------------------------------------------------------------------------------------------------------
    Read in and create the directory from which the model will be run.
    --------------------------------------------------------------------------------------------------------------------------
    */
    rewind(pInFile);
    if ((CheckToken(pInFile, "ModelSubdir", inFileName) == true) ||
        (CheckToken(pInFile, "ModelSubDir", inFileName) == true)) {
        line = GetCurDataLine();
        MyTrim(line);
        if (strlen(line) < 12) {
            LogError(ERR_IN_PARSE, "Bad ModelSubdir");
            ExitProgram(1);
        }
        strcpy(tmp1, &line[11]);
        
        //strip whitespace
        MyTrim(tmp1);

        //strip quotes
        if (tmp1[0] == '"') { tmp1[0] = ' '; }
        if (tmp1[strlen(tmp1) - 1] == '"') { tmp1[strlen(tmp1) - 1] = ' '; }
        MyTrim(tmp1);
        
        strcpy(pDirName, tmp1);
        strcpy(m_DirPrefix, tmp1);
    }


    /*
    --------------------------------------------------------------------------------------------------------------------------
    Setup the executable to be called during the OSTRICH run
    --------------------------------------------------------------------------------------------------------------------------
    */
    rewind(pInFile);
    FindToken(pInFile, "ModelExecutable", inFileName);
    line = GetCurDataLine();

    i = ExtractString(line, tmp2);
    i = ValidateExtraction(i, 1, 1, "Model()");
    i = ExtractFileName(&(line[i]), tmp1);

    int len = (int)strlen(tmp1) + 1;
    m_ExecCmd = new char[len];
    MEM_CHECK(m_ExecCmd);

    strcpy(m_ExecCmd, tmp1);


    /*
    --------------------------------------------------------------------------------------------------------------------------
    Read in the 'File Pairs': a set of template files and their model equivalents.
    --------------------------------------------------------------------------------------------------------------------------
    */
    rewind(pInFile);
    FindToken(pInFile, "BeginFilePairs", inFileName);
    line = GetNxtDataLine(pInFile, inFileName);
    while (strstr(line, "EndFilePairs") == NULL) {
        if ((strstr(line, ";") == NULL) && (strstr(line, "\t") == NULL)) {
            LogError(ERR_FILE_IO, "Model::CTOR(): missing separator (;) in file pair.");
        }

        /*
        ----------------------------------------------------------------------------------------------------------------------
        Read in file pairs, taking care to preserve full path, even in the presence of long and space-separated filenames.
        ----------------------------------------------------------------------------------------------------------------------
        */
        //first file in pair.....
        i = ExtractFileName(line, tmp1);
        i = ExtractFileName(&(line[i]), tmp2);

        // Push back into the vector of file pairs
        std::vector<std::string> filePair{ tmp1, tmp2 };
        fileListPairs.push_back(filePair);

        // Get the next line of data
        line = GetNxtDataLine(pInFile, inFileName);
    }

    /*
    --------------------------------------------------------------------------------------------------------------------------
    Read in any extra model files, these will need to be copied to the model subdirectory.
    --------------------------------------------------------------------------------------------------------------------------
    */
    rewind(pInFile);
    if (CheckToken(pInFile, "BeginExtraFiles", inFileName) == true) {
        //make sure end token exists
        FindToken(pInFile, "EndExtraFiles", inFileName);
        rewind(pInFile);
        FindToken(pInFile, "BeginExtraFiles", inFileName);

        line = GetNxtDataLine(pInFile, inFileName);
        while (strstr(line, "EndExtraFiles") == NULL) {
            //extra file
            ExtractFileName(line, tmp1);

            // add to cleanup list
            if (m_pFileCleanupList == NULL) {
                m_pFileCleanupList = new FileList(tmp1);
            }
            else {
                m_pFileCleanupList->Insert(tmp1);
            }
            
            // Get the next data line
            line = GetNxtDataLine(pInFile, inFileName);
        }
    }

    /*
    --------------------------------------------------------------------------------------------------------------------------
    Read in any extra model folders, these will need to be copied (along with their contects) to the model subdirectory.
    --------------------------------------------------------------------------------------------------------------------------
    */
    rewind(pInFile);
    if (CheckToken(pInFile, "BeginExtraDirs", inFileName) == true) {
        //make sure end token exists
        FindToken(pInFile, "EndExtraDirs", inFileName);
        rewind(pInFile);
        FindToken(pInFile, "BeginExtraDirs", inFileName);

        line = GetNxtDataLine(pInFile, inFileName);
        while (strcmp(line, "EndExtraDirs") != 0) {
            //extra dir
            ExtractFileName(line, tmp1);

            // Move all entries in the directory
            if (pDirName[0] != '.') {
                // Construct the path to the source source
                std::filesystem::path sourcePath = std::filesystem::current_path() /= tmp1;

                // Get a list of files in the source directory
                std::vector<std::string> sourceFiles;
                for (const auto& entry : std::filesystem::recursive_directory_iterator(sourcePath)) {
                    if (!std::filesystem::is_directory(entry.path())) {
                        std::filesystem::path partialPath = entry.path().lexically_relative(std::filesystem::current_path());
                        sourceFiles.push_back(partialPath.string());
                    }
                }

                // Add files within the directory to the cleanup list
                for (int entry = 0; entry < sourceFiles.size(); entry++) {
                    // add to cleanup list
                    if (m_pFileCleanupList == NULL) {
                        m_pFileCleanupList = new FileList(&sourceFiles[entry][0]);
                    }
                    else {
                        m_pFileCleanupList->Insert(&sourceFiles[entry][0]);
                    }
                }
            }

            // Get the next data line
            line = GetNxtDataLine(pInFile, inFileName);
        }
    }

    /*
    --------------------------------------------------------------------------------------------------------------------------
    Read in diskless model
    --------------------------------------------------------------------------------------------------------------------------
    */
    // todo: reenable this functionality
    rewind(pInFile);
    if (CheckToken(pInFile, "DisklessModel", inFileName) == true) {
        line = GetCurDataLine();
        sscanf(line, "%s %s", tmp1, tmp2);
        MyStrLwr(tmp2);
        if (strncmp(tmp2, "yes", 3) == 0) { m_bDiskless = true; }
    }

    /*
    --------------------------------------------------------------------------------------------------------------------------
    Read in solve on primary options
    --------------------------------------------------------------------------------------------------------------------------
    */
    rewind(pInFile);
    if (CheckToken(pInFile, "SolveOnPrimary", inFileName) == true) {
        line = GetCurDataLine();
        sscanf(line, "%s %s", tmp1, tmp2);
        MyStrLwr(tmp2);
        if (strncmp(tmp2, "yes", 3) == 0) { m_bSolveOnPrimary = true; }
    }

    /*
    --------------------------------------------------------------------------------------------------------------------------
    Read in the name of the batch file to run whenever a new 'best' solution is found.
    --------------------------------------------------------------------------------------------------------------------------
    */
    rewind(pInFile);
    if (CheckToken(pInFile, "PreserveBestModel", inFileName) == true) {
        line = GetCurDataLine();
        //line format = 'PreserveBestModel   <var>'
        /*--------------------------------------------------------------------------------------------------------------------
        Read in executable, taking care to preserve full path, even in the presence of long and space-separated filenames.
        ----------------------------------------------------------------------------------------------------------------------
        */
        line = GetCurDataLine();
        sscanf(line, "%s %s", tmp1, tmp2);
        MyStrLwr(tmp2);
        MyTrim(tmp2);

        if (strcmp(tmp2, "yes") == 0) {
            m_bPreserveModelBest = true;
        }
        else if (strcmp(tmp2, "no") == 0) {
            m_bPreserveModelBest = false;
        }
        else {
            // Convert the input command to a path
            std::filesystem::path preserveInput = tmp2;

            // Check that the script exists
            if (!std::filesystem::exists(preserveInput)) {
                sprintf(tmp1, "File for preserving model best output (|%s|) not found", tmp2);
                std::cout << tmp1 << std::endl;
                LogError(ERR_FILE_IO, tmp1);
                ExitProgram(1);
            }

            // Set the variables into the class fields
            m_bPreserveModelBest = true;
            m_PreserveBestCmd = preserveInput;
        }
    }

    /*
    --------------------------------------------------------------------------------------------------------------------------
    Check for alternate objective function (default is WSSE)
    --------------------------------------------------------------------------------------------------------------------------
    */
    // todo: reenable this functionality
    rewind(pInFile);
    if (CheckToken(pInFile, "ObjectiveFunction", inFileName) == true) {
        line = GetCurDataLine();
        sscanf(line, "%s %s", tmp1, tmp2);
        MyStrLwr(tmp2);
        if (strstr(tmp2, "user") != NULL) { m_ObjFuncId = OBJ_FUNC_USER; }
        else if (strstr(tmp2, "sawe") != NULL) { m_ObjFuncId = OBJ_FUNC_SAWE; }
        else if (strstr(tmp2, "wsse") != NULL) { m_ObjFuncId = OBJ_FUNC_WSSE; }
        else if (strstr(tmp2, "pato") != NULL) { m_ObjFuncId = OBJ_FUNC_PATO; }
        else if (strstr(tmp2, "gcop") != NULL) { m_ObjFuncId = OBJ_FUNC_GCOP; }
    } 

    /*
    --------------------------------------------------------------------------------------------------------------------------
    Read in flag to check sensitivities
    --------------------------------------------------------------------------------------------------------------------------
    */
    // todo: reenable this functionality
    /*rewind(pInFile);
    if (CheckToken(pInFile, "CheckSensitivities", inFileName) == true) {
        line = GetCurDataLine();
        sscanf(line, "%s %s", tmp1, tmp2);
        MyStrLwr(tmp2);
        if (strncmp(tmp2, "yes", 3) == 0) { m_bCheckGlobalSens = true; }
    }*/

    /*
    --------------------------------------------------------------------------------------------------------------------------
    Read in flag to use surrogate models
    --------------------------------------------------------------------------------------------------------------------------
    */
    // todo: reenable this functionality
    /*rewind(pInFile);
    if (CheckToken(pInFile, "SurrogateApproach", inFileName) == true) {
        line = GetCurDataLine();
        sscanf(line, "%s %s", tmp1, tmp2);
        MyStrLwr(tmp2);
        if (strncmp(tmp2, "yes", 3) == 0) { m_bUseSurrogates = true; }
    }*/

    /*
    --------------------------------------------------------------------------------------------------------------------------
    Read in flag to use SuperMUSE
    --------------------------------------------------------------------------------------------------------------------------
    */
    // todo: reenable this functionality
    /*rewind(pInFile);
    bool bSMUSE = false;
    if (CheckToken(pInFile, "SuperMUSE", inFileName) == true) {
        line = GetCurDataLine();
        sscanf(line, "%s %s", tmp1, tmp2);
        MyStrLwr(tmp2);
        if (strncmp(tmp2, "yes", 3) == 0) {
            EnableSuperMUSE();
            InitSuperMUSE(pInFile, (ModelABC*)this);
            bSMUSE = true;
        }
    }*/

    /*
    --------------------------------------------------------------------------------------------------------------------------
    Read in flag to preserve model output files. 
    --------------------------------------------------------------------------------------------------------------------------
    */
    rewind(pInFile);
    if (CheckToken(pInFile, "PreserveModelOutput", inFileName) == true) {
        line = GetCurDataLine();
        sscanf(line, "%s %s", tmp1, tmp2);
        MyStrLwr(tmp2);
        MyTrim(tmp2);

        if (strcmp(tmp2, "yes") == 0) {
            m_bPreserveModelOutput = true;

        } else if (strcmp(tmp2, "no") == 0) {
            m_bPreserveModelOutput = false;

        } else { 
            // Convert the input command to a path
            std::filesystem::path preserveInput = tmp2;

            // Check that the script exists
            if (!std::filesystem::exists(preserveInput)) {
                sprintf(tmp1, "File for preserving model output (|%s|) not found", tmp2);
                std::cout << tmp1 << std::endl;
                LogError(ERR_FILE_IO, tmp1);
                ExitProgram(1);
            }

            // Set the variables into the class fields
            m_bPreserveModelOutput = true;
            m_PreserveOutputCmd = preserveInput;
        } 
    }

    /*
    --------------------------------------------------------------------------------------------------------------------------
    Read in warm start flag. This only applies if we want to restart Ostrich from a previously aborted or interrupted run in 
    which the OstModel0.txt file has been preserved.
    --------------------------------------------------------------------------------------------------------------------------
    */
    // TODO: rework this process
    rewind(pInFile);
    if (CheckToken(pInFile, "WarmStart", inFileName) == true) {
        line = GetCurDataLine();
        sscanf(line, "%s %s", tmp1, tmp2);
        MyStrLwr(tmp2);
        if (strncmp(tmp2, "yes", 3) == 0) {
            std::cout << "Warm Start has been activated" << std::endl;;
            std::cout << "Ostrich will resume a pervious search" << std::endl;

            m_bWarmStart = true;
        }
    }

    /*
    --------------------------------------------------------------------------------------------------------------------------
    Read in caching flag. This only applies if we want to search the history of model evaluations (stored in OstModel0.txt) 
    prior to running the model. If the candidate solution has already been evaluated, it will be found in the cache and so we 
    won't have to run the model.
    --------------------------------------------------------------------------------------------------------------------------
    */
    rewind(pInFile);
    if (CheckToken(pInFile, "Caching", inFileName) == true) {
        line = GetCurDataLine();
        sscanf(line, "%s %s", tmp1, tmp2);
        MyStrLwr(tmp2);
        if (strncmp(tmp2, "yes", 3) == 0) {
            m_bCaching = true;
        }
    }

    /*
    --------------------------------------------------------------------------------------------------------------------------
    Read in number of digits of precision in I/O.
    --------------------------------------------------------------------------------------------------------------------------
    */
    rewind(pInFile);
    if (CheckToken(pInFile, "NumDigitsOfPrecision", inFileName) == true) {
        line = GetCurDataLine();
        sscanf(line, "%s %d", tmp1, &m_Precision);
        if ((m_Precision < 1) || (m_Precision > 32)) {
            LogError(ERR_FILE_IO, "Invalid precision setting - defaulting to 6 digits.");
            m_Precision = 6;
        }
    }

    /*
    --------------------------------------------------------------------------------------------------------------------------
    Read in telescoping strategy, if any.
    --------------------------------------------------------------------------------------------------------------------------
    */
    // todo: reenable this functionality
    /*char tscp[DEF_STR_SZ];
    rewind(pInFile);
    if (CheckToken(pInFile, "TelescopingStrategy", inFileName) == true) {
        line = GetCurDataLine();
        sscanf(line, "%s %s", tmp1, tscp);

        MyStrLwr(tscp);
        if (strcmp(tscp, "convex-power") == 0) m_Telescope = TSCOPE_PVEX;
        if (strcmp(tscp, "convex") == 0) m_Telescope = TSCOPE_CVEX;
        if (strcmp(tscp, "linear") == 0) m_Telescope = TSCOPE_LINR;
        if (strcmp(tscp, "concave") == 0) m_Telescope = TSCOPE_CAVE;
        if (strcmp(tscp, "delayed-concave") == 0) m_Telescope = TSCOPE_DCVE;
    }*/

    /*
    --------------------------------------------------------------------------------------------------------------------------
    Read in BoxCox transformation, if any.
    --------------------------------------------------------------------------------------------------------------------------
    */
    // todo: reenable this functionality
    bool bBoxCox = false;
    double boxCoxValue = 1.00;
    /*char boxcox[DEF_STR_SZ];
    rewind(pInFile);
    if (CheckToken(pInFile, "BoxCoxTransformation", inFileName) == true) {
        line = GetCurDataLine();
        sscanf(line, "%s %s", tmp1, boxcox);
        MyStrLwr(boxcox);
        if (strcmp(boxcox, "extract") == 0) {
            boxCoxValue = ExtractBoxCoxValue();
        }
        else {
            boxCoxValue = atof(boxcox);
        }
        bBoxCox = true;
    }

    fclose(pInFile);

    if (bSMUSE) CleanSuperMUSE();*/



    /*
    --------------------------------------------------------------------------------------------------------------------------
    Read in the observation group information
    --------------------------------------------------------------------------------------------------------------------------
    */
    NEW_PRINT("ParameterGroup", 1);
    m_pParamGroup = new ParameterGroup(true);
    MEM_CHECK(m_pParamGroup);

    //check bounds on parameter group
    m_pParamGroup->CheckBounds();

    if ((m_ObjFuncId == OBJ_FUNC_WSSE) || (m_ObjFuncId == OBJ_FUNC_SAWE)) {
        NEW_PRINT("ObservationGroup", 1);
        m_pObsGroup = new ObservationGroup();
        MEM_CHECK(m_pObsGroup);
    }

    //setup the objective function
    if (m_ObjFuncId == OBJ_FUNC_WSSE) {
        NEW_PRINT("WSSE", 1);
        m_pObjFunc = new WSSE(m_pObsGroup, bBoxCox, boxCoxValue);
    }
    else if (m_ObjFuncId == OBJ_FUNC_SAWE) {
        NEW_PRINT("SAWE", 1);
        m_pObjFunc = new SAWE(m_pObsGroup);
    }
    else if (m_ObjFuncId == OBJ_FUNC_PATO) {
        NEW_PRINT("PATO", 1);
        m_pObjFunc = new PATO(m_pParamGroup);
    }
    else if (m_ObjFuncId == OBJ_FUNC_GCOP) {
        NEW_PRINT("GCOP", 1);
        m_pObjFunc = new GCOP(m_pParamGroup);
    }
    else {/* OBJ_FUNC_USER) */ 
        NEW_PRINT("USER", 1);
        m_pObjFunc = new UserObjFunc(GetOstExeOut());
    }
    MEM_CHECK(m_pObjFunc);

    if (m_pObjFunc->CalcMultiObjFunc(NULL, -1) > 1) {
        m_bMultiObjProblem = true;
    }





    /*
    --------------------------------------------------------------------------------------------------------------------------
    Initialize surrogate models, if the surrogate-based approach is enabled.
    --------------------------------------------------------------------------------------------------------------------------
    */
    /*if (m_bUseSurrogates == true) {
        m_pDecision = new DecisionModule((ModelABC*)this);
    }

    CheckGlobalSensitivity();

    pInFile = fopen(inFileName, "r");
    if (pInFile == NULL) {
        FileOpenFailure("Model::CTOR", inFileName);
    }

    if (CheckToken(pInFile, "BeginParameterCorrection", inFileName) == true) {
        m_pParameterCorrection = new ParameterCorrection(m_pParamGroup);
    }
    fclose(pInFile);
    */

    IncCtorCount();

} 

/***************************************************************************************************************************************************************
Free up memory.
**************************************************************************************************************************************************************/
void Algorithm::Destroy(void) {
    delete m_pObsGroup;
    delete m_pParamGroup;
    delete m_pParameterCorrection;
    delete m_pObjFunc;
    delete[] m_ExecCmd;
    delete[] m_CurMultiObjF;
    m_bPreserveModelBest = false;
    delete m_pDecision;

    if (m_pFileCleanupList != NULL) {
        IroncladString dirName = GetExeDirName();
        if (dirName[0] != '.') {
            //m_pFileCleanupList->Cleanup(dirName, m_DirPrefix, 0);
        }
        delete m_pFileCleanupList;
    }

    //cleanup diskless data if needed
    DisklessIsotherm(NULL, NULL);
    DisklessMcCammon(NULL, NULL);

    IncDtorCount();
}


/**************************************************************************************************************************************************************
PerformParameterCorrections()
**************************************************************************************************************************************************************/
void Algorithm::PerformParameterCorrections(void) {
    if (m_pParameterCorrection != NULL) m_pParameterCorrection->Execute();
}/* end PerformParameterCorrections() */


/**************************************************************************************************************************************************************
GetObjectiveFunctionString()

Returns the objective function string indicating its type
**************************************************************************************************************************************************************/
UnchangeableString Algorithm::GetObjectiveFunctionString(void) {
    return m_pObjFunc->GetObjFuncStr();
}

/**************************************************************************************************************************************************************
GetObjFuncPtr()

Returns a pointer to the objective function.
**************************************************************************************************************************************************************/
ObjectiveFunction* Algorithm::GetObjFuncPtr(void) {
    return m_pObjFunc;
} 

/**************************************************************************************************************************************************************
GetObsGroupPtr()

Returns the observation group pointer
**************************************************************************************************************************************************************/
ObservationGroup* Algorithm::GetObsGroupPtr(void) {
    return m_pObsGroup;
}

/**************************************************************************************************************************************************************
GetParamGroupPtr()

Returns the parameter group pointer.
**************************************************************************************************************************************************************/
ParameterGroup* Algorithm::GetParamGroupPtr(void) {
    return m_pParamGroup;
}


/**************************************************************************************************************************************************************
MPI Communication - ConfigureWorkerDirectory()

Sends the worker subdirectory stem from the primary to the secondary worker
**************************************************************************************************************************************************************/
void Algorithm::ConfigureWorkerDirectory(int workerRank) {
    
    // Send the worker file stem to the secondary worker
    MPI_Send(m_DirPrefix, strlen(m_DirPrefix) + 1, MPI_CHAR, workerRank, tag_directory, MPI_COMM_WORLD);
}

/**************************************************************************************************************************************************************
MPI Communication - ConfigureWorkerSolveCommand()

Sends the worker solve command from the primary to the secondary worker
**************************************************************************************************************************************************************/
void Algorithm::ConfigureWorkerSolveCommand(int workerRank) {

    // Send the worker file stem to the secondary worker
    MPI_Send(m_ExecCmd, strlen(m_ExecCmd) + 1, MPI_CHAR, workerRank, tag_solve, MPI_COMM_WORLD);
    
}

/**************************************************************************************************************************************************************
MPI Communication - TerminateWorkers()

Sends the worker solve command from the primary to the secondary worker
**************************************************************************************************************************************************************/
void  Algorithm::TerminateWorkers() {
    // todo: build this function

    // Get the number of workers in the MPI space
    int numberOfMpiProcesses;
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfMpiProcesses);

    for (int entryWorker = 1; entryWorker < numberOfMpiProcesses; entryWorker++) {
        SendWorkerContinue(entryWorker, false);
    }
}

/**************************************************************************************************************************************************************
MPI Communication - ConfigureWorkerArchiveCommand()

Sends the worker archive command from the primary to the secondary worker
**************************************************************************************************************************************************************/
void Algorithm::ConfigureWorkerArchiveCommand(int workerRank, bool bMPI) {

    // Define the variable that will be sent to the workers
    int sendValue;

    // Send information to preserve the individual model runs
    if (!m_bPreserveModelOutput) {
        // Model will not be preserved
        sendValue = 0;
        
        if (bMPI) {
            // Send configuration to MPI worker
            MPI_Send(&sendValue, 1, MPI_INT, workerRank, tag_archive, MPI_COMM_WORLD);

        } else {
            // Set configuration to the local primary worker
            m_primaryWorker.SetWorkerPreserveArchiveCommand(sendValue, m_PreserveOutputCmd.string());
        }
        
    } else if (m_bPreserveModelOutput && m_PreserveOutputCmd.string().length() == 0) {
        // Model will be preserved with the the standard full archive method
        sendValue = 1;

        if (bMPI) {
            // Send configuration to MPI worker
            MPI_Send(&sendValue, 1, MPI_INT, workerRank, tag_archive, MPI_COMM_WORLD);

        } else {
            // Set configuration to the local primary worker
            m_primaryWorker.SetWorkerPreserveArchiveCommand(sendValue, m_PreserveOutputCmd.string());
        } 

    } else if (m_bPreserveModelOutput && m_PreserveOutputCmd.string().length() > 0) {
        // Model will be preserved with a custom archive command
        sendValue = 2;

        // Send the worker file stem to the secondary worker
        if (bMPI) {
            // Send configuration to MPI worker
            MPI_Send(&sendValue, 1, MPI_INT, workerRank, tag_archive, MPI_COMM_WORLD);
            MPI_Send(&m_PreserveOutputCmd.string()[0], m_PreserveOutputCmd.string().length() + 1, MPI_CHAR, workerRank, tag_archive, MPI_COMM_WORLD);

        } else {
            // Set configuration to the local primary worker
            m_primaryWorker.SetWorkerPreserveArchiveCommand(sendValue, m_PreserveOutputCmd.string());
        }
    }

    // Send information to preserve the best of the model run
    if (!m_bPreserveModelBest) {
        // Model will not be preserved
        sendValue = 0;
        
        if (bMPI) {
            // Send configuration to MPI worker
            MPI_Send(&sendValue, 1, MPI_INT, workerRank, tag_archive, MPI_COMM_WORLD);

        } else {
            // Set configuration to the local primary worker
            m_primaryWorker.SetWorkerPreserveBestCommand(sendValue, m_PreserveBestCmd.string());
        }
        
    } else if (m_bPreserveModelBest && m_PreserveBestCmd.string().length() == 0) {
        // Model will be preserved with the the standard full archive method
        sendValue = 1;

        if (bMPI) {
            // Send configuration to MPI worker
            MPI_Send(&sendValue, 1, MPI_INT, workerRank, tag_archive, MPI_COMM_WORLD);

        } else {
            // Set configuration to the local primary worker
            m_primaryWorker.SetWorkerPreserveBestCommand(sendValue, m_PreserveBestCmd.string());
        }
        
    } else if (m_bPreserveModelBest && m_PreserveBestCmd.string().length() > 0) {
        // Model will be preserved with a custom archive command
        sendValue = 2;

        if (bMPI) {
            // Send configuration to MPI worker
            MPI_Send(&sendValue, 1, MPI_INT, workerRank, tag_archive, MPI_COMM_WORLD);

            // Send the worker file stem to the secondary worker
            MPI_Send(&m_PreserveBestCmd.string()[0], m_PreserveBestCmd.string().length() + 1, MPI_CHAR, workerRank, tag_archive, MPI_COMM_WORLD);

        } else {
            // Set configuration to the local primary worker
            m_primaryWorker.SetWorkerPreserveBestCommand(sendValue, m_PreserveBestCmd.string());
        }
    }
}

/**************************************************************************************************************************************************************
MPI Communication - ConfigureWorkerExtraFiles()

Sends the extra filenames for the model from the primary to the secondary worker
**************************************************************************************************************************************************************/
void  Algorithm::ConfigureWorkerExtraFiles(int workerRank, bool bMPI) {

    if (bMPI) {
        // Send the number of extra files to the worker
        // Create a holder for the number of files
        int numberOfFiles = 0;

        // Loop and count the files
        FileList* pCur;
        for (pCur = m_pFileCleanupList; pCur != NULL; pCur = pCur->GetNext()) {
            numberOfFiles++;
        }

        // Send the number
        MPI_Send(&numberOfFiles, 1, MPI_INT, workerRank, tag_textLength, MPI_COMM_WORLD);

        // Send each file to the worker
        for (pCur = m_pFileCleanupList; pCur != NULL; pCur = pCur->GetNext()) {
            // Get the filename from the file list
            std::string fileName = pCur->GetName();

            // Send it to the secondary worker
            MPI_Send(&fileName[0], fileName.length() + 1, MPI_CHAR, workerRank, tag_textFile, MPI_COMM_WORLD);
        }

    } else {
        // Setup and set values in the local primary worker
        // Define the list
        std::vector<std::string> files;

        // Loop and count the files
        FileList* pCur;
        for (pCur = m_pFileCleanupList; pCur != NULL; pCur = pCur->GetNext()) {
            // Get the filename from the file list
            std::string fileName = pCur->GetName();

            // Append to the list
            files.push_back(fileName);
        }

        // Set into the primary worker
        m_primaryWorker.SetWorkerExtraFiles(files);
    }
}


/**************************************************************************************************************************************************************
MPI Communication - ConfigureWorkerFilePairs()

Sends the file pairs from the primary to the secondary worker
**************************************************************************************************************************************************************/
void Algorithm::ConfigureWorkerFilePairs(int workerRank, bool bMPI) {

    if (bMPI) {
        // Send the number of extra files to the worker
        // Create a holder for the number of files
        int numberOfFiles = fileListPairs.size();

        // Send the number
        MPI_Send(&numberOfFiles, 1, MPI_INT, workerRank, tag_fileLength, MPI_COMM_WORLD);

        // Send each file to the worker
        for (int entryPair = 0; entryPair < numberOfFiles; entryPair++) {
            // Get the filename from the file list
            std::string templateName = fileListPairs[entryPair][0];
            std::string destinationName = fileListPairs[entryPair][1];

            // Send it to the secondary worker
            MPI_Send(&templateName[0], templateName.length() + 1, MPI_CHAR, workerRank, tag_filePairs, MPI_COMM_WORLD);
            MPI_Send(&destinationName[0], destinationName.length() + 1, MPI_CHAR, workerRank, tag_filePairs, MPI_COMM_WORLD);
        }

    } else {
        // Setup and set values in the local primary worker
        m_primaryWorker.SetWorkerFilePairs(fileListPairs);
    }
};

/**************************************************************************************************************************************************************
MPI Communication - ConfigureWorkerObservations()

Sends the observations from from the primary to the secondary worker
**************************************************************************************************************************************************************/
void Algorithm::ConfigureWorkerObservations(int workerRank, bool bMPI) {

    if (bMPI) {
        // Send the number of extra files to the worker
        // Get the count variables
        int numberOfObservations = m_pObsGroup->GetNumObs();
        int numberOfGroups = m_pObsGroup->GetNumGroups();

        // Send the count variables
        MPI_Send(&numberOfObservations, 1, MPI_INT, workerRank, tag_obsLengthNum, MPI_COMM_WORLD);
        MPI_Send(&numberOfGroups, 1, MPI_INT, workerRank, tag_obsLengthGroup, MPI_COMM_WORLD);

        // Loop over the observations
        for (int entryObservation = 0; entryObservation < numberOfObservations; entryObservation++) {
            // Get the observations from the list
            Observation* obs = m_pObsGroup->GetObsPtr(entryObservation);

            // Get the values from the observation
            std::string obsName = obs->GetName();
            double obsValue = obs->GetMeasuredValueUntransformed();
            double obsWeight = obs->GetWeightFactor();
            std::string obsFile = obs->GetFileName();
            std::string obsKeyword = obs->GetKeyword();
            int obsLine = obs->GetLine();
            int obsColumn = obs->GetColumn();
            char obsToken = obs->GetToken();

            bool obsAug = obs->IsAugmented();
            int obsAugValue = 0;
            if (obsAug) {
                obsAugValue = 1;
            }
            std::string obsGroup = obs->GetGroup();

            // Send each observation component to the secondary worker
            MPI_Send(&obsName[0], obsName.length() + 1, MPI_CHAR, workerRank, tag_obsName, MPI_COMM_WORLD);
            MPI_Send(&obsValue, 1, MPI_DOUBLE, workerRank, tag_obsValue, MPI_COMM_WORLD);
            MPI_Send(&obsWeight, 1, MPI_DOUBLE, workerRank, tag_obsWeight, MPI_COMM_WORLD);
            MPI_Send(&obsFile[0], obsFile.length() + 1, MPI_CHAR, workerRank, tag_obsFile, MPI_COMM_WORLD);
            MPI_Send(&obsKeyword[0], obsKeyword.length() + 1, MPI_CHAR, workerRank, tag_obsKeyword, MPI_COMM_WORLD);
            MPI_Send(&obsLine, 1, MPI_INT, workerRank, tag_obsLine, MPI_COMM_WORLD);
            MPI_Send(&obsColumn, 1, MPI_INT, workerRank, tag_obsColumn, MPI_COMM_WORLD);
            MPI_Send(&obsToken, 1, MPI_CHAR, workerRank, tag_obsToken, MPI_COMM_WORLD);
            MPI_Send(&obsAugValue, 1, MPI_INT, workerRank, tag_obsAugmented, MPI_COMM_WORLD);
            MPI_Send(&obsGroup[0], obsGroup.length() + 1, MPI_CHAR, workerRank, tag_obsGroup, MPI_COMM_WORLD);

        }

    } else {
        // Set the observation group into the workder
        m_primaryWorker.SetWorkerObservations(m_pObsGroup);
    }
    
}

/**************************************************************************************************************************************************************
MPI Communication - ConfigureWorkerParameterGroups()

Sends the parameters from from the primary to the secondary worker
**************************************************************************************************************************************************************/
void Algorithm::ConfigureWorkerParameterGroups(int workerRank, bool bMPI) {

    if (bMPI) {
        // Send the total number of parameters
        int numberOfTotalParameters = m_pParamGroup->GetNumParams() + m_pParamGroup->GetNumTiedParams();
        MPI_Send(&numberOfTotalParameters, 1, MPI_INT, workerRank, tag_paramTotalNum, MPI_COMM_WORLD);

        // Real parameters
        int numberOfRealParameters = m_pParamGroup->GetNumRealParams();
        MPI_Send(&numberOfRealParameters, 1, MPI_INT, workerRank, tag_paramTotalReal, MPI_COMM_WORLD);

        for (int entryParameter = 0; entryParameter < numberOfRealParameters; entryParameter++) {
            // Get the parameter from the group, casing to the correct type
            RealParam* param = (RealParam*)m_pParamGroup->GetParamPtr(entryParameter);

            // Get the information from the parameter
            std::string name = param->GetName();
            double sampleValue = param->GetEstimatedValueTransformed();
            sampleValue = param->ConvertOutVal(sampleValue);
            std::string fixFmt = param->GetFixFmt();

            // Send the values to the worker
            MPI_Send(&name[0], name.length() + 1, MPI_CHAR, workerRank, tag_paramRealName, MPI_COMM_WORLD);
            MPI_Send(&sampleValue, 1, MPI_DOUBLE, workerRank, tag_paramRealInit, MPI_COMM_WORLD);
            MPI_Send(&fixFmt[0], fixFmt.length() + 1, MPI_CHAR, workerRank, tag_paramRealFmt, MPI_COMM_WORLD);
        }

        // Integer parameters
        // Value is offset by the number of real parameters to start indexing the pList correctly
        int numberOfIntegerParameters = m_pParamGroup->GetNumIntParams();
        MPI_Send(&numberOfIntegerParameters, 1, MPI_INT, workerRank, tag_paramTotalInt, MPI_COMM_WORLD);

        for (int entryParameter = numberOfRealParameters; entryParameter < numberOfRealParameters + numberOfIntegerParameters; entryParameter++) {
            // Get the parameter from the group, casing to the correct type
            IntParam* param = (IntParam*)m_pParamGroup->GetParamPtr(entryParameter);

            // Get the information from the parameter
            std::string name = param->GetName();
            int sampleValue = param->GetEstimatedValueTransformed();
            sampleValue = param->ConvertOutVal(sampleValue);

            // Send the values to the worker
            MPI_Send(&name[0], name.length() + 1, MPI_CHAR, workerRank, tag_paramInitName, MPI_COMM_WORLD);
            MPI_Send(&sampleValue, 1, MPI_INT, workerRank, tag_paramIntInit, MPI_COMM_WORLD);
        }

        // Tied parameters 
        // These will appear as more real parameters on the secondary worker
        int numberOfTiedParameters = m_pParamGroup->GetNumTiedParams();
        MPI_Send(&numberOfTiedParameters, 1, MPI_INT, workerRank, tag_paramTotalReal, MPI_COMM_WORLD);

        for (int entryParameter = 0; entryParameter < numberOfTiedParameters; entryParameter++) {
            // Get the parameter from the group
            TiedParamABC* param = m_pParamGroup->GetTiedParamPtr(entryParameter);

            // Get the information from the parameter
            std::string name = param->GetName();
            double sampleValue = param->GetEstimatedValueTransformed();
            std::string fixFmt = param->GetFixFmt();

            // Send the values to the worker
            MPI_Send(&name[0], name.length() + 1, MPI_CHAR, workerRank, tag_paramRealName, MPI_COMM_WORLD);
            MPI_Send(&sampleValue, 1, MPI_DOUBLE, workerRank, tag_paramRealInit, MPI_COMM_WORLD);
            MPI_Send(&fixFmt[0], fixFmt.length() + 1, MPI_CHAR, workerRank, tag_paramRealFmt, MPI_COMM_WORLD);
        }

        // Geom parameters
        // These will appear as more real parameters on the secondary worker
        /*int numberOfGeomParameters = m_pParamGroup->m_NumGeom;
        MPI_Send(&numberOfGeomParameters, 1, MPI_INT, workerRank, tag_paramTotalReal, MPI_COMM_WORLD);

        for (int entryParameter = 0; entryParameter < numberOfGeomParameters; entryParameter++) {
            // Get the parameter from the group
            GeomParamABC* param = m_pParamGroup->GetGeomParamPtr(entryParameter);

            // Get the information from the parameter
            std::string name = param->GetName();
            double sampleValue = param->GetEstimatedValueTransformed();
            std::string fixFmt = param->GetFixFmt();

            // Send the values to the worker
            MPI_Send(&name[0], name.length() + 1, MPI_CHAR, workerRank, tag_paramRealName, MPI_COMM_WORLD);
            MPI_Send(&sampleValue, 1, MPI_DOUBLE, workerRank, tag_paramRealInit, MPI_COMM_WORLD);
            MPI_Send(&fixFmt[0], fixFmt.length() + 1, MPI_CHAR, workerRank, tag_paramRealFmt, MPI_COMM_WORLD);
        }*/
    } else {
        // Set values into the local primary worker
        m_primaryWorker.SetWorkerParameters(m_pParamGroup);
    }
    
}

/**************************************************************************************************************************************************************
MPI Communication - ConfigureWorkers()

Sends all necessary information from the primary to the secondary workers
**************************************************************************************************************************************************************/
void Algorithm::ConfigureWorkers() {

    // Get the number of workers in the MPI space
    int numberOfMpiProcesses;
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfMpiProcesses);

    // Loop on the workers and send the information
    for (int entryWorker = 1; entryWorker < numberOfMpiProcesses; entryWorker++) {
        
        // Configure the worker directory
        ConfigureWorkerDirectory(entryWorker);

        // Configure the worker solve command
        ConfigureWorkerSolveCommand(entryWorker);

        // Configure the worker archive command
        ConfigureWorkerArchiveCommand(entryWorker, true);

        // Get the extra files from the directory
        ConfigureWorkerExtraFiles(entryWorker, true);

        // Get the file pairs for the templates and working files
        ConfigureWorkerFilePairs(entryWorker, true);

        // Get the observations for comparison
        ConfigureWorkerObservations(entryWorker, true);

        // Get the parameters for the analysis
        ConfigureWorkerParameterGroups(entryWorker, true);

    }

    // Setup the primary worker if using
    if (m_bSolveOnPrimary) {
        // Construct a local model worker
        m_primaryWorker = ModelWorker(false);

        // Configure the worker directory
        m_primaryWorker.SetWorkerDirectory(m_DirPrefix);

        // Configure the worker solve command
        m_primaryWorker.SetWorkerSolveCommand(m_ExecCmd);

        // Configure the worker archive command
        ConfigureWorkerArchiveCommand(0, false);

        // Get the extra files from the directory
        ConfigureWorkerExtraFiles(0, false);

        // Get the file pairs for the templates and working files
        ConfigureWorkerFilePairs(0, false);

        // Get the observations for comparison
        ConfigureWorkerObservations(0, false);

        // Get the parameters for the analysis
        ConfigureWorkerParameterGroups(0, false);
        
        // Configure the worker directory
        m_primaryWorker.SetupWorker();
    }
    
}

/**************************************************************************************************************************************************************
MPI Communication - SendWorkerParameters()

Sends parameters from the primary to the secondary workers
**************************************************************************************************************************************************************/
void Algorithm::SendWorkerParameters(int workerRank, int alternativeIndex, std::vector<double> parametersRegular) {

    // Create a vector to hold the parameters
    std::vector<double> parametersTemp;

    // Put the index in as the first entry
    parametersTemp.push_back((double)alternativeIndex);

    // Convert the parameters from transformation to normal space
    for (int entryParameter = 0; entryParameter < m_pParamGroup->m_NumParams; entryParameter++) {
        // Extract the pointer
        ParameterABC *param = m_pParamGroup->GetParamPtr(entryParameter);

        // Update the parameter with the converted value
        param->SetEstimatedValueTransformed(parametersRegular[entryParameter]);

        // Get the converted value
        parametersTemp.push_back(param->ConvertOutVal(parametersRegular[entryParameter]));
    }

    // Update the tied parameters
    for (int entryParameter = 0; entryParameter < m_pParamGroup->m_NumTied; entryParameter++) {
        // Extract the pointer
        TiedParamABC* param = m_pParamGroup->GetTiedParamPtr(entryParameter);

        // Get the updated value
        parametersTemp.push_back(param->GetEstimatedValueTransformed());
    }
 
    // Get the total number of parameters
    int numberOfParamters = m_pParamGroup->m_NumParams + m_pParamGroup->m_NumTied;   

    // Send the array to the secondary worker
    //MPI_Send(&parametersTemp[0], numberOfParamters + 1, MPI_DOUBLE, workerRank, tag_data, MPI_COMM_WORLD);
    MPI_Ssend(parametersTemp.data(), numberOfParamters + 1, MPI_DOUBLE, workerRank, tag_data, MPI_COMM_WORLD);

}

/**************************************************************************************************************************************************************
MPI Communication - SendWorkerContinue()

Sends parameters from the primary to the secondary workers
**************************************************************************************************************************************************************/
void Algorithm::SendWorkerContinue(int workerRank, bool workerContinue) {

    int transferValue = 1;
    if (workerContinue) {
        MPI_Ssend(&transferValue, 1, MPI_INT, workerRank, tag_continue, MPI_COMM_WORLD);

    } else {
        transferValue = 0;
        MPI_Ssend(&transferValue, 1, MPI_INT, workerRank, tag_continue, MPI_COMM_WORLD);
    }

}

/**************************************************************************************************************************************************************
MPI Communication - SendWorkerPreservebest()

Sends parameters from the primary to the secondary workers
**************************************************************************************************************************************************************/
void Algorithm::SendWorkerPreserveBest(int workerRank, bool preserveModel) {

    int transferValue = 1;
    if (preserveModel) {
        MPI_Ssend(&transferValue, 1, MPI_INT, workerRank, tag_preserveBest, MPI_COMM_WORLD);

    } else {
        transferValue = 0;
        MPI_Ssend(&transferValue, 1, MPI_INT, workerRank, tag_preserveBest, MPI_COMM_WORLD);
    }
}


/**************************************************************************************************************************************************************
MPI Communication - ReceiveWorkerPreserveBest()

Sends from worker that the model preserve is complete
**************************************************************************************************************************************************************/
void Algorithm::ReceiveWorkerPreserveBest() {

    int transferValue;
    MPI_Recv(&transferValue, 1, MPI_INT, MPI_ANY_SOURCE, tag_preserveBest, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}


/**************************************************************************************************************************************************************
MPI Communication - ManagePreserveBest()

Manages the process for preserving the best model
**************************************************************************************************************************************************************/
void Algorithm::ManagePreserveBest(double &solutionObjective, double alternativeObjective, MPI_Status mpiStatus) {


    if (m_bPreserveModelBest) {
        // Determine if the worker should archive the model
        if (alternativeObjective < solutionObjective) {
            // Solution value is better the previous best solution.
            // Replace the objective value
            solutionObjective = alternativeObjective;

            // Tell the worker to preserve the model
            SendWorkerPreserveBest(mpiStatus.MPI_SOURCE, true);

            // Wait for confirmation
            ReceiveWorkerPreserveBest();

        }
        else {
            // Tell the worker not to preserve the output
            SendWorkerPreserveBest(mpiStatus.MPI_SOURCE, false);
        }
    }
}


/**************************************************************************************************************************************************************
MPI Communication - ManageSingleObjectiveIterations()

Manages the solution of a set of parameter alternatives when using a single objective function
**************************************************************************************************************************************************************/
void Algorithm::ManageSingleObjectiveIterations(std::vector<std::vector<double>> parameters, int startingIndex, int numberOfParameters, 
                                                std::vector<double>& returnArray) {

    // Get the number of workers in the MPI space
    int numberOfMpiProcesses;
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfMpiProcesses);

    // Create a copy of the current best objective function
    double bestObjectiveIteration = m_BestObjective;

    // If caching is enabled, check the solutions against the list of values to be created
    std::vector<int> indicesToSolve;
    for (int entryAlternative = startingIndex; entryAlternative < parameters.size(); entryAlternative++) {
        if (m_bCaching) {
            // Caching is enabled. Check that the member is not in the cache
            std::vector<std::vector<double>> ::iterator itr = std::find(m_CacheMembers.begin(), m_CacheMembers.end(), parameters[entryAlternative]);

            if (itr == m_CacheMembers.end()) {
                // Alternative is not present in the cache. Add it to be solved.
                indicesToSolve.push_back(entryAlternative);
            }
            else {
                // Log the cache hit
                m_NumCacheHits++;

                // Add the objective to the result vector
                returnArray[entryAlternative] = m_CacheObjectives[std::distance(m_CacheMembers.begin(), itr)];
            }

        } else {
            // Caching is not enabled. Add it to the solve list
            indicesToSolve.push_back(entryAlternative);
        }
    }

    // Setup the request space
    MPI_Message mpiMessage;
    MPI_Status mpiStatus;

    int requestFlag = 1;
    int sendFlag = 1;

    // Set the counter variables         
    int sendCounter = 0;
    int receiveCounter = 0;

    // Loop and send work to the secondary workers as they become available
    while (sendCounter < indicesToSolve.size()){  
        // Get the rank number of a secondary worker waiting for work
        int workerRank, workerValue;

        // Receive until no more receives are available
        while (requestFlag) {
            // Attempt to receive from the workers
            MPI_Improbe(MPI_ANY_SOURCE, tag_data, MPI_COMM_WORLD, &requestFlag, &mpiMessage, &mpiStatus);
            
            // Handle any available receives
            if (requestFlag) {
                try {
                    // Get the alternative index and objective function from the secondary worker
                    double data[2];
                    MPI_Mrecv(&data, 2, MPI_DOUBLE, &mpiMessage, &mpiStatus);

                    // Set the objective value into the objective array
                    returnArray[(int)data[0]] = data[1];

                    // Store the alternative into the cache, if enabled
                    if (m_bCaching) {
                        m_CacheMembers.push_back(parameters[(int)data[0]]);
                        m_CacheObjectives.push_back(data[1]);
                    }

                    // Check if the alternative should be preserved
                    ManagePreserveBest(bestObjectiveIteration, data[1], mpiStatus);

                } catch (...) {
                    // Read from the secondary worker failed. Fail by keeping the infinite value in the array
                }

                // Increment the receive counter
                receiveCounter++;

                // Increment the solve counter
                m_NumSolves++;
            }
        }

        // Send until no more workers are available
        while (sendFlag) {
            // Attempt to send to the workers
            MPI_Improbe(MPI_ANY_SOURCE, tag_continue, MPI_COMM_WORLD, &sendFlag, &mpiMessage, &mpiStatus);

            // Handle any available sends
            if (sendFlag) {
                // Worker is available to accept work. Go through the steps to send work to it.
                // Receive the continue flag from the worker
                MPI_Mrecv(&workerValue, 1, MPI_INT, &mpiMessage, &mpiStatus);

                // Extract the worker rank from the status request
                workerRank = mpiStatus.MPI_SOURCE;

                // Transmit a continue flag to the idle secondary worker
                SendWorkerContinue(workerRank, true);

                // Include the tied parameters
                std::vector<double> parametersWithTied = AddTiedParametersToAlternative(parameters[indicesToSolve[sendCounter]]);

                // Transmit the parameteres to the worker to solve
                SendWorkerParameters(workerRank, indicesToSolve[sendCounter], parametersWithTied);

                // Increment the send counter
                sendCounter++;

                // Stop sending if all values have been completed
                if (sendCounter == indicesToSolve.size()) {
                    break;
                }
            }
        }

        // Solve on the primary worker if enabled
        if (m_bSolveOnPrimary && sendCounter < indicesToSolve.size()) {
            // Include the tied parameters
            std::vector<double> parametersWithTied = AddTiedParametersToAlternative(parameters[indicesToSolve[sendCounter]]);

            // Set the parameters on the primary worker
            m_primaryWorker.SetStandardParameters(parametersWithTied);

            // Attempt to solve the model
            try {
                // Solve the model
                double objective =  m_primaryWorker.StdExecute(0);
                returnArray[indicesToSolve[sendCounter]] = objective;

                // Store the alternative into the cache, if enabled
                if (m_bCaching) {
                    m_CacheMembers.push_back(parameters[indicesToSolve[sendCounter]]);
                    m_CacheObjectives.push_back(objective);
                }

                // Check if the alternative should be preserved
                if (m_bPreserveModelBest && objective < bestObjectiveIteration) {
                    // Objective has improved. Swap with the current value
                    bestObjectiveIteration = objective;

                    // Preserve the best model
                    m_primaryWorker.PreserveBestModel(true);
                }
              
                if (m_bPreserveModelOutput){
                    // Preserve the model
                    m_primaryWorker.PreserveArchiveModel();
                }
                
            }
            catch (...){
                // Solution failed. Fail by keeping the infinite value in the array
            }

            // Increment the counters
            sendCounter++;
            receiveCounter++;

        }
        
        // Reset the flags for the next loop
        requestFlag = 1;
        sendFlag = 1;

    }

    // All alternatives have been sent. Receive all of the solutions
    while (receiveCounter < indicesToSolve.size()) {
        // Determine if we can receive any values
        MPI_Improbe(MPI_ANY_SOURCE, tag_data, MPI_COMM_WORLD, &requestFlag, &mpiMessage, &mpiStatus);
        if (requestFlag) {
            // Get the alternative index and objective function from the secondary worker
            double data[2];
            MPI_Mrecv(&data, 2, MPI_DOUBLE, &mpiMessage, &mpiStatus);

            // Set the objective value into the objective array
            returnArray[(int)data[0]] = data[1];

            // Store the alternative into the cache, if enabled
            if (m_bCaching) {
                m_CacheMembers.push_back(parameters[(int)data[0]]);
            }

            // Check if the alternative should be preserved
            ManagePreserveBest(bestObjectiveIteration, data[1], mpiStatus);

            // Increment the receive counter
            receiveCounter++;

            // Increment the solve counter
            m_NumSolves++;
        }
    }    
}

std::vector<double> Algorithm::AddTiedParametersToAlternative(std::vector<double> parameters) {

    // Need to introduce tied parameters here
// Set the alternative values into the parameter object
    for (int entryParameter = 0; entryParameter < m_pParamGroup->GetNumParams(); entryParameter++) {
        // Get the parameter
        ParameterABC* temp = m_pParamGroup->GetParamPtr(entryParameter);

        // Set the value into the parameter
        temp->SetEstimatedValueTransformed(parameters[entryParameter]);
    }

    // Get the tied parameter values
    std::vector<double> parametersWithTied = parameters;
    for (int entryParameterTied = 0; entryParameterTied < m_pParamGroup->GetNumTiedParams(); entryParameterTied++) {
        // Get the pointer to the tied parameter
        TiedParamABC* temp = m_pParamGroup->GetTiedParamPtr(entryParameterTied);

        // Get the value and pust back to the 
        parametersWithTied.push_back(temp->GetEstimatedValueTransformed());
    }


    return parametersWithTied;
}



/**************************************************************************************************************************************************************
GetBestSingleObjective()

Retrieves the chromosome value and index that has the best fitness value.
**************************************************************************************************************************************************************/
void Algorithm::GetBestSingleObjective(std::vector<double> objectives, double& bestObjective, int& bestIndex) {

    bestIndex = std::min_element(objectives.begin(), objectives.end()) - objectives.begin();
    bestObjective = *std::min_element(objectives.begin(), objectives.end());

}

/**************************************************************************************************************************************************************
CheckGlobalSensitivity()

Checks that each observation is sensitive to at least one parameter over the range of possible parameter values. If an observation is not sensitive to any
parameters, a warning will be reported and the observation will be ignored in the calibration.

Also checks the sensitivity of each parameter. If a given parameter does not affect the objective function over the entire parameter range, a warning will
be reported and the parameter will be ignored in the calibration.
**************************************************************************************************************************************************************/
/*void Algorithm::CheckGlobalSensitivity(void)
{
    if ((m_pObsGroup == NULL) || (m_pParamGroup == NULL)) return;
    if (m_bCheckGlobalSens == false) return;

    int i, j, nobs, nprm;
    double upr, lwr, Fupr, Flwr;
    double* pInit, * obs_sum, * prm_sum;
    double* ObsUpr, * ObsLwr;
    char tmp1[DEF_STR_SZ];
    UnchangeableString* obs_names, * prm_names;

    nobs = m_pObsGroup->GetNumObs();
    obs_sum = new double[nobs];
    ObsUpr = new double[nobs];
    ObsLwr = new double[nobs];
    obs_names = new UnchangeableString[nobs];

    nprm = m_pParamGroup->GetNumParams();
    pInit = new double[nprm];
    prm_sum = new double[nprm];
    prm_names = new UnchangeableString[nprm];

    //save intitial parameters
    m_pParamGroup->ReadParams(pInit);

    for (i = 0; i < nobs; i++)
    {
        obs_sum[i] = 0.00;
        obs_names[i] = m_pObsGroup->GetObsPtr(i)->GetName();
    }

    for (j = 0; j < nprm; j++)
    {
        prm_names[j] = m_pParamGroup->GetParamPtr(j)->GetName();
        upr = m_pParamGroup->GetParamPtr(j)->GetUprBnd();
        m_pParamGroup->GetParamPtr(j)->SetEstVal(upr);
        Fupr = Execute();

        //store computed observation values
        for (i = 0; i < nobs; i++)
        {
            ObsUpr[i] = m_pObsGroup->GetObsPtr(i)->GetComputedVal(true, true);
        }

        lwr = m_pParamGroup->GetParamPtr(j)->GetLwrBnd();
        m_pParamGroup->GetParamPtr(j)->SetEstVal(lwr);
        Flwr = Execute();

        //store computed observation values
        for (i = 0; i < nobs; i++)
        {
            ObsLwr[i] = m_pObsGroup->GetObsPtr(i)->GetComputedVal(true, true);
        }

        //store change in objective function
        prm_sum[j] = fabs(Fupr - Flwr);

        //accumulate changes in observation values
        for (i = 0; i < nobs; i++)
        {
            obs_sum[i] += fabs(ObsUpr[i] - ObsLwr[i]);
        }

        //restore parameter values
        m_pParamGroup->WriteParams(pInit);
    }

    //perform parameter sensitivity checks
    for (j = 0; j < nprm; j++)
    {
        if (prm_sum[j] <= NEARLY_ZERO)
        {
            sprintf(tmp1, "%s appears to be insensitive and has been set to a constant value", prm_names[j]);
            LogError(ERR_INS_PARM, tmp1);
            m_pParamGroup->ExcludeParam(prm_names[j]);
        }
    }

    //perform observation sensitivity checks: phase 1, log the warning messages
    for (i = 0; i < nobs; i++)
    {
        if (obs_sum[i] <= NEARLY_ZERO)
        {
            sprintf(tmp1, "%s appears to be insensitive and has been excluded from the calibration", obs_names[i]);
            LogError(ERR_INS_OBS, tmp1);
            m_pObsGroup->ExcludeObs(obs_names[i]);
        }
    }

    delete[] pInit;
    delete[] obs_sum;
    delete[] prm_sum;
    delete[] ObsUpr;
    delete[] ObsLwr;
    delete[] prm_names;
    delete[] obs_names;
} /* end CheckGlobalSensitivity() */

