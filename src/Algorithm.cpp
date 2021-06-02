/******************************************************************************

The Model class encapsulates the interaction of the Ostrich optimization tools
with the externally executed modeling program. The class divides
model components into three groups: the parameter group, the observation group
and the objective function group. In addition to being able to execute the
model, the Model class provides Ostrich algorithms with access to
these groups.

******************************************************************************/
// TODO: updaate doc string
#include <mpi.h>
#include <math.h>
#include <string>
#include <iostream>

#include "Algorithm.h"
#include "ObservationGroup.h"
#include "Observation.h"
#include "ResponseVarGroup.h"
#include "SurrogateParameterGroup.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "TiedParamABC.h"
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


/******************************************************************************
 default CTOR
******************************************************************************/
Algorithm::Algorithm(void) {
    FilePair* pFilePair;
    char* line;
    char tmp1[DEF_STR_SZ];
    char tmp2[DEF_STR_SZ];
    char tmp3[DEF_STR_SZ];
    FILE* pInFile;
    int id, i, j;
    bool quoteWrap; //if true, then executable must be wrapped in quotes


    //RegisterModelPtr(this);
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
    Read in the model executable, modifying so that output is redirected to a file. Also, copy the executable file to the 
    directory of execution.
    --------------------------------------------------------------------------------------------------------------------------
    */
    rewind(pInFile);
    FindToken(pInFile, "ModelExecutable", inFileName);
    line = GetCurDataLine();


    /*
    --------------------------------------------------------------------------------------------------------------------------
    Read in executable, taking care to preserve full path, even in the presence of long and space-separated filenames.
    --------------------------------------------------------------------------------------------------------------------------
    */
    i = ExtractString(line, tmp2);
    i = ValidateExtraction(i, 1, 1, "Model()");
    i = ExtractFileName(&(line[i]), tmp1);

    //must wrap in quotes if there is whitespace in the execuable path
    //quoteWrap = false;
    //j = (int)strlen(tmp1);
    //for (i = 0; i < j; i++) { if (tmp1[i] == ' ') { quoteWrap = true; } }
    //if (quoteWrap == true) { tmp1[j++] = '"'; }
    //tmp1[j] = (char)NULL;
    //if (quoteWrap == true) {
    //    MyStrRev(tmp1);
    //    tmp1[j++] = '"';
    //    tmp1[j] = (char)NULL;
    //    MyStrRev(tmp1);
    //}

    /*
    --------------------------------------------------------------------------------------------------------------------------
    Extract OSTRICH executable file name from rest of path and save it. This is so that Ostrich can cleanup after itself when 
    copying files around.
    --------------------------------------------------------------------------------------------------------------------------
    */
    // TODO: make this process more robust
    /*#ifdef _WIN32
        m_pFileCleanupList = new FileList("ostrich.exe");
    #else
        m_pFileCleanupList = new FileList("Ostrich");
    #endif*/

    /*
    --------------------------------------------------------------------------------------------------------------------------
    Extract executable file name from rest of path and save it. This is so that Ostrich can cleanup after itself when copying 
    files around.
    --------------------------------------------------------------------------------------------------------------------------
    */
   /* if (m_InternalModel == false) {
        j = 0;
        for (i = (int)(strlen(tmp1)) - 1; i > 0; i--) {
            if ((tmp1[i] != '\\') && (tmp1[i] != '/')) {
                tmp2[j] = tmp1[i];
                j++;
            }
            else {
                if (quoteWrap == true) {
                    tmp2[j] = '"';
                    j++;
                }
                tmp2[j] = (char)NULL;
                MyStrRev(tmp2);
                m_pFileCleanupList->Insert(tmp2);
                break;
            }
        }
    }*/

    /*
    --------------------------------------------------------------------------------------------------------------------------
    Check to see if the model is internal. If it is, it will be a function call and end in a '()'. Setup the operations for an 
    internal model.
    --------------------------------------------------------------------------------------------------------------------------
    */
    if ((strlen(tmp1) > 1) && (tmp1[strlen(tmp1) - 2] == '(') && (tmp1[strlen(tmp1) - 1] == ')')) {
        m_InternalModel = true;
    }

    /*
    --------------------------------------------------------------------------------------------------------------------------
    Check to see if the model is external. Setup operations for an external model.
    --------------------------------------------------------------------------------------------------------------------------
    */
    if (m_InternalModel == false) {
        // Store the run command
        MyTrim(tmp1);
        int len = (int)strlen(tmp1) + 1;

        m_ExecCmd = new char[len];
        strcpy(m_ExecCmd, tmp1);
    }


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
        }/* end if() */

        /*
        ----------------------------------------------------------------------------------------------------------------------
        Read in file pairs, taking care to preserve full path,
        even in the presence of long and space-separated filenames.
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
    // TODO: Update with new filesystem syntx
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
    rewind(pInFile);
    if (CheckToken(pInFile, "DisklessModel", inFileName) == true) {
        line = GetCurDataLine();
        sscanf(line, "%s %s", tmp1, tmp2);
        MyStrLwr(tmp2);
        if (strncmp(tmp2, "yes", 3) == 0) { m_bDiskless = true; }
    }

    /*
    --------------------------------------------------------------------------------------------------------------------------
    Read in the name of the batch file to run whenever a new 'best' solution is found.
    --------------------------------------------------------------------------------------------------------------------------
    */
    rewind(pInFile);
    if (CheckToken(pInFile, "PreserveBestModel", inFileName) == true)
    {
        line = GetCurDataLine();
        //line format = 'PreserveBestModel   <var>'
        /*--------------------------------------------------------------------------------------------------------------------
        Read in executable, taking care to preserve full path, even in the presence of long and space-separated filenames.
        ----------------------------------------------------------------------------------------------------------------------
        */
        i = ExtractString(line, tmp2);
        i = ValidateExtraction(i, 1, 1, "Model()");
        i = ExtractFileName(&(line[i]), tmp1);

        //must wrap in quotes if there is whitespace in the execuable path
        quoteWrap = false;
        j = (int)strlen(tmp1);
        for (i = 0; i < j; i++) { if (tmp1[i] == ' ') { quoteWrap = true; } }
        if (quoteWrap == true) { tmp1[j++] = '"'; }
        tmp1[j] = (char)NULL;
        if (quoteWrap == true) {
            MyStrRev(tmp1);
            tmp1[j++] = '"';
            tmp1[j] = (char)NULL;
            MyStrRev(tmp1);
        }

        if (pDirName[0] != '.') {
            #ifdef _WIN32
                sprintf(tmp2, "copy %s %s", tmp1, pDirName);
            #else
                sprintf(tmp2, "cp %s %s", tmp1, pDirName);
            #endif
            system(tmp2);
        } 

        //make sure the executable exists
        strcpy(tmp2, tmp1);

        if (tmp2[0] == '"') {
            tmp2[0] = ' ';
            tmp2[strlen(tmp2) - 1] = ' ';
            MyTrim(tmp2);
        }

        if (MY_ACCESS(tmp2, 0) == -1){
            sprintf(tmp1, "File for saving best solution (|%s|) not found", tmp2);
            LogError(ERR_FILE_IO, tmp1);
            ExitProgram(1);
        }

        #ifdef _WIN32 //windows version
            strcat(tmp1, " > OstSaveOut.txt");
        #else //Linux (bash, dash, and csh)
            strcpy(tmp2, tmp1);
            // '>&' redircts both output and error
            strcat(tmp2, " > OstSaveOut.txt 2>&1");
            strcpy(tmp1, tmp2);
        #endif

        m_bSave = true;
        m_SaveCmd = new char[(int)strlen(tmp1) + 1];
        strcpy(m_SaveCmd, tmp1);
    }/* end if(PreserveBestModel) */

    /*
    --------------------------------------------------------------------------------------------------------------------------
    Check for alternate objective function (default is WSSE)
    --------------------------------------------------------------------------------------------------------------------------
    */
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
    rewind(pInFile);
    if (CheckToken(pInFile, "CheckSensitivities", inFileName) == true) {
        line = GetCurDataLine();
        sscanf(line, "%s %s", tmp1, tmp2);
        MyStrLwr(tmp2);
        if (strncmp(tmp2, "yes", 3) == 0) { m_bCheckGlobalSens = true; }
    }

    /*
    --------------------------------------------------------------------------------------------------------------------------
    Read in flag to use surrogate models
    --------------------------------------------------------------------------------------------------------------------------
    */
    rewind(pInFile);
    if (CheckToken(pInFile, "SurrogateApproach", inFileName) == true) {
        line = GetCurDataLine();
        sscanf(line, "%s %s", tmp1, tmp2);
        MyStrLwr(tmp2);
        if (strncmp(tmp2, "yes", 3) == 0) { m_bUseSurrogates = true; }
    }

    /*
    --------------------------------------------------------------------------------------------------------------------------
    Read in flag to use SuperMUSE
    --------------------------------------------------------------------------------------------------------------------------
    */
    rewind(pInFile);
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
    }

    /*
    --------------------------------------------------------------------------------------------------------------------------
    Read in flag to preserve model output files. This only applies if a model subdirectory is used.
    --------------------------------------------------------------------------------------------------------------------------
    */
    // TODO: rework this process
    rewind(pInFile);
    if (CheckToken(pInFile, "PreserveModelOutput", inFileName) == true) {
        line = GetCurDataLine();
        sscanf(line, "%s %s", tmp1, tmp2);
        MyStrLwr(tmp2);
        MyTrim(tmp2);
        if (strcmp(tmp2, "yes") == 0) {
            m_bPreserveModelOutput = true;
        }
        else if (strcmp(tmp2, "no") == 0) {
            m_bPreserveModelOutput = false;
        }
        else { //a valid PreserveModelOutput script?
            //line format = 'PreserveModelOutput  <var>'
            /*--------------------------------------------------------
            Read in executable, taking care to preserve full path,
            even in the presence of long and space-separated filenames.
            --------------------------------------------------------*/
            i = ExtractString(line, tmp2);
            i = ValidateExtraction(i, 1, 1, "Model()");
            i = ExtractFileName(&(line[i]), tmp1);

            //must wrap in quotes if there is whitespace in the execuable path
            quoteWrap = false;
            j = (int)strlen(tmp1);
            for (i = 0; i < j; i++) { if (tmp1[i] == ' ') { quoteWrap = true; } }
            if (quoteWrap == true) { tmp1[j++] = '"'; }
            tmp1[j] = (char)NULL;
            if (quoteWrap == true) {
                MyStrRev(tmp1);
                tmp1[j++] = '"';
                tmp1[j] = (char)NULL;
                MyStrRev(tmp1);
            }

            //stage to workdir if needed
            if (pDirName[0] != '.') {
                #ifdef _WIN32
                    sprintf(tmp2, "copy %s %s", tmp1, pDirName);
                #else
                    sprintf(tmp2, "cp %s %s", tmp1, pDirName);
                #endif
                system(tmp2);
            } 

            //make sure the executable exists
            strcpy(tmp2, tmp1);

            if (tmp2[0] == '"') {
                tmp2[0] = ' ';
                tmp2[strlen(tmp2) - 1] = ' ';
                MyTrim(tmp2);
            }
            if (MY_ACCESS(tmp2, 0) == -1) {
                sprintf(tmp1, "File for preserving model output (|%s|) not found", tmp2);
                LogError(ERR_FILE_IO, tmp1);
                ExitProgram(1);
            }

            m_bPreserveModelOutput = true;
            m_PreserveCmd = new char[(int)strlen(tmp1) + 1];
            strcpy(m_PreserveCmd, tmp1);
        } 
    }

    /*
    --------------------------------------------------------------------------------------------------------------------------
    Read in warm start flag. This only applies if we want to restart Ostrich from a previously aborted or interrupted run in 
    which the OstModel0.txt file has been preserved.
    --------------------------------------------------------------------------------------------------------------------------
    */
    rewind(pInFile);
    if (CheckToken(pInFile, "OstrichWarmStart", inFileName) == true) {
        line = GetCurDataLine();
        sscanf(line, "%s %s", tmp1, tmp2);
        MyStrLwr(tmp2);
        if (strncmp(tmp2, "yes", 3) == 0) {
            printf("Warm Start has been activated\n");
            printf("Ostrich will resume a pervious search.\n");
            m_bWarmStart = true;
            RestoreRandomSeed();
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
    if (CheckToken(pInFile, "OstrichCaching", inFileName) == true) {
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
    char tscp[DEF_STR_SZ];
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
    }

    /*
    --------------------------------------------------------------------------------------------------------------------------
    Read in BoxCox transformation, if any.
    --------------------------------------------------------------------------------------------------------------------------
    */
    char boxcox[DEF_STR_SZ];
    bool bBoxCox = false;
    double boxCoxValue = 1.00;
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

    if (bSMUSE) CleanSuperMUSE();

    /*
    --------------------------------------------------------------------------------------------------------------------------
    Read in database conversion information.
    --------------------------------------------------------------------------------------------------------------------------
    */
    DatabaseABC* access_DbaseList = new AccessConverter();
    DatabaseABC* netcdf_DbaseList = new NetCDFConverter();
    MEM_CHECK(access_DbaseList);
    MEM_CHECK(netcdf_DbaseList);
    if (access_DbaseList->ReadFromFile()) {
        NEW_PRINT("AccessConverter", 1);
        m_DbaseList = access_DbaseList;
        delete netcdf_DbaseList;
        netcdf_DbaseList = NULL;
    }
    else if (netcdf_DbaseList->ReadFromFile()) {
        NEW_PRINT("NetCDFConverter", 1);
        m_DbaseList = netcdf_DbaseList;
        delete access_DbaseList;
        access_DbaseList = NULL;
    }

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
    }/* end if() */

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
    Read in special parameters.
    --------------------------------------------------------------------------------------------------------------------------
    */
    m_pParamGroup->InitSpecialParams(inFileName);

    /*
    --------------------------------------------------------------------------------------------------------------------------
    Check template files against parameters, each parameter should appear in
    at least one template file or at least one database entry.
    --------------------------------------------------------------------------------------------------------------------------
    */
    //m_pParamGroup->CheckTemplateFiles(m_FileList);
    //m_pParamGroup->CheckDbaseFiles(m_DbaseList);

    /*
    --------------------------------------------------------------------------------------------------------------------------
    Check parameters for uniqueness, each parameter should be unique and
    should not be a substring of another parameter.
    --------------------------------------------------------------------------------------------------------------------------
    */
    m_pParamGroup->CheckMnemonics();

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

    if (m_bSolveOnPrimary) {
        // Setup a secondary worker object on the primary worker
        // TODO: create this

    }

    IncCtorCount();

} 

/*
------------------------------------------------------------------------------------------------------------------------------
Free up memory.
------------------------------------------------------------------------------------------------------------------------------
*/
void Algorithm::Destroy(void) {
    delete m_pObsGroup;
    delete m_pParamGroup;
    delete m_pParameterCorrection;
    delete m_pObjFunc;
    delete m_DbaseList;
    delete[] m_ExecCmd;
    delete[] m_SaveCmd;
    delete[] m_PreserveCmd;
    delete[] m_CurMultiObjF;
    m_bSave = false;
    delete m_pDecision;
    delete m_BestAlternative;

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


/*
------------------------------------------------------------------------------------------------------------------------------
PerformParameterCorrections()
------------------------------------------------------------------------------------------------------------------------------
*/
void Algorithm::PerformParameterCorrections(void) {
    if (m_pParameterCorrection != NULL) m_pParameterCorrection->Execute();
}/* end PerformParameterCorrections() */

/*
------------------------------------------------------------------------------------------------------------------------------
 SaveBest()
------------------------------------------------------------------------------------------------------------------------------
*/
void Algorithm::SaveBest(int id) {
    char saveDir[DEF_STR_SZ];
    if (m_bSave) {
        //cd to model subdirectory, if needed
        if (m_DirPrefix[0] != '.') {
            sprintf(saveDir, "%s%d", m_DirPrefix, id);
            MY_CHDIR(saveDir);
        }

        system(m_SaveCmd);

        if (m_DirPrefix[0] != '.'){
            //cd out of model subdirectory, if needed
            MY_CHDIR("..");
        }
    }
}/* end SaveBest() */


/*
------------------------------------------------------------------------------------------------------------------------------
GetObjFuncStr()
------------------------------------------------------------------------------------------------------------------------------
*/
UnchangeableString Algorithm::GetObjFuncStr(void) {
    return m_pObjFunc->GetObjFuncStr();
}

/*
------------------------------------------------------------------------------------------------------------------------------
GetObjFuncPtr()
   Returns a pointer to the objective function.
------------------------------------------------------------------------------------------------------------------------------
*/
ObjectiveFunction* Algorithm::GetObjFuncPtr(void) {
    return m_pObjFunc;
} /* end GetObjFuncPtr() */

/*
------------------------------------------------------------------------------------------------------------------------------
GetCounter()
   Returns the number of times the model has been executed
------------------------------------------------------------------------------------------------------------------------------
*/
int Algorithm::GetCounter(void) {
    return m_Counter;
}

/*
------------------------------------------------------------------------------------------------------------------------------
SetCounter()
------------------------------------------------------------------------------------------------------------------------------
*/
void Algorithm::SetCounter(int count) {
    m_Counter = count;
} 


/*
------------------------------------------------------------------------------------------------------------------------------
AddDatabase()
   Adds a database conversion to the list.
------------------------------------------------------------------------------------------------------------------------------
*/
void Algorithm::AddDatabase(DatabaseABC* pDbase) {
    if (m_DbaseList == NULL) { m_DbaseList = pDbase; }
    else { m_DbaseList->InsertDbase(pDbase); }
}

/*
------------------------------------------------------------------------------------------------------------------------------
GetObsGroupPtr()
   Returns the observation group pointer
------------------------------------------------------------------------------------------------------------------------------
*/
ObservationGroup* Algorithm::GetObsGroupPtr(void) {
    return m_pObsGroup;
}

/*
------------------------------------------------------------------------------------------------------------------------------
GetParamGroupPtr()
   Returns the parameter group pointer.
------------------------------------------------------------------------------------------------------------------------------
*/
ParameterGroup* Algorithm::GetParamGroupPtr(void) {
    return m_pParamGroup;
}



/*
------------------------------------------------------------------------------------------------------------------------------
MPI Communication - ConfigureWorkerDirectory()
   Transfers the worker subdirectory stem from the primary to the secondary worker
------------------------------------------------------------------------------------------------------------------------------
*/
void Algorithm::ConfigureWorkerDirectory(int workerRank) {
    
    // Send the worker file stem to the secondary worker
    MPI_Send(m_DirPrefix, strlen(m_DirPrefix) + 1, MPI_CHAR, workerRank, tag_directory, MPI_COMM_WORLD);
}

/*
------------------------------------------------------------------------------------------------------------------------------
MPI Communication - ConfigureWorkerSolveCommand()
   Transfers the worker solve command from the primary to the secondary worker
------------------------------------------------------------------------------------------------------------------------------
*/
void Algorithm::ConfigureWorkerSolveCommand(int workerRank) {

    // Send the worker file stem to the secondary worker
    MPI_Send(m_ExecCmd, strlen(m_ExecCmd) + 1, MPI_CHAR, workerRank, tag_solve, MPI_COMM_WORLD);
}

/*
------------------------------------------------------------------------------------------------------------------------------
MPI Communication - ConfigureWorkerExtraFiles()
   Transfers the extra filenames for the model from the primary to the secondary worker
------------------------------------------------------------------------------------------------------------------------------
*/
void  Algorithm::ConfigureWorkerExtraFiles(int workerRank) {

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
}


/*
------------------------------------------------------------------------------------------------------------------------------
MPI Communication - ConfigureWorkerFilePairs()
   Transfers the file pairs from the primary to the secondary worker
------------------------------------------------------------------------------------------------------------------------------
*/
void Algorithm::ConfigureWorkerFilePairs(int workerRank) {

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
};

/*
------------------------------------------------------------------------------------------------------------------------------
MPI Communication - ConfigureWorkerObservations()
   Transfers the observations from from the primary to the secondary worker
------------------------------------------------------------------------------------------------------------------------------
*/
void Algorithm::ConfigureWorkerObservations(int workerRank) {

    // Get the count variables
    int numberOfObservations = m_pObsGroup->GetNumObs();
    int numberOfGroups = m_pObsGroup->GetNumGroups();

    // Send the count variables
    MPI_Send(&numberOfObservations, 1, MPI_INT, workerRank, tag_obsLengthNum, MPI_COMM_WORLD);
    MPI_Send(&numberOfGroups, 1, MPI_INT, workerRank, tag_obsLengthGroup, MPI_COMM_WORLD);

    // Loop over the observations
    for (int entryObservation = 0; entryObservation < numberOfObservations; entryObservation++) {
        // Get the observations from the list
        Observation *obs = m_pObsGroup->GetObsPtr(entryObservation);

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
}

/*
------------------------------------------------------------------------------------------------------------------------------
MPI Communication - ConfigureWorkerParameterGroups()
   Transfers the parameters from from the primary to the secondary worker
------------------------------------------------------------------------------------------------------------------------------
*/
void Algorithm::ConfigureWorkerParameterGroups(int workerRank) {

    
    // Send the total number of parameters
    int numberOfTotalParameters = m_pParamGroup->GetNumParams() + m_pParamGroup->GetNumTiedParams();
    MPI_Send(&numberOfTotalParameters, 1, MPI_INT, workerRank, tag_paramTotalNum, MPI_COMM_WORLD);

    // Real parameters
    int numberOfRealParameters = m_pParamGroup->GetNumRealParams();
    MPI_Send(&numberOfRealParameters, 1, MPI_INT, workerRank, tag_paramTotalReal, MPI_COMM_WORLD);

    for (int entryParameter = 0; entryParameter < numberOfRealParameters; entryParameter++) {
        // Get the parameter from the group, casing to the correct type
        RealParam*param = (RealParam*) m_pParamGroup->GetParamPtr(entryParameter);

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


    // Special parameters
    // TODO: Write special parameters

    // Tied parameters 
    // These will appear as more real parameters on the secondary worker
    int numberOfTiedParameters = m_pParamGroup->GetNumTiedParams();
    MPI_Send(&numberOfTiedParameters, 1, MPI_INT, workerRank, tag_paramTotalReal, MPI_COMM_WORLD);

    for (int entryParameter = 0; entryParameter < numberOfTiedParameters; entryParameter++) {
        // Get the parameter from the group
        TiedParamABC *param = m_pParamGroup->GetTiedParamPtr(entryParameter);

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
    // TODO: write geom parameter

}

/*
------------------------------------------------------------------------------------------------------------------------------
MPI Communication - ConfigureWorkers()
   Transfers all necessary information from the primary to the secondary workers
------------------------------------------------------------------------------------------------------------------------------
*/
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

        // Get the extra files from the directory
        ConfigureWorkerExtraFiles(entryWorker);

        // Get the file pairs for the templates and working files
        ConfigureWorkerFilePairs(entryWorker);

        // Get the observations for comparison
        ConfigureWorkerObservations(entryWorker);

        // Get the parameters for the analysis
        ConfigureWorkerParameterGroups(entryWorker);

    }

    // Do any adjustments on the primary for the secondary workers
    
}

/*
------------------------------------------------------------------------------------------------------------------------------
MPI Communication - SendWorkerParameters()
   Transfers parameters from the primary to the secondary workers
------------------------------------------------------------------------------------------------------------------------------
*/
void Algorithm::SendWorkerParameters(int workerRank, int alternativeIndex, double parametersRegular[]) {


    // Convert the parameters from transformation to log space
    double *parametersConverted = new double[m_pParamGroup->m_NumParams];
    for (int entryParameter = 0; entryParameter < m_pParamGroup->m_NumParams; entryParameter++) {
        // Extract the pointer
        ParameterABC *param = m_pParamGroup->GetParamPtr(entryParameter);

        // Update the parameter with the converted value
        param->SetEstimatedValueTransformed(parametersRegular[entryParameter]);

        // Get the converted value
        parametersConverted[entryParameter] = param->GetTransformedVal();
    }

    // Update the tied parameters
    double *parametersTied = new double[m_pParamGroup->m_NumTied];
    for (int entryParameter = 0; entryParameter < m_pParamGroup->m_NumTied; entryParameter++) {
        // Extract the pointer
        TiedParamABC* param = m_pParamGroup->GetTiedParamPtr(entryParameter);

        // Get the updated value
        parametersTied[entryParameter] = param->GetEstimatedValueTransformed();
    }
 

    // Create a copy of the vector to add the alternative number as the lead value, joining the data together
    int numberOfParamters = m_pParamGroup->m_NumParams + m_pParamGroup->m_NumTied;
    double *parametersTemp = new double[numberOfParamters + 1];
    
    // Fill the array
    parametersTemp[0] = (double)alternativeIndex;
    
    for (int entryParameter = 0; entryParameter < m_pParamGroup->m_NumParams; entryParameter++) {
        parametersTemp[entryParameter + 1] = parametersConverted[entryParameter];
    }

    for (int entryParameter = 0; entryParameter < m_pParamGroup->m_NumTied; entryParameter++) {
        parametersTemp[m_pParamGroup->m_NumParams + entryParameter + 1] = parametersTied[entryParameter];
    }    

    // Send the array to the secondary worker
    MPI_Send(&parametersTemp[0], numberOfParamters + 1, MPI_DOUBLE, workerRank, tag_data, MPI_COMM_WORLD);
}

/*
------------------------------------------------------------------------------------------------------------------------------
MPI Communication - SendWorkerContinue()
   Transfers parameters from the primary to the secondary workers
------------------------------------------------------------------------------------------------------------------------------
*/
void Algorithm::SendWorkerContinue(int workerRank, bool workerContinue) {

    int transferValue = 1;
    if (workerContinue) {
        MPI_Send(&transferValue, 1, MPI_INT, workerRank, tag_continue, MPI_COMM_WORLD);
    }
    else {
        transferValue = 0;
        MPI_Send(&transferValue, 1, MPI_INT, workerRank, tag_continue, MPI_COMM_WORLD);
    }

}


/*
------------------------------------------------------------------------------------------------------------------------------
MPI Communication - ManageSingleObjectiveIterations()
   Manages the solution of a set of parameter alternatives
------------------------------------------------------------------------------------------------------------------------------
*/
void Algorithm::ManageSingleObjectiveIterations(double **parameters, int numberOfParameters, int numberOfAlternatives, double* returnArray) {

    // Get the number of workers in the MPI space
    int numberOfMpiProcesses;
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfMpiProcesses);

    // Begin solving the alternatives
    if (m_bSolveOnPrimary) {
        // Include the primary worker as a solution
        // TODO: create this
    }
    else {
        // Do not solve models on the primary worker
        // Setup the request space
        MPI_Request* request = new MPI_Request[numberOfMpiProcesses - 1];
        MPI_Message mpiMessage;
        MPI_Status mpiStatus;

        int requestFlag = 0;
        int sourceFlag = 0;
        int message = 0;

        // Set the counter variables         
        int sendCounter = 0;
        int receiveCounter = 0;

        // Loop and send work to the secondary workers as they become available
        while (sendCounter < numberOfAlternatives){
            // Get the rank number of a secondary worker waiting for work
            int workerRank, workerValue;

            // TODO: Break each of these into separate functions for reuse
            // Check for any workers with a continue flag 
            MPI_Improbe(MPI_ANY_SOURCE, tag_continue, MPI_COMM_WORLD, &sourceFlag, &message, &mpiStatus);
            if (sourceFlag) {
                // Worker is available to accept work. Go through the steps to send work to it.
                // Receive the continue flag from the worker
                MPI_Mrecv(&workerValue, 1, MPI_INT, &message, &mpiStatus);

                // Extract the worker rank from the status request
                workerRank = mpiStatus.MPI_SOURCE;

                // Transmit a continue flag to the idle secondary worker
                SendWorkerContinue(workerRank, true);

                // Transmit the parameteres to the worker to solve
                SendWorkerParameters(workerRank, sendCounter, parameters[sendCounter]);

                // Increment the send counter
                sendCounter++;
            }

            // Determine if we can receive any values               
            MPI_Improbe(MPI_ANY_SOURCE, tag_data, MPI_COMM_WORLD, &requestFlag, &message, &mpiStatus);
            if (requestFlag) {
                // Get the alternative index and objective function from the secondary worker
                double data[2];
                MPI_Mrecv(&data, 2, MPI_DOUBLE, &message, &mpiStatus);

                // Set the objective value into the objective array
                returnArray[(int)data[0]] = data[1];

                // Increment the solution counter
                receiveCounter++;
            }
        }

        // All alternatives have been sent. Receive all of the solutions
        while (receiveCounter < numberOfAlternatives) {
            // Determine if we can receive any values               
            MPI_Improbe(MPI_ANY_SOURCE, tag_data, MPI_COMM_WORLD, &requestFlag, &message, &mpiStatus);
            if (requestFlag) {
                // Get the alternative index and objective function from the secondary worker
                double data[2];
                MPI_Mrecv(&data, 2, MPI_DOUBLE, &message, &mpiStatus);

                // Set the objective value into the objective array
                returnArray[(int)data[0]] = data[1];

                // Increment the solution counter
                receiveCounter++;

            }
        }    
    }
}


/*
------------------------------------------------------------------------------------------------------------------------------
MPI Communication - TerminateWorkers()
   Tells the workers to stop and begin cleanup
------------------------------------------------------------------------------------------------------------------------------
*/
void Algorithm::TerminateWorkers() {

    // Get the number of workers in the MPI space
    int numberOfMpiProcesses;
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfMpiProcesses);

    for (int entryWorker = 1; entryWorker < numberOfMpiProcesses; entryWorker++) {
        SendWorkerContinue(entryWorker, false);
    }
}





/*****************************************************************************
Bookkeep()
   Performs bookkeeping operations related to parallel executing.
   if bFinal == true, iteration is complete so collect metrics
   if bFinal == false, in the middle of iteration, share information between
    processors
******************************************************************************/
/*void Algorithm::Bookkeep(bool bFinal)
{
    int id, nprocs, temp, i;

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (nprocs == 1) return;

    if (m_bUseSurrogates == true)
    {
        m_pDecision->Bookkeep(bFinal);
    }

    /* ---------------------------------------------------------------
    Collect total evals
    --------------------------------------------------------------- 
    if (bFinal == true)
    {
        for (i = 1; i < nprocs; i++)
        {
            temp = m_Counter;
            MPI_Bcast(&temp, 1, MPI_INTEGER, i, MPI_COMM_WORLD);
            if (id == 0) m_Counter += temp;
        }
    }
} /* end Bookkeep() */


/******************************************************************************
GatherTask()
   Read output file of a SuperMUSE task (stored in pDir directory) and compute
   the associated objective function.
******************************************************************************/
/*double Algorithm::GatherTask(char* pDir)
{
    double val;

    //inc. number of times model has been executed
    m_Counter++;

    //cd to task subdirectory
    MY_CHDIR(pDir);

    //extract computed observations from model output file(s)
    if (m_pObsGroup != NULL) { m_pObsGroup->ExtractVals(); }

    //compute obj. func.
    val = m_pObjFunc->CalcObjFunc();

    //cd out of task subdirectory
    MY_CHDIR("..");

    //ouput results
    Write(val);

    m_CurObjFuncVal = val;

    return (val);
}/* end GatherTask() */


/******************************************************************************
WriteMetrics()
******************************************************************************/
/*void Algorithm::WriteMetrics(FILE* pFile)
{
    if (m_bUseSurrogates == true)
    {
        m_pDecision->WriteMetrics(pFile);
    }
    else
    {
        fprintf(pFile, "Total Evals             : %d\n", m_Counter);
        fprintf(pFile, "Telescoping Strategy    : ");
        switch (m_Telescope)
        {
        case(TSCOPE_PVEX): fprintf(pFile, "convex-power\n"); break;
        case(TSCOPE_CVEX): fprintf(pFile, "convex\n"); break;
        case(TSCOPE_LINR): fprintf(pFile, "linear\n"); break;
        case(TSCOPE_CAVE): fprintf(pFile, "concave\n"); break;
        case(TSCOPE_DCVE): fprintf(pFile, "delayed-concave\n"); break;
        default: fprintf(pFile, "none\n"); break;
        }
        if (m_bCaching == true)
            fprintf(pFile, "Cache Hits              : %d\n", m_NumCacheHits);
        if (m_pParameterCorrection != NULL)
            m_pParameterCorrection->WriteMetrics(pFile);
    }
} /* end WriteMetrics() */

/******************************************************************************
CheckGlobalSensitivity()

Checks that each observation is sensitive to at least one parameter over the
range of possible parameter values. If an observation is not sensitive to any
parameters, a warning will be reported and the observation will be ignored in
the calibration.

Also checks the sensitivity of each parameter. If a given parameter does not
affect the objective function over the entire parameter range, a warning will
be reported and the parameter will be ignored in the calibration.
******************************************************************************/
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

