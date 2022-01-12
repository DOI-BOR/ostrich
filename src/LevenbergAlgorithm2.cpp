/******************************************************************************
The Levenberg-Marquardt algorithm is a hybrid numerical optimization method 
that initially uses the Steepest-Descent technique. However, since it is 
known that the Steepest-Descent algorithm converges very slowly near the 
optimum point, it is desirable to smoothly transition to a polynomial 
approximation method near the optimum. This implementation of the Levenberg 
algorithm is based upon the description/solution provided in the WinPEST 
user's manual, pages 9-42.

******************************************************************************/

// Include the class header
#include "LevenbergAlgorithm2.h"

/******************************************************************************
CTOR()

Initialize data members to reasonable defaults. Certain defaults are then
overridden by the contents of the given confiugration fileName.
******************************************************************************/
LevenbergAlgorithm::LevenbergAlgorithm() {

    // todo: define each of these values
    /*m_Alpha = 0.00;
    m_Beta = 0.00;
    m_Phi = 1E6;
    m_PhiRatio = 1000.00;
    m_PhiRelRed = 1000.00;
    m_Converge = 1E-4;
    m_RatioConv = 0.30;
    m_RelRedConv = 0.01;
    m_Lambda = 10.00;
    m_LamSF = 1.1;
    m_MaxLambdas = 10;
    m_MaxIter = 30;
    m_MoveLimit = 0.10; //limit parameter moves to 10% of current value
    m_NumMS = 1;*/
    //m_pList = NULL;

    // Define each of these parameters
    //m_NumObs = pObsGroup->GetNumObs();
    //m_NumParams = pParamGroup->GetNumParams();

    /*/NEW_PRINT("double", m_NumParams);
    m_pUpgrade = new double[m_NumParams];

    NEW_PRINT("double", m_NumParams);
    m_pTmpVec = new double[m_NumParams];

    NEW_PRINT("double", m_NumObs);
    m_pGamma = new double[m_NumObs];

    NEW_PRINT("double *", m_NumParams);
    m_pPbyP1 = new double* [m_NumParams];

    NEW_PRINT("double *", m_NumParams);
    m_pPbyP2 = new double* [m_NumParams];

    NEW_PRINT("double *", m_NumParams);
    m_pPbyO1 = new double* [m_NumParams];

    NEW_PRINT("double *", m_NumParams);
    m_pPbyO2 = new double* [m_NumParams];

    NEW_PRINT("double *", m_NumParams);
    m_pScale = new double* [m_NumParams];
    MEM_CHECK(m_pScale);*/

    /*for (i = 0; i < m_NumParams; i++) {
        NEW_PRINT("double", m_NumParams);
        m_pPbyP1[i] = new double[m_NumParams];
    }
    for (i = 0; i < m_NumParams; i++) {
        NEW_PRINT("double", m_NumParams);
        m_pPbyP2[i] = new double[m_NumParams];
    }
    for (i = 0; i < m_NumParams; i++) {
        NEW_PRINT("double", m_NumObs);
        m_pPbyO1[i] = new double[m_NumObs];
    }
    for (i = 0; i < m_NumParams; i++) {
        NEW_PRINT("double", m_NumObs);
        m_pPbyO2[i] = new double[m_NumObs];
    }
    for (i = 0; i < m_NumParams; i++) {
        NEW_PRINT("double", m_NumParams);
        m_pScale[i] = new double[m_NumParams];
    }
    MEM_CHECK(m_pScale[i - 1]);*/

    //fill in non-diagonal scale elements
    /*for (i = 0; i < m_NumParams; i++) {
        for (j = 0; j < m_NumParams; j++) {
            m_pScale[i][j] = 0.00;
        }
    }*/


    //set up statistics class   
    //NEW_PRINT("StatsClass", 1);
    //m_pStats = new StatsClass();
    //MEM_CHECK(m_pStats);
    //RegisterStatsPtr(m_pStats);

    //m_pJacob = NULL;
    //m_pJacobUW = NULL;
    //m_pJacobT = NULL;
    //m_pNormal = NULL;


    FILE* pFile;
    char* line;
    char tmp[DEF_STR_SZ];
    char tmp2[DEF_STR_SZ];
    IroncladString inputFilename = GetInFileName();

    //read in algorithm parameters
    pFile = fopen(inputFilename, "r");
    if (pFile != NULL) {
        FindToken(pFile, "EndLevMarAlg", inputFilename);
        rewind(pFile);

        FindToken(pFile, "BeginLevMarAlg", inputFilename);
        line = GetNxtDataLine(pFile, inputFilename);
        while (strstr(line, "EndLevMarAlg") == NULL) {            
            if (strstr(line, "ConvergenceTolerance") != NULL) {
                // Convergence tolerance 
                sscanf(line, "%s %lf", tmp, &m_ObjectiveToleranceMaximum);

            } else if (strstr(line, "StepSize") != NULL) {
                // Initial step size for algorithm
                sscanf(line, "%s %lf", tmp, &m_StepSize);

            } else if (strstr(line, "StepScale") != NULL) {
                // Step size adjustment factor for the algorithm
                sscanf(line, "%s %lf", tmp, &m_StepSizeScaleFactor);

            } else if (strstr(line, "MaximumIterations") != NULL) {
                //maximum number of iterations
                sscanf(line, "%s %d", tmp, &m_NumIterationMaximum);

            } else if (strstr(line, "StepTolerance") != NULL) {
                //maximum number of iterations
                sscanf(line, "%s %lf", tmp, &m_StepToleranceMaximum);

            } else if (strstr(line, "MinimumLambdas") != NULL) {
                //maximum number of iterations
                sscanf(line, "%s %d", tmp, &m_LambdasMinimum);

            } else if (strstr(line, "DerivativeType") != NULL) {
                // Type of derivative to be used
                sscanf(line, "%s %s", tmp, tmp2);

                if (strcmp(tmp2, "FirstForward") == 0) { 
                    m_DerivativeType = FIRST_FORWARD;
                } else if (strcmp(tmp2, "FirstBackward") == 0) { 
                    m_DerivativeType = FIRST_BACKWARD;
                } else if (strcmp(tmp2, "FirstCental") == 0) { 
                    m_DerivativeType = FIRST_CENTRAL; 
                }

            } else {
                sprintf(tmp, "Unknown token: %s", line);
                LogError(ERR_FILE_IO, tmp);
            }

            line = GetNxtDataLine(pFile, inputFilename);
        } 
    }
    else {
        LogError(ERR_FILE_IO, "Using default algorithm setup.");
    }

    fclose(pFile);


    IncCtorCount();
}/* end CTOR */

/******************************************************************************
Destroy()

Frees up the various matrices and vectors used by the algorithm.
******************************************************************************/
void LevenbergAlgorithm::Destroy(void) {
    //ParameterList* pList, * pNext;
    //int i;

    //delete[] m_pUpgrade;
    //delete[] m_pTmpVec;
    //delete[] m_pGamma;

    //for (i = 0; i < m_NumParams; i++) { delete[] m_pPbyP1[i]; }
    //delete[] m_pPbyP1;

    //for (i = 0; i < m_NumParams; i++) { delete[] m_pPbyP2[i]; }
    //delete[] m_pPbyP2;

    //for (i = 0; i < m_NumParams; i++) { delete[] m_pPbyO1[i]; }
    //delete[] m_pPbyO1;

    //for (i = 0; i < m_NumParams; i++) { delete[] m_pPbyO2[i]; }
    //delete[] m_pPbyO2;

    //for (i = 0; i < m_NumParams; i++) { delete[] m_pScale[i]; }
    //delete[] m_pScale;

    ////free up the model backups used in lambda trials
    //delete m_pNonBkup;
    //delete m_pInitBkup;
    //delete m_pDecBkup;
    //delete m_pIncBkup;
    //delete m_pStats;

    //if (m_pList != NULL) {
    //    pList = m_pList;
    //    while (pList->pNxt != NULL) {
    //        pNext = pList->pNxt; //isolate parameter
    //        pList->pNxt = pNext->pNxt; //unlink paramter
    //        //free up parameter
    //        delete[] pNext->p.v;
    //        delete pNext;
    //    }
    //    //free up head of list
    //    delete[] m_pList->p.v;
    //    delete m_pList;
    //}

    IncDtorCount();
}/* end Destroy() */



/******************************************************************************
WarmStart()

Read the best solution from a previous run.
******************************************************************************/
void LevenbergAlgorithm::WarmStart(void) {
    // todo: create this
}

/******************************************************************************
Optimize()

Perform calibration using either Levenberg-Marquardt or GML-MS.
******************************************************************************/
void LevenbergAlgorithm::Optimize(void) {

    // Initialize the workers
    ConfigureWorkers();

    // Define the initial population to solve
    std::vector<std::vector<double>> samples;
    std::vector<bool> lockedParameters;

    // Define the initial jacobian
    std::vector<std::vector<double>> jacobian;
    std::vector<double> lowerValues;
    std::vector<double> upperValues;

    // Create a flag to indicate if in Jacobian or linear search mode
    bool jacobianSolve = false;

    //handle warm start
    if(m_bWarmStart) {
        //WarmStart();
        // Load a previous analysis that was interrupted
        std::cout << "Warm start has not been configured for the Levenberg Algorithm. Exiting the analysis..." << std::endl;
        throw std::invalid_argument("Warm start has not been configured for the Levenberg Algorithm");

    } else {
        // Start a clean analysis
        WriteSetup2("LevenbergAlgorithm", m_ExecCmd, m_pObjFunc->GetObjFuncStr(), m_pParamGroup->GetNumParams(), m_pParamGroup->GetNumTiedParams());
        WriteStartingMetrics(); // todo: write this function

        // Get the values out of the parameter group
        std::vector<double> currentValues;

        for (int entry = 0; entry < m_pParamGroup->GetNumParams(); entry++) {
            ParameterABC* temp = m_pParamGroup->GetParamPtr(entry);
            double test = temp->GetEstimatedValueTransformed();

            currentValues.push_back(temp->GetEstimatedValueTransformed());
            upperValues.push_back(temp->GetUpperBoundTransformed());
            lowerValues.push_back(temp->GetLowerBoundTransformed());
        }

        // Initialize the sample
        CalculateJacobianParameters(currentValues, lowerValues, upperValues, m_StepSize, m_StepToleranceMaximum, m_DerivativeType, samples, lockedParameters);

        // Set the initial condition as the first solve and as the current values
        samples.push_back(currentValues);

        // Need to set the initial condition as the best best alternative
        m_BestAlternative = currentValues;

        // Toggle the jacobian solve flag to indicate to start with jacobian evaluation
        jacobianSolve = true;
    }

    // Enter the solution loop
    std::vector<double> objectivesJacobian;
    bool continueIterations = true;

    while (continueIterations) {

        // Solve for the Jacobian
        if (jacobianSolve) {
            // Initialize the objectives for the solve
            objectivesJacobian = std::vector<double>(samples.size(), INFINITY);

            // Solve the samples
            ManageSingleObjectiveIterations(samples, 0, m_pParamGroup->GetNumParams(), objectivesJacobian);

            // Need to set the first entry as the best if it is iteration zero
            if (m_Iteration == 0) {
                m_BestObjective = objectivesJacobian.back();

                objectivesJacobian.pop_back();
                samples.pop_back();
            }

            // Evaluate for the Jacobian matrix
            CalculateJacobian2(m_BestAlternative, m_BestObjective, objectivesJacobian, samples, m_DerivativeType, lockedParameters, jacobian);

            // Handle any adjustements needed to the jacobian samples
            if (m_DerivativeType == FIRST_CENTRAL) {
                // Keep the upper location as the gradient location
                std::vector<double> objectivesJacobianReduced;
                std::vector<std::vector<double>> samplesLocationsReduced;

                for (int entry = 1; entry < samples.size(); entry += 2) {
                    objectivesJacobianReduced.push_back(objectivesJacobian[entry]);
                    samplesLocationsReduced.push_back(samples[entry]);
                }

                objectivesJacobian = objectivesJacobianReduced;
                samples = samplesLocationsReduced;
            }

            // Write the gradient information to file

            // Toggle the jacobian solve to indicate lambda compute
            jacobianSolve = false;

        }

        // Keep the best output from the gradient for later comparison
        double bestObjectiveGradient;
        int bestObjectiveIndexGradient;

        GetBestSingleObjective(objectivesJacobian, bestObjectiveGradient, bestObjectiveIndexGradient);
        std::vector<double> bestAlternativeGradient = samples[bestObjectiveIndexGradient];


        // Go into a line search for the solution while lambdas produce improvement
        int lambdaCounter = 0;
        double bestLambda;
        int lambdaDirection;
        std::vector<double> objectivesLambda;

        while (!jacobianSolve) {
            // Clear the existing samples
            samples.clear();

            // Create the lambdas
            std::vector<double> lambdas;
            if (lambdaCounter == 0) {
                // Use the initial lamabdas
                lambdas = CreateInitialLambdas();
            } else {
                lambdas = CreateAdditionalLambdas(bestLambda, lambdaDirection);
            }

            for (int entryLambda = 0; entryLambda < lambdas.size(); entryLambda++) {
                // Create a delta vector to hold the parameter correction at the specific delta
                std::vector<double> delta = std::vector<double>(m_pParamGroup->GetNumParams(), 0);

                // Calculate the parameter correction
                for (int entryParameter = 0; entryParameter < m_pParamGroup->GetNumParams(); entryParameter++) {
                    // Calculate the delta for the parameter
                    double parameterDelta = (jacobian[0][entryParameter] * (m_BestObjective - objectivesJacobian[entryParameter])) /
                        (pow(jacobian[0][entryParameter], 2) * (1 + lambdas[entryLambda]));

                    // Append into the delta vector
                    delta[entryParameter] += parameterDelta + m_BestAlternative[entryParameter];
                }

                // Enforce lower parameter bounds
                for (int entry = 0; entry < delta.size(); entry++) {
                    if (delta[entry] < lowerValues[entry]) {
                        delta[entry] = lowerValues[entry];
                    }
                }

                // Enforce upper parameter bounds
                for (int entry = 0; entry < delta.size(); entry++) {
                    if (delta[entry] > upperValues[entry]) {
                        delta[entry] = upperValues[entry];
                    }
                }

                // Push into the sample vector for evalation
                samples.push_back(delta);
            }

            // Evaluate the line search across the deltas
            // Initialize the objectives for the solve
            objectivesLambda = std::vector<double>(samples.size(), INFINITY);

            // Solve the samples
            ManageSingleObjectiveIterations(samples, 0, m_pParamGroup->GetNumParams(), objectivesLambda);

            // Adjust additional lambdas if greater than a single lambda is utilized
            if (lambdas.size() > 1 || m_LambdasMinimum > 1) {
                // Determine if more Lambdas are necessary
                // Find the index of the minimum lambda
                int minLambdaIndex = std::min_element(lambdas.begin() + 1, lambdas.end()) - (lambdas.begin() + 1) + 1;

                // Find the index of the maximum lambda
                int maxLambdaIndex = std::max_element(lambdas.begin() + 1, lambdas.end()) - (lambdas.begin() + 1) + 1;

                // Get the index of the current minimum
                double minimumObjectiveIndex = std::min_element(objectivesLambda.begin(), objectivesLambda.end()) - objectivesLambda.begin();
                int minimumObjectiveValue = *std::min_element(objectivesLambda.begin(), objectivesLambda.end());

                // Set the behavior for the next iteration
                if (*std::min_element(objectivesLambda.begin(), objectivesLambda.end()) > m_BestObjective) {
                    // All solutions are greater than the current best. Break from the lambda loop
                    jacobianSolve = true;
                } else {
                    // Solution has improved. 
                    if (minLambdaIndex == maxLambdaIndex) {
                        // Solution degenerates in the case of 1 element, but the direction has been previously established.
                        // Keep the direction and upate the best lambda value.
                        bestLambda = lambdas[0];

                    } else if (minimumObjectiveIndex == minLambdaIndex && minLambdaIndex != 0) {
                        // Solution is bound by the minimum lambda
                        lambdaDirection = -1;
                        bestLambda = lambdas[minLambdaIndex];

                    } else if (minimumObjectiveIndex == maxLambdaIndex && maxLambdaIndex != 0) {
                        // Solution is bound by the maximum lambda
                        lambdaDirection = 1;                        
                        bestLambda = lambdas[maxLambdaIndex];

                    } else {
                        // Solution is not bound by lambda. Break from the lambda loop.
                        jacobianSolve = true;
                    }
                } 

                // Increment the lambda counter
                lambdaCounter++;

            } else {
                // Direciton has not been established. Transition back to the Jacobian solve.
                jacobianSolve = true;
            }
        }

        // Get the best values from the search
        double bestObjectiveIteration;
        int bestObjectiveIndexIteration;
        
        GetBestSingleObjective(objectivesLambda, bestObjectiveIteration, bestObjectiveIndexIteration);
        std::vector<double> bestAlterativeIteration = samples[bestObjectiveIndexIteration];

        // Compare to the best result from the gradient operation swap, if needed
        if (bestObjectiveGradient < bestObjectiveIteration) {
            bestObjectiveIteration = bestObjectiveGradient;
            bestAlterativeIteration = bestAlternativeGradient;
        }

        // Calculate the tolerance
        m_ObjectiveTolerance = std::abs((bestObjectiveIteration - m_BestObjective) / m_BestObjective);


        // Set the best alternative from the index
        if (bestObjectiveIteration < m_BestObjective) {
            m_BestObjective = bestObjectiveIteration;
            m_BestAlternative = bestAlterativeIteration;
        }
        
        // Adjust the step size if necessary
        if (m_ObjectiveTolerance <= m_ObjectiveToleranceMaximum && m_StepSize >= m_StepToleranceMaximum) {
            // Reduce the step size
            m_StepSize /= m_StepSizeScaleFactor;
        } 


        // Log the iteration
        // Set the values into the pointer groups
        for (int entryParam = 0; entryParam < m_pParamGroup->GetNumParams(); entryParam++) {
            ParameterABC* temp = m_pParamGroup->GetParamPtr(entryParam);
            temp->SetEstimatedValueTransformed(m_BestAlternative[entryParam]);
        }

        // Write iterationr result to the log file
        WriteRecord2(m_Iteration, m_BestObjective, m_ObjectiveTolerance, m_pObsGroup, m_pObjFunc, m_pParamGroup);

        // Increment the generation variable
        m_Iteration++;

        // Determine if the solve should continue
        if (!(m_Iteration < m_NumIterationMaximum && m_StepSize >= m_StepToleranceMaximum)) {
            continueIterations = false;
        } 

        // Regenerate the jacobian sample matrix
        if (continueIterations) {
            // Clear the sample vectors
            samples.clear();

            // Clear the objective vectors
            lockedParameters.clear();
            objectivesJacobian.clear();
            objectivesLambda.clear();

            // Update the samples for the new jacobian location
            CalculateJacobianParameters(m_BestAlternative, lowerValues, upperValues, m_StepSize, m_StepToleranceMaximum, m_DerivativeType, samples, lockedParameters);

            // Toggle the Jacobian solve to true to continue normal solve proceedures
            jacobianSolve = true;
        } 
    }

    // Terminate the secondary workers and let the clean up.
    TerminateWorkers();

    //write algorithm metrics
    WriteEndingMetrics();
}

/*
Create the initial lambdas
*/
std::vector<double> LevenbergAlgorithm::CreateInitialLambdas() {
    // Calculate the lambdas based on the number of compute slots
    int numberOfLambdas;
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfLambdas);

    if (!m_bSolveOnPrimary) {
        numberOfLambdas--;
    }

    // Override the compute slots with the mimimum number, if greater than the sltos
    if (numberOfLambdas < m_LambdasMinimum) {
        numberOfLambdas = m_LambdasMinimum;
    }

    std::vector<double> lambdas = { 0, 0 };
    int counter = 1;
    while (lambdas.size() < numberOfLambdas) {
        // Append the positive version of the counter into the vector
        lambdas.push_back(counter);

        // Append the negative version of the coutner into the vector
        lambdas.push_back(-counter);

        // Increment the counter
        counter++;
    }

    // Adjust the size of the lambdas to handle the one overage that may result from above
    while (lambdas.size() > numberOfLambdas) {
        lambdas.pop_back();
    }

    // Apply powers of two to the lambdas
    for (int entryLambda = 1; entryLambda < lambdas.size(); entryLambda++) {
        lambdas[entryLambda] = pow(2, lambdas[entryLambda]);
    }
    
    // Return to the calling functions
    return lambdas;
}


/*
Create the supplemental lambdas
*/
std::vector<double> LevenbergAlgorithm::CreateAdditionalLambdas(double lambda, int direction) {

    // Get the factor applied to the current lambda as an initial seed
    int factor = (int)(log(lambda) / log(2));

    // Calculate the number of new lambdas to generate basec on the compute slots
    int numberOfLambdas;
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfLambdas);

    if (!m_bSolveOnPrimary) {
        numberOfLambdas--;
    }

    // Make an initial adustment to the factor
    if (direction == -1) {
        factor--;
    } else {
        factor++;
    }

    // Create the additional lambdas
    std::vector<double> lambdas;
    for (int entry = 0; entry < numberOfLambdas; entry++) {
        lambdas.push_back(factor);

        if (direction == -1) {
            factor--;
        } else {
            factor++;
        }
    }

    // Apply powers of two to the lambdas
    for (int entryLambda = 0; entryLambda < lambdas.size(); entryLambda++) {
        lambdas[entryLambda] = pow(2, lambdas[entryLambda]);
    }

    // Return to the calling function
    return lambdas;

}



/******************************************************************************
CopyPoint()

Copy data from one point to another.
******************************************************************************/
/*void LevenbergAlgorithm::CopyPoint(MyPoint* from, MyPoint* to, int np) {
   int i;

   for(i = 0; i < np; i++)
   {
      to->v[i] = from->v[i];
   }
}*//* end CopyPoint() */

/******************************************************************************
GetMinDist()

Compute minimum distance from prospective point to previously evaluated
points.
******************************************************************************/
//double LevenbergAlgorithm::GetMinDist(MyPoint * Point, ParameterList * List)
//{
//   ParameterList * pCur;
//   int i, np;
//   double d, dmin;
//   double v1, v2;
//
//   dmin = NEARLY_HUGE;
//   np = m_pModel->GetParamGroupPtr()->GetNumParams();
//
//   for(pCur = List; pCur != NULL; pCur = pCur->pNxt)
//   {
//      d = 0.00;
//      for(i = 0; i < np; i++)
//      {
//         v1 = Point->v[i];
//         v2 = pCur->p.v[i];
//         d += (v1-v2)*(v1-v2);
//      }/* end for() */
//      d = sqrt(d);
//
//      if(d < dmin) dmin = d;
//   }/* end for() */   
//
//   return dmin;
//}/* end GetMinDist() */

/******************************************************************************
GetRndParamSet()

Generate a random parameter set.
******************************************************************************/
/*void LevenbergAlgorithm::GetRndParamSet(MyPoint * Point)
{
   int i, np;
   double lwr, upr, range, r;

   np = m_pModel->GetParamGroupPtr()->GetNumParams();

   for(i = 0; i < np; i++)
   {
      lwr = m_pModel->GetParamGroupPtr()->GetParamPtr(i)->GetLowerBoundTransformed();
      upr = m_pModel->GetParamGroupPtr()->GetParamPtr(i)->GetUpperBoundTransformed();
      range = upr - lwr;
      r = (double)MyRand() / (double)MY_RAND_MAX;
      Point->v[i] = (r * range) + lwr;
   }
}*//* end GetRndParamSet() */

/******************************************************************************
CalibrateGML()

Perform calibration using Levenberg-Marquardt, equivalent to one application
of GML-MS.
******************************************************************************/
//void LevenbergAlgorithm::CalibrateGML(void)
//{
//   static int GMLcount = 0;
//   StatusStruct pStatus;
//   double oldPhi;
//   double * pMinJac;
//   int done;
//   int id;
//   
//   MPI_Comm_rank(MPI_COMM_WORLD, &id);
//
//   m_CurIter = 0;
//
//   //write banner
//   WriteBanner(m_pModel, "iter  obj. function  ", "lambda");
//
//   m_Phi = m_pModel->Execute();
//   m_pModel->SaveBest(0); //save the input and output files of the best configuration
//   m_BestSavedPhi = m_Phi;
//   InsertParamSet();
//   m_NumEvals++;
// 
//   //write iteration data
//   WriteRecord(m_pModel, 0, m_Phi, m_Lambda);
//   pStatus.curIter = 0;
//   pStatus.maxIter = m_MaxIter;
//   pStatus.pct = (((float)100.00*(float)GMLcount)/(float)m_NumMS);
//   pStatus.numRuns = m_pModel->GetCounter();
//   WriteStatus(&pStatus);
//
//   //main loop, iterates using Levenberg-Marquardt alg.
//   done = 0;
//   while(done == 0)
//   {
//      if(IsQuit() == true)
//      { 
//         done = 1;
//         MPI_Bcast(&done, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
//         break;
//      }
//
//      m_CurIter++;
//
//      /* ----------------------------------------------------
//      Calculate the Jacobian matrix, possibly in parallel.
//      ------------------------------------------------------ */
//      CalcJacobian();
//      pMinJac = m_pStats->GetMinJac();
//
//      if(id == 0)
//      {
//         CalcNormal(); //calculate the normal matrix
//         CalcScale(); //calculate the scale matrix
//         //determine best lambda for current iteration
//         oldPhi = m_Phi;
//         AdjustLambda();
//
//         //check move against the best Jacobian evaluation
//         if(pMinJac[0] < m_Phi)
//         {
//            m_Phi = pMinJac[0];
//            //use the min Jacobian data to adjust model and parameter and observation groups
//            m_pModel->SetObjFuncVal(pMinJac[0]);
//            m_pModel->GetParamGroupPtr()->WriteParams(&(pMinJac[1]));
//            m_pModel->GetObsGroupPtr()->WriteObservations(&(pMinJac[1+m_NumParams]));
//         }/* end if() */
//            
//         //check for convergence
//         m_PhiRatio = m_Phi / oldPhi;
//         m_PhiRelRed = fabs(1.00 - m_PhiRatio);
//
//         //write iteration data
//         WriteRecord(m_pModel, m_CurIter, m_Phi, m_Lambda);
//         pStatus.curIter = m_CurIter;
//         pStatus.pct = (((float)100.00*(float)GMLcount)/(float)m_NumMS) + 
//                       (((float)100.00*(float)m_CurIter)/((float)m_MaxIter*(float)m_NumMS));
//         pStatus.numRuns = m_pModel->GetCounter();
//         WriteStatus(&pStatus);
//
//         if ((m_CurIter >= m_MaxIter) ||
//            (m_PhiRelRed < m_Converge) ||
//            ((oldPhi - m_Phi) < m_Converge))
//         {
//            done = 1;
//            pStatus.pct = (((float)100.00*(float)GMLcount+1)/(float)m_NumMS);
//         }/* end if() */
//      }/* end if(master processor) */
//      MPI_Bcast(&done, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
//
//      //perform intermediate bookkeeping
//      m_pModel->Bookkeep(false);
//   }/* end while() */
//   m_NumIters = m_CurIter;
//   GMLcount++;
//}/* end CalibrateGML() */

/******************************************************************************
CalcJacobian()

Calculate the Jacobian matrix (matrix of dObs/dParams) using the algorithm
in the Statistics Class.
******************************************************************************/
/*void LevenbergAlgorithm::CalcJacobian(void)
{
   m_pJacob   = m_pStats->CalcJacobian(&m_BestSavedPhi);
   m_pJacobT  = m_pStats->GetJacobT();
   m_pJacobUW = m_pStats->GetJacobUW();
   m_pStats->AdjustJacobian();
}*//* end CalcJacobian() */

/******************************************************************************
CalcNormal()

Calculate what is known as the 'normal' regression matrix:
   (J^T)*Q*J.
******************************************************************************/
/*void LevenbergAlgorithm::CalcNormal(void)
{   
   m_pNormal = m_pStats->CalcNormal();
}*//* end CalcNormal() */

/******************************************************************************
CalcScale()

Calculate the scaling matrix. It is an all diagonal matrix that scales the
'normal' matrix so as to avoid numerical round-off errors and instability 
problems.
******************************************************************************/
//void LevenbergAlgorithm::CalcScale(void)
//{
//   int i, p;   
//
//   p = m_NumParams - m_pStats->GetNumHeldParams();
//   for(i = 0; i < p; i++)
//   {
//      m_pScale[i][i] = (1.00 / sqrt(m_pNormal[i][i]));
//   }/* end for() */  
//}/* end CalcScale() */

/******************************************************************************
AdjustLambda()

Modify lambda in various ways to determine the best lambda for the current 
iteration. Reducing lambda is favored, but increasing lambda and leaving lambda
at its current value are also tested.
******************************************************************************/
//void LevenbergAlgorithm::AdjustLambda(void)
//{   
//   int iter;
//   double oldPhi, phiConst, phiDec, phiInc, phiTry;
//   double lamConst, lamDec, lamInc, lamTry;   
//      
//   lamConst = m_Lambda;
//   lamDec   = m_Lambda / m_LamSF;
//   lamInc   = m_Lambda * m_LamSF;
//
//   //display banner
//   WriteInnerEval(WRITE_LEV, m_MaxLambdas, '.');
//   
//   //Compute initial lambda effects
//   WriteInnerEval(1, m_MaxLambdas, '.');
//   m_pInitBkup->Store();
//   phiConst =  TryLambda(lamConst); //non-adjusted lambda trial
//   m_pNonBkup->Store();
//
//   WriteInnerEval(2, m_MaxLambdas, '-');
//   m_pInitBkup->SemiRestore();
//   phiDec = TryLambda(lamDec); //decreased lambda trial
//   m_pDecBkup->Store(); 
//
//   WriteInnerEval(3, m_MaxLambdas, '+');
//   m_pInitBkup->SemiRestore();
//   phiInc = TryLambda(lamInc); //increased lambda trial
//   m_pIncBkup->Store();
//   m_pInitBkup->SemiRestore();
//
//   iter = 3;
//
//   //check to see if none of the lambda adjustments were effective
//   if((m_Phi < phiConst) && (m_Phi < phiDec) && (m_Phi < phiInc))
//   {
//      m_Lambda /= m_LamSF;
//      WriteInnerEval(WRITE_ENDED, m_MaxLambdas, 'n');
//      return;
//   }/* end if() */
//
//   /*------------------------------------------------
//   Decreasing lambda caused obj. func. to decrease...
//   and more so than a constant or increasing lambda
//   ------------------------------------------------*/
//   if((phiDec < m_Phi) && (phiDec <= phiInc) && (phiDec <= phiConst))
//   {         
//      oldPhi = m_Phi;
//      m_pDecBkup->SemiRestore();
//      while(iter < m_MaxLambdas)
//      {
//         //converged?
//         m_PhiRatio = phiDec / oldPhi;
//         m_PhiRelRed = 1.00 - m_PhiRatio;
//         if((m_PhiRatio < m_RatioConv) || 
//            ((m_PhiRelRed < m_RelRedConv) && (m_PhiRelRed > 0.00)))
//         {
//            //update phi, lambda
//			m_pDecBkup->SemiRestore();
//            m_Phi = phiDec;
//            m_Lambda = lamDec;
//            WriteInnerEval(WRITE_ENDED, m_MaxLambdas, 'c');
//            return;
//         }/* end if() */
//
//         //try decreasing lambda
//         WriteInnerEval(iter+1, m_MaxLambdas, '-');
//         lamTry = lamDec / m_LamSF;
//         oldPhi = phiDec;
//         phiTry = TryLambda(lamTry);
//         
//         if(phiTry < oldPhi) //obj. reduced, accept trial move
//         {   
//            phiDec = phiTry;
//            lamDec = lamTry;            
//            m_pDecBkup->Store();            
//         }/* end if() */
//         else { break;} //failed to reduce, reject move and exit loop.
//         iter++;         
//      }/* end while */      
//   }/* end if() */
//
//   /*------------------------------------------------
//   Increasing lambda caused obj. func. to decrease...
//   and more so than a constant or decreasing lambda
//   ------------------------------------------------*/
//   if((phiInc < m_Phi) && (phiInc <= phiDec) && (phiInc <= phiConst))
//   {         
//      oldPhi = m_Phi;
//      m_pIncBkup->SemiRestore();
//      while(iter < m_MaxLambdas)
//      {
//         //converged?
//         m_PhiRatio = phiInc / oldPhi;
//         m_PhiRelRed = 1.00 - m_PhiRatio;
//         if((m_PhiRatio < m_RatioConv) || 
//            ((m_PhiRelRed < m_RelRedConv) && (m_PhiRelRed > 0.00)))
//         {
//            //update phi, lambda            
//			m_pIncBkup->SemiRestore();
//            m_Phi = phiInc;
//            m_Lambda = lamInc;
//            WriteInnerEval(WRITE_ENDED, m_MaxLambdas, 'c');
//            return;
//         }/* end if() */
//
//         //try increasing lambda
//         WriteInnerEval(iter+1, m_MaxLambdas, '+');
//         lamTry = lamInc * m_LamSF;
//         oldPhi = phiInc;
//         phiTry = TryLambda(lamTry);
//         
//         if(phiTry < oldPhi) //obj. reduced, accept trial move
//         {
//            phiInc = phiTry;
//            lamInc = lamTry;            
//            m_pIncBkup->Store();            
//         }/* end if() */
//         else { break;}//failed to reduce, reject move and exit loop.
//         iter++;         
//      }/* end while */      
//   }/* end if() */
//
//   /*
//   Didn't converge on a lambda, but some lambda(s) did reduce obj. func.
//   Use the lambda that had the best result.
//   */
//   if((phiDec <= phiConst) && (phiDec <= phiInc))
//   {     
//      m_pDecBkup->SemiRestore();
//      m_Phi = phiDec;
//      m_Lambda = lamDec;
//   }/* end if() */
//   else if((phiConst <= phiDec) && (phiConst <= phiInc))
//   {      
//      m_pNonBkup->SemiRestore();
//      m_Phi = phiConst;
//      m_Lambda = lamConst;
//   }/* end else if() */
//   else
//   {    
//      m_pIncBkup->SemiRestore();
//      m_Phi = phiInc;
//      m_Lambda = lamInc;
//   }/* end else() */
//
//   WriteInnerEval(WRITE_ENDED, m_MaxLambdas, 'y');
//}/* end AdjustLambda() */

/******************************************************************************
TryLambda()

Using the lambda argument, compute alpha (the Marquardt parameter) and the 
associated upgrade vector. Then apply this upgrade vector to the current set
of model parameters and execute the model.

Returns: The value of the obj. function, as returned by the Model class.
******************************************************************************/
//double LevenbergAlgorithm::TryLambda(double lambda)
//{
//   double phi;
//
//   //fill in vector of residuals [needed by CalcUpgrade() and CalcBeta()]
//   m_pResid = m_pStats->CalcResiduals();
//   m_pStats->AdjustResiduals();
//
//   //solve for alpha (Marquardt parameter)
//   CalcAlpha(lambda);
//   //calculate the upgrade vector (the optimal direction)
//   CalcUpgrade();
//   //calculate the gamma vector (used in beta calculation)
//   CalcGamma();
//   /* Calculate beta, a factor that scales the upgrade vector 
//   to achieve optimal step size */
//   CalcBeta();
//   //Update model parameters using beta and upgrade vector
//   AdjModelParams();
//   //Evaluate objective function at revised location
//   phi = m_pModel->Execute();
//   if(phi < m_BestSavedPhi)
//   {
//      m_pModel->SaveBest(0);
//      m_BestSavedPhi = phi;
//   }
//   InsertParamSet();
//   m_NumEvals++;
//   return phi;
//}/* end TryLambda() */

/******************************************************************************
CalcAlpha()

Compute alpha, the Marquardt parameter, based on the given lambda value. 
By definition, alpha = lambda/max(Si), where max(Si) is the maximum value of 
the diagonal elements of the scale matrix.
******************************************************************************/
//void LevenbergAlgorithm::CalcAlpha(double lambda)
//{
//   int i, p;
//   double max;
//   double cur;
//
//   max = m_pScale[0][0]*m_pScale[0][0];
//   p = m_NumParams - m_pStats->GetNumHeldParams();
//   for(i = 0; i < p; i++)
//   {
//      cur = m_pScale[i][i]*m_pScale[i][i];
//      if(cur > max) {max = cur;}
//   }/* end for() */
//   
//   m_Alpha = lambda/max;
//}/* end CalcAlpha() */

/******************************************************************************
CalcUpgrade()

Compute the upgrade vector using a sequence of matrix and vector 
multiplication and inversion, as shown on page 20 of the WinPEST Manual.

Adjusted for Box-Cox Transformations ---- weights are incorporated into 
residuals prior to transformation.
******************************************************************************/
//void LevenbergAlgorithm::CalcUpgrade(void)
//{
//   int i, n, p;
//
//   p = m_NumParams - m_pStats->GetNumHeldParams();
//   n = m_NumObs - m_pStats->GetNumHeldObs();
//
//   //compute (JS)^T=S^T*J^T=S*J^T, where S^T = S, since S is a diag. matrix.
//   //S is [pXp], J^T is [pXo], and S*J^T is [pXo]
//   MatMult(m_pScale, m_pJacobT, m_pPbyO1, p, p, n);
//   MatMult(m_pPbyO1, m_pJacob, m_pPbyP1, p, n, p);
//   //multiply result by S [pxp], result is [pxp]
//   MatMult(m_pPbyP1, m_pScale, m_pPbyP2, p, p, p);
//   //add alpha*S^T*S to the result of previous step
//   for(i = 0; i < p; i++)
//   { m_pPbyP2[i][i] += (m_pScale[i][i]*m_pScale[i][i]*m_Alpha); }
//   //invert the result of the previous step   
//   MatInv(m_pPbyP2, m_pPbyP1, p);
//   //multiply iverse by S*J^T, which is presently stored in m_pPbyO1
//   MatMult(m_pPbyP1, m_pPbyO1, m_pPbyO2, p, p, n);
//   VectMult(m_pPbyO2, m_pResid, m_pTmpVec, p, n);
//
//   //premultiply result by S [pXp] and store in upgrade vector [pX1]
//   VectMult(m_pScale, m_pTmpVec, m_pUpgrade, p, p);  
//}/* end CalcUpgrade() */

/******************************************************************************
CalcGamma()

Gamma is used in the computation of beta, the optimum step size. Formulation 
of gamma is given on page 21 of the WinPEST Manual.
******************************************************************************/
//void LevenbergAlgorithm::CalcGamma(void)
//{
//   int p, n;
//   p = m_NumParams - m_pStats->GetNumHeldParams();
//   n = m_NumObs - m_pStats->GetNumHeldObs();
//
//   //gamma [ox1] = J [oxp] * u [px1]
//   VectMult(m_pJacobUW, m_pUpgrade, m_pGamma, n, p);
//}/* end CalcGamma() */

/******************************************************************************
CalcBeta()

Beta is the optimum step size for the direction specified by the upgrade 
vector. Formulation of beta is given on page 21 of the WinPEST Manual.
******************************************************************************/
//void LevenbergAlgorithm::CalcBeta(void)
//{
//   int i, n;
//   double numer;
//   double denom;
//   double wt;
//
//   n = m_NumObs - m_pStats->GetNumHeldObs();
//   numer = 0.00;
//   denom = 0.00;
//
//   for(i = 0; i < n; i++)
//   {
//      wt=GetObsWeight(m_pModel->GetObsGroupPtr()->GetObsPtr(i));
//      numer += (m_pResid[i]*m_pGamma[i]*wt);
//      denom += (m_pGamma[i]*m_pGamma[i]*wt*wt);
//   }/* end for() */
//
//   m_Beta = (numer/denom);   
//}/* end CalcBeta() */

/******************************************************************************
AdjModelParams()

Modify model parameters using beta (step size) and upgrade vector (direction).
******************************************************************************/
//void LevenbergAlgorithm::AdjModelParams(void)
//{
//   ParameterGroup * pGroup;
//   ParameterABC * pParam;
//   double upr, lwr;
//   double oldVal;
//   double curVal;
//   double adjust;
//   double range;
//   double maxAdj;
//   int i;
//
//   //must adjust the upgrade vector to account for held parameters
//   m_pStats->AdjustVector(m_pUpgrade, false);
//
//   pGroup = m_pModel->GetParamGroupPtr();
//   for(i = 0; i < m_NumParams; i++)
//   {
//      pParam = pGroup->GetParamPtr(i);
//      oldVal = pParam->GetEstimatedValueTransformed();
//      upr = pParam->GetUpperBoundTransformed();
//      lwr = pParam->GetLowerBoundTransformed();
//      range = upr - lwr;
//      maxAdj = range * m_MoveLimit;
//      /*
//      The optimal parameter adj. (not considering move limits) is 
//      defined as: beta * (ugrade vector).
//      */
//      adjust = m_Beta * m_pUpgrade[i];
//
//      /* 
//      Check the optimal adjustment against move limits to prevent large 
//      changes in parameters. We want to prevent large jumps because the 
//      numerical solution relies on a Taylor Series expansion about the 
//      current set of parameters. This is a linear approximation and 
//      will only be valid in the proximity of the current paramter set.
//      */      
//      if(fabs(adjust) > maxAdj)
//      {
//         if(adjust > 0) {adjust = maxAdj;}
//         else {adjust = -1.00 * maxAdj;}
//         m_NumMoveViols++;
//      } /* end if() */
//
//      curVal = oldVal + adjust; 
//   
//      //if move exceeds limits, only go half the distance to the limit
//      if(curVal <= lwr){ curVal = (oldVal+lwr)/2.00; m_NumLwrViols++;}
//      if(curVal >= upr){ curVal = (oldVal+upr)/2.00; m_NumUprViols++;}
//
//      pParam->SetEstimatedValueTransformed(curVal);
//   }/* end for() */   
//}/* end AdjModelParams() */

/******************************************************************************
InsertParamSet()

Insert the most recently evaluated parameter set into the list.
******************************************************************************/
//void LevenbergAlgorithm::InsertParamSet(void)
//{
//   if(m_bMS == false) return;
//
//   ParameterList * pTmp;
//
//   if(m_pList == NULL)
//   {
//      //allocate space for new list entry
//      NEW_PRINT("ParameterList", 1);
//      m_pList = new ParameterList;
//      MEM_CHECK(m_pList);
//
//      NEW_PRINT("double", m_NumParams);
//      m_pList->p.v = new double[m_NumParams];
//      MEM_CHECK(m_pList->p.v)
//
//      m_pModel->GetParamGroupPtr()->ReadParams(m_pList->p.v);
//      m_pList->pNxt = NULL;
//   }/* end if(first entry) */
//   else
//   {
//      for(pTmp = m_pList; pTmp->pNxt != NULL; pTmp = pTmp->pNxt)
//      {
//         //no-op, just advancing to end of list
//      }
//
//      //allocate space for new list entry
//      NEW_PRINT("ParameterList", 1);
//      pTmp->pNxt = new ParameterList;
//      MEM_CHECK(pTmp->pNxt);
//      pTmp = pTmp->pNxt;
//
//      NEW_PRINT("double", m_NumParams);
//      pTmp->p.v = new double[m_NumParams];
//      MEM_CHECK(pTmp->p.v)
//
//      m_pModel->GetParamGroupPtr()->ReadParams(pTmp->p.v);
//      pTmp->pNxt = NULL;
//   }/* end else(not first entry) */   
//}/* end InsertParamSet() */


/******************************************************************************
WriteMetrics()

Write out algorithm metrics and setup.
******************************************************************************/
void LevenbergAlgorithm::WriteStartingMetrics(void) {

    // Open the log file
    char fileName[DEF_STR_SZ];
    FILE* pFile;
    sprintf(fileName, "OstOutput%d.txt", 0);
    pFile = fopen(fileName, "a");

    fprintf(pFile, "Algorithm         : Levenberg-Marquardt\n");
    fprintf(pFile, "Max Iterations    : %d\n", m_NumIterationMaximum);
    fprintf(pFile, "Objective Tolerance : %d\n", m_ObjectiveToleranceMaximum);
    fprintf(pFile, "Step Tolerance   : %lf\n", m_StepToleranceMaximum);

    // Write the parameter to the file
    fprintf(pFile, "Run   "); // TODO: align columns 
    fprintf(pFile, "%-12s  ", "Objective");
    m_pParamGroup->Write(pFile, WRITE_BNR);
    fprintf(pFile, "Convergence\n");

    // Close the log file
    fclose(pFile);
   
}

/**************************************************************************************************************************************************************
WriteEndingMetrics()

Write out setup and metrics for the algorithm.
**************************************************************************************************************************************************************/
void LevenbergAlgorithm::WriteEndingMetrics(void) {
    // TODO: Move htis to the write utility class
    // TODO: Update to C++

    // Open the log file
    char fileName[DEF_STR_SZ];
    FILE* pFile;
    sprintf(fileName, "OstOutput%d.txt", 0);
    pFile = fopen(fileName, "a");

    // Write the result information
    WriteOptimalToFileWithGroup2(pFile, m_pParamGroup, m_BestObjective);

    // Write the algorithm information
    fprintf(pFile, "\nAlgorithm Metrics\n");
    fprintf(pFile, "Algorithm               : Levenberg-Marquardt\n");
    fprintf(pFile, "Actual Objective Convergence  : %E\n", m_ObjectiveTolerance);
    fprintf(pFile, "Actual Step Convergence       : %lf\n", m_StepSize * m_StepSizeScaleFactor);
    fprintf(pFile, "Actual Iterations      : %d\n", m_Iteration);

    if (m_ObjectiveTolerance <= m_ObjectiveToleranceMaximum) {
        fprintf(pFile, "Algorithm successfully converged on a solution\n");
    } else {
        fprintf(pFile, "Algorithm failed to converge on a solution, more iterations may be needed\n");
    }

    // Close the log file
    fclose(pFile);
}

/******************************************************************************
LEV_Program()

Create a model and solve using the Levenber-Marquardt algorithm.
******************************************************************************/
void LEV_Program(int argc, StringType argv[]) {
   
    // Get the rank of the process 
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Determine the behavior based on rank
    if (rank == 0) {
        // This is the primary worker. Setup the algorithm and call the processing routine
        // Create the LM algorithm
        LevenbergAlgorithm* LM = new LevenbergAlgorithm();

        // Begin the primary worker operations. This transfers information to the secondary workers and manages the solution process
        LM->Optimize();

    } else {
        // This is a secondary worker. Create the secondary worker and call the work routine
        // Create and setup the secondary worker, telling it to configure from MPI
        ModelWorker worker = ModelWorker(true);

        // Start the work in the processor 
        worker.WorkMPI();

        // Tear-down is triggered by the primary worker. Perform any additional cleanup actions
    }

} 

/******************************************************************************
GMLMS_Program()

Create a model and solve using the Levenber-Marquardt algorithm with 
multi-starts.
******************************************************************************/
// Todo: update this function
//void GMLMS_Program(int argc, StringType argv[])
//{
//   NEW_PRINT("Model", 1);
//   ModelABC * model = new Model;
//   MEM_CHECK(model);
//
//   if(model->GetObjFuncId() == OBJ_FUNC_WSSE) 
//   { 
//      NEW_PRINT("LevenbergAlgorithm", 1);
//      LevenbergAlgorithm * LA = new LevenbergAlgorithm(model, true);
//      MEM_CHECK(LA);
//
//      LA->Calibrate(); 
//      delete LA;
//   }
//   else 
//   { 
//      printf("GML-MS algorithm can only be used for regression.\n");
//   }/* end else() */
//   
//   delete model;
//} /* end GMLMS_Program() */
