/**************************************************************************************************************************************************************

Computes the gradient of the objective function using forward finite differences, for use by gradient-based optimization algorithms. Parameter sets are
evaluated through a user-supplied evaluator so that the perturbations of all parameters can be solved concurrently using the asynchronous MPI work-farming
functionality provided by the Algorithm base class.

**************************************************************************************************************************************************************/

// Include the class header
#include "OptMathClass2.h"

/**************************************************************************************************************************************************************
CTOR

Initializes the finite difference settings. All settings are assigned default values (range-relative increments of 0.001) and then the user input file is
checked for overriding values.
**************************************************************************************************************************************************************/
OptMathClass2::OptMathClass2(ParameterGroup* pParamGroup, Evaluator evaluator) {

    m_pParamGroup = pParamGroup;
    m_Evaluate = evaluator;
    m_NumParams = m_pParamGroup->GetNumParams();
    m_DiffInc = std::vector<double>(m_NumParams, 0.001);

    InitFromFile(GetInFileName());

    IncCtorCount();
}

/**************************************************************************************************************************************************************
Destroy()

Free up memory used by the class.
**************************************************************************************************************************************************************/
void OptMathClass2::Destroy(void) {

    m_DiffInc.clear();

    IncDtorCount();
}

/**************************************************************************************************************************************************************
WriteMetrics()

Reports on the setup of the finite difference calculations and also various run-time metrics.
**************************************************************************************************************************************************************/
void OptMathClass2::WriteMetrics(FILE* pFile) {

    fprintf(pFile, "\nFinite Difference Metrics\n");
    fprintf(pFile, "Difference Type    : Forward\n");
    fprintf(pFile, "Increment Type     : ");
    if (m_DiffIncType == FD_VALUE_REL) { fprintf(pFile, "Value-Relative\n"); }
    else if (m_DiffIncType == FD_ABSOLUTE) { fprintf(pFile, "Absolute\n"); }
    else { fprintf(pFile, "Range-Relative\n"); }

    fprintf(pFile, "Finite Difference Increments\n");
    for (int entryParameter = 0; entryParameter < m_NumParams; entryParameter++) {
        fprintf(pFile, "%-12s : %E\n", m_pParamGroup->GetParamPtr(entryParameter)->GetName(), m_DiffInc[entryParameter]);
    }
    fprintf(pFile, "Finite Difference Mimumum Increment : %E\n", m_MinInc);
    fprintf(pFile, "Gradient Evals     : %d\n", m_GradCount);
    fprintf(pFile, "Derivative Evals   : %d\n", m_DiffCount);
}

/**************************************************************************************************************************************************************
InitFromFile()

Reads the finite difference configuration from the BeginMathAndStats group of the input file.
**************************************************************************************************************************************************************/
void OptMathClass2::InitFromFile(IroncladString pFileName) {

    FILE* pFile;
    char* line;
    char tmp[DEF_STR_SZ];
    char tmp2[DEF_STR_SZ];
    char* pTmp;
    int i, j;

    // Read in the configuration
    pFile = fopen(pFileName, "r");
    if (pFile == NULL) {
        // Couldn't open file, use defaults and log the error
        LogError(ERR_FILE_IO, "Couldn't open math config. file. Using Defaults");
        return;
    }

    // Make sure the correct tokens are present
    if (CheckToken(pFile, "BeginMathAndStats", pFileName) == true) {
        // Identify the end of the block
        FindToken(pFile, "EndMathAndStats", pFileName);
        rewind(pFile);

        // Import data starting from the beginning of the block
        FindToken(pFile, "BeginMathAndStats", pFileName);
        line = GetNxtDataLine(pFile, pFileName);

        while (strstr(line, "EndMathAndStats") == NULL) {
            if (strstr(line, "DiffIncType") != NULL) {
                // Set the increment type
                sscanf(line, "%s %s", tmp, tmp2);
                MyStrLwr(tmp2);
                if (strstr(tmp2, "range-relative") != NULL) { m_DiffIncType = FD_RANGE_REL; }
                else if (strstr(tmp2, "value-relative") != NULL) { m_DiffIncType = FD_VALUE_REL; }
                else if (strstr(tmp2, "absolute") != NULL) { m_DiffIncType = FD_ABSOLUTE; }
                else if (strstr(tmp2, "optimal") != NULL) {
                    // Optimal increments are not supported by the parallel implementation
                    LogError(ERR_FILE_IO, "Optimal finite difference increments are not supported. Using range-relative increments.");
                    m_DiffIncType = FD_RANGE_REL;
                }
            }
            else if (strstr(line, "DiffType") != NULL) {
                // Only forward differences are supported. The entry is accepted for compatibility.
            }
            else if (strstr(line, "DiffRelIncrement") != NULL) {
                // Set the range-relative increments
                MyStrLwr(line);
                pTmp = line + (int)strlen("DiffRelIncrement");

                for (i = 0; i < m_NumParams; i++) {
                    j = ExtractString(pTmp, tmp);
                    m_DiffInc[i] = atof(tmp);
                    if (j == -1) { break; }
                    pTmp += j;
                }

                // If only a partial list, apply the first value to the remaining parameters
                for (i = i; i < m_NumParams; i++) { m_DiffInc[i] = m_DiffInc[0]; }

                // This keyword (DiffRelIncrement) is range-relative
                m_DiffIncType = FD_RANGE_REL;
            }
            else if (strstr(line, "DiffIncrement") != NULL) {
                // Set the increments
                MyStrLwr(line);
                pTmp = line + (int)strlen("DiffIncrement");

                for (i = 0; i < m_NumParams; i++) {
                    j = ExtractString(pTmp, tmp);
                    m_DiffInc[i] = atof(tmp);
                    if (j == -1) { break; }
                    pTmp += j;
                }

                // If only a partial list, apply the first value to the remaining parameters
                for (i = i; i < m_NumParams; i++) { m_DiffInc[i] = m_DiffInc[0]; }
            }
            else if (strstr(line, "DiffMinIncrement") != NULL) {
                // Set the minimum increment
                sscanf(line, "%s %lf", tmp, &m_MinInc);
            }
            else {
                sprintf(tmp, "Unknown token: %s", line);
                LogError(ERR_FILE_IO, tmp);
            }

            // Get the next line
            line = GetNxtDataLine(pFile, pFileName);
        }
    }

    fclose(pFile);
}

/**************************************************************************************************************************************************************
CalcIncrement()

Computes the finite difference increment for the parameter at <index>, which currently has a value of <value>. A forward step is used unless it would exceed
the upper bound, in which case the step is reversed. Returns the signed increment.
**************************************************************************************************************************************************************/
double OptMathClass2::CalcIncrement(int index, double value, FiniteDiffIncType incrementType) {

    double upr = m_pParamGroup->GetParamPtr(index)->GetUpperBoundTransformed();
    double lwr = m_pParamGroup->GetParamPtr(index)->GetLowerBoundTransformed();
    double dx;

    if (incrementType == FD_VALUE_REL) { dx = MyMax(fabs(m_DiffInc[index] * value), m_MinInc); }
    else if (incrementType == FD_ABSOLUTE) { dx = fabs(m_DiffInc[index]); }
    else { dx = fabs(m_DiffInc[index] * (upr - lwr)); }

    // Trick from NR in C to ensure the increment is exactly representable
    double next = value + dx;
    dx = next - value;

    // Avoid exceeding the upper parameter limit by switching direction
    if ((value + dx) > upr) { dx *= -1.00; }

    return dx;
}

/**************************************************************************************************************************************************************
CalcDerivatives()

Computes the partial derivatives of the parameters listed in <indices> at <point>, which has an objective of <objective>. The perturbed parameter sets are
evaluated together and the results are stored into <gradient>. If a better minimum than fmin is found, fmin and xmin are replaced with the new minimum value
and corresponding parameter values.
**************************************************************************************************************************************************************/
void OptMathClass2::CalcDerivatives(std::vector<int> indices, std::vector<double> point, double objective, FiniteDiffIncType incrementType,
                                    std::vector<double>& gradient, double& fmin, std::vector<double>& xmin) {

    // Create the perturbed parameter sets. A zero increment cannot be used, so the derivative is set to zero and the set is flagged as invalid (empty).
    std::vector<double> increments;
    std::vector<std::vector<double>> samples;
    for (int entryIndex = 0; entryIndex < (int)indices.size(); entryIndex++) {
        int index = indices[entryIndex];
        double dx = CalcIncrement(index, point[index], incrementType);
        increments.push_back(dx);

        if (dx == 0.00) {
            samples.push_back(std::vector<double>());
        }
        else {
            std::vector<double> sample = point;
            sample[index] = point[index] + dx;
            samples.push_back(sample);
        }
    }

    // Evaluate the perturbed parameter sets together
    std::vector<double> objectives = m_Evaluate(samples);

    // Compute the derivatives and update the optimal, if appropriate
    for (int entryIndex = 0; entryIndex < (int)indices.size(); entryIndex++) {
        int index = indices[entryIndex];

        if (samples[entryIndex].empty()) {
            gradient[index] = 0.00;
            continue;
        }

        m_DiffCount++;
        gradient[index] = (objectives[entryIndex] - objective) / increments[entryIndex];

        if (objectives[entryIndex] < fmin) {
            fmin = objectives[entryIndex];
            xmin = samples[entryIndex];
        }
    }
}

/**************************************************************************************************************************************************************
CalcGradient()

Calculates the gradient at <point>, which has an objective of <objective>, using forward differences. The perturbations of all parameters are evaluated
together. Any derivative that is ~0.00 is recomputed using a range-relative increment (if a different increment type was in use).

If a better minimum than fmin is found, fmin and xmin are replaced with the new minimum value and corresponding parameter values.
**************************************************************************************************************************************************************/
std::vector<double> OptMathClass2::CalcGradient(std::vector<double> point, double objective, double& fmin, std::vector<double>& xmin) {

    std::vector<double> gradient(m_NumParams, 0.00);

    // Compute the partial derivatives of all parameters together
    std::vector<int> indices;
    for (int entryParameter = 0; entryParameter < m_NumParams; entryParameter++) { indices.push_back(entryParameter); }
    CalcDerivatives(indices, point, objective, m_DiffIncType, gradient, fmin, xmin);

    // If an FD calculation is ~0.00, retry using a range-relative increment
    if (m_DiffIncType != FD_RANGE_REL) {
        std::vector<int> retryIndices;
        for (int entryParameter = 0; entryParameter < m_NumParams; entryParameter++) {
            if (fabs(gradient[entryParameter]) <= NEARLY_ZERO) { retryIndices.push_back(entryParameter); }
        }

        if (!retryIndices.empty()) {
            CalcDerivatives(retryIndices, point, objective, FD_RANGE_REL, gradient, fmin, xmin);
        }
    }

    m_GradCount++;

    return gradient;
}
