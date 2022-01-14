
#include "MathClass2.h"

void FirstOrderBackwardDifference(double currentValue, double lowerValue, double stepSize, double &lowerAdjusted) {
    // todo: doc string

    // Calculate the adjustment from the current value
    double adjustment = currentValue * stepSize;

    // Calculate the lower bound value
    lowerAdjusted = currentValue - adjustment;

    if (lowerAdjusted < lowerValue) {
        lowerAdjusted = lowerValue;
    }
    
}

void FirstOrderForwardDifference(double currentValue, double upperValue, double stepSize, double &upperAdjusted) {
    // todo: doc string

    // Calculate the adjustment from the current value
    double adjustment = currentValue * stepSize;

    // Calculate the upper bound value
    upperAdjusted = currentValue + adjustment;

    if (upperAdjusted > upperValue) {
        upperAdjusted = upperValue;
    }

}

void FirstOrderCentralDifference(double currentValue, double lowerValue, double upperValue, double stepSize, double& lowerAdjusted, double& upperAdjusted) {
    // todo: doc string

    FirstOrderBackwardDifference(currentValue, lowerValue, stepSize, lowerAdjusted);
    FirstOrderForwardDifference(currentValue, upperValue, stepSize, upperAdjusted);
}

void EvaluateTolerance(double currentValue, double upperValue, double lowerValue, double &percentRange) {
    // todo: doc string

    percentRange = std::abs((upperValue - lowerValue) / currentValue);
}


std::vector<double> CreateParameterSet(std::vector<double> currentValues, double adjustedValue, int index) {
    // todo: doc string

    // Copy the current parameter vector
    std::vector<double> updatedParameterSet = currentValues;

    // Set the new parameter value into the vector
    updatedParameterSet[index] = adjustedValue;

    // Return to the calling function
    return updatedParameterSet;

}


std::vector<std::vector<double>> Tranpose(std::vector<std::vector<double>> inputMatrix) {
    // todo: doc string

    std::vector<std::vector<double>> outputMatrix(inputMatrix[0].size(), std::vector<double>(inputMatrix.size()));
    for (size_t i = 0; i < inputMatrix.size(); ++i)
        for (size_t j = 0; j < inputMatrix[0].size(); ++j)
            outputMatrix[j][i] = inputMatrix[i][j];

    return outputMatrix;

}

void CalculateJacobianParameters(std::vector<double> currentValues, std::vector<double> lowerValues, std::vector<double> upperValues, double stepSize, double tol, int order,
                                 std::vector<std::vector<double>> &jacobianLocations, std::vector<bool> &lockedParameter) {
    // todo: doc string

    // Loop on each parameter and perturb it
    if (order == FIRST_BACKWARD) {
        // First order backward difference
        for (int entryParameter = 0; entryParameter < currentValues.size(); entryParameter++) {
            double lowerAdjusted, percentRange;

            // Find the backward location
            FirstOrderBackwardDifference(currentValues[entryParameter], lowerValues[entryParameter], stepSize, lowerAdjusted);

            // Test against the tolerance
            EvaluateTolerance(currentValues[entryParameter], currentValues[entryParameter], lowerAdjusted, percentRange);

            // Determine action based on whether the tolerance is met
            if (percentRange >= tol){
                // Parameter is active. Set into the location vector.
                jacobianLocations.push_back(CreateParameterSet(currentValues, lowerAdjusted, entryParameter));
                lockedParameter.push_back(false);

            } else {
                // Lock the parameter
                lockedParameter.push_back(true);
            }
        }

    } else if (order == FIRST_CENTRAL) {
        // First order central difference
        for (int entryParameter = 0; entryParameter < currentValues.size(); entryParameter++) {
            double lowerAdjusted, upperAdjusted, percentRange;

            // Find the backward and forward locations
            FirstOrderBackwardDifference(currentValues[entryParameter], lowerValues[entryParameter], stepSize, lowerAdjusted);
            FirstOrderForwardDifference(currentValues[entryParameter], upperValues[entryParameter], stepSize, upperAdjusted);

            // Calculate parameter range as a fraction of the current value
            EvaluateTolerance(currentValues[entryParameter], upperAdjusted, lowerAdjusted, percentRange);

            // If the parameter range is larger than the tolerance, compute the parameter. Otherwise lock the parameter from the calcuation
            if (percentRange >= tol) {
                // Construct the parameter sets
                jacobianLocations.push_back(CreateParameterSet(currentValues, lowerAdjusted, entryParameter));
                jacobianLocations.push_back(CreateParameterSet(currentValues, upperAdjusted, entryParameter));

                // Indicate that the parameter is unlocked
                lockedParameter.push_back(false);

            } else {
                // Indicate that the parameter is locked
                lockedParameter.push_back(true);
            }
        }

    } else if (order == FIRST_FORWARD) {
        // First order forward difference
        for (int entryParameter = 0; entryParameter < currentValues.size(); entryParameter++) {
            double upperAdjusted, percentRange;

            // Find the forward location
            FirstOrderForwardDifference(currentValues[entryParameter], upperValues[entryParameter], stepSize, upperAdjusted);

            // Test against the tolerance
            EvaluateTolerance(currentValues[entryParameter], currentValues[entryParameter], upperAdjusted, percentRange);

            // Determine action based on whether the tolerance is met
            if (percentRange >= tol) {
                // Parameter is active. Set into the location vector.
                jacobianLocations.push_back(CreateParameterSet(currentValues, upperAdjusted, entryParameter));
                lockedParameter.push_back(false);

            } else {
                // Lock the parameter
                lockedParameter.push_back(true);
            }
        }
    }

    // todo: add second order differences
}

void CalculateJacobian2(std::vector<double> currentValues, double currentObjective, std::vector<double> objectives, std::vector<std::vector<double>> jacobianLocations, int order, 
                        std::vector<bool> &lockedParameter, std::vector<std::vector<double>> &jacobian) {
    
    int indexPosition = 0;

    if (order == FIRST_BACKWARD) {
        // First order backward difference
        std::vector<double> functionDerivatives;

        for (int entryParameter = 0; entryParameter < lockedParameter.size(); entryParameter++) {
            if (!lockedParameter[entryParameter]) {
                // Get the partial derivative of the parameter. Order is swapped to implicitly introduce the negative direction.
                double partialDerivatives = (objectives[indexPosition] - currentObjective) / (currentValues[entryParameter] - jacobianLocations[indexPosition][entryParameter]);

                // Add to the derivative vector
                functionDerivatives.push_back(partialDerivatives);

                // If gradietn is zero, lock the parameter
                if (partialDerivatives == 0) {
                    lockedParameter[entryParameter] = true;
                }

                // Increment the index counter 
                indexPosition++;

            } else {
                // Parameter is locked. Add a zero into the derivative vector as a placeholder
                functionDerivatives.push_back(0);
            }
        }

        // Append the function derivatives into the jacobian
        jacobian.push_back(functionDerivatives);


    } else if (order == FIRST_CENTRAL) {
        // First order central difference
        std::vector<double> functionDerivatives;

        for (int entryParameter = 0; entryParameter < lockedParameter.size(); entryParameter++) {
            if (!lockedParameter[entryParameter]) {
                // Get the partial derivative of the parameter
                double partialDerivatives = (objectives[indexPosition * 2 + 1] - objectives[indexPosition * 2] ) / 
                                            (jacobianLocations[indexPosition * 2 + 1][entryParameter] - jacobianLocations[indexPosition * 2][entryParameter]);

                // Add to the derivative vector
                functionDerivatives.push_back(partialDerivatives);

                // If gradietn is zero, lock the parameter
                if (partialDerivatives == 0) {
                    lockedParameter[entryParameter] = true;
                }

                // Increment the index counter 
                indexPosition++;

            } else {
                // Parameter is locked. Add a zero into the derivative vector as a placeholder
                functionDerivatives.push_back(0);
            }
        }

        // Append the function derivatives into the jacobian
        jacobian.push_back(functionDerivatives);

    } else if (order == FIRST_FORWARD) {
        // First order forward difference
        std::vector<double> functionDerivatives;

        for (int entryParameter = 0; entryParameter < lockedParameter.size(); entryParameter++) {
            if (!lockedParameter[entryParameter]) {
                // Get the partial derivative of the parameter
                double partialDerivatives = (objectives[indexPosition] - currentObjective) / (jacobianLocations[indexPosition][entryParameter] - currentValues[entryParameter]);

                // Add to the derivative vector
                functionDerivatives.push_back(partialDerivatives);

                // If gradietn is zero, lock the parameter
                if (partialDerivatives == 0) {
                    lockedParameter[entryParameter] = true;
                }

                // Increment the index counter 
                indexPosition++;

            } else {
                // Parameter is locked. Add a zero into the derivative vector as a placeholder
                functionDerivatives.push_back(0);
            }
        }

        // Append the function derivatives into the jacobian
        jacobian.push_back(functionDerivatives);

    }
}