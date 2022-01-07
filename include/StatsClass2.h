
#ifndef STATS_CLASS_H
#define STATS_CLASS_H

// Include standard classes
#include <vector>

// Include custom classes
#include "MyHeaderInc.h"
#include "ParameterGroup.h"
#include <ParameterABC.h>


void CalculateJacobianParameters(std::vector<double> currentValues, std::vector<double> lowerValues, std::vector<double> upperValues, 
                                 double stepSize, double tol, std::vector<std::vector<double>> &jacobianLocations, 
                                 std::vector<bool> &lockedParameter);

void CalculateJacobian2(std::vector<double> currentValues, double currentObjective, std::vector<double> objectives, std::vector<std::vector<double>> jacobianLocations, 
                       std::vector<bool> lockedParameter, std::vector<std::vector<double>>& jacobian);

std::vector<std::vector<double>> Tranpose(std::vector<std::vector<double>> inputMatrix);



#endif /* STATS_CLASS_H */