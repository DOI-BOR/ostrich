#ifndef READ_UTILITY_H
#define READ_UTILITY_H


#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <filesystem>



std::vector<std::vector<std::vector<double>>> ReadWorkerFiles(int numberOfWorkers, bool solveOnPrimary);
std::vector<std::vector<double>> ReadOutputFile();

#endif 