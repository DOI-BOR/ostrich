
#include "ReadUtility.h"

std::vector<double> ParseLineToDoubles(std::string line) {

	// Create a vector to hold the converted values
	std::vector<double> parsedValues;

	// Create the offset value for the line
	std::string::size_type offset_line = 0;

	// Remove trailling whitespaces
	while (line.back() == ' ') {
		line.pop_back();
	}

	// Parse the values from the line
	while (offset_line < line.length()) {

		// Create the offset value for the substring
		std::string::size_type offset_substring = 0;

		// Read a line from the file
		double value = std::stod(line.substr(offset_line), &offset_substring);

		// Append the line into the holder
		parsedValues.push_back(value);

		// Add the subtring offset to the line offset to shift the line
		offset_line += offset_substring;

	}

	// Return to the calling function
	return parsedValues;

}


std::vector<std::vector<std::vector<double>>> ReadWorkerFiles(int numberOfWorkers, bool solveOnPrimary) {

	// Set the starting worker number
	int startingWorker = 0;

	// Define the vector to hold the data
	std::vector<std::vector<std::vector<double>>> previousSolves;

	if (!solveOnPrimary) {
		// Increment the starting counter
		startingWorker = 1;

		// Push back an empty value
		std::vector<std::vector<double>> temp;
		previousSolves.push_back(temp);
	}

	// Iterate on the worker files
	for (int entryWorker = startingWorker; entryWorker < numberOfWorkers; entryWorker++) {
		// Create the name of the file
		// TODO: This should be combined into a function from ModelWorker::Write if it works
		std::filesystem::path filePath = std::filesystem::current_path();
		filePath /= "OstModel" + std::to_string(entryWorker) + ".txt";
		filePath = filePath.lexically_normal();

		// Create the vector to hold the data
		std::vector<std::vector<double>> workerSolves;

		// Open the file
		std::fstream workerFile;
		workerFile.open(filePath, std::ios::in);
		
		// Read each line and extract information
		int lineCounter = 0;
		std::string line;
		while (std::getline(workerFile, line)) {
			if (lineCounter > 0) {				
				// Append the line into the vector for the file
				workerSolves.push_back(ParseLineToDoubles(line));

			}

			// Increment the line counter
			lineCounter++;

		}

		// Close the worker file
		workerFile.close();

		// Append the worker information into the final data vector
		previousSolves.push_back(workerSolves);

	}

	// Return the worker solves to the calling function
	return previousSolves;

}


std::vector<std::vector<double>> ReadOutputFile() {

	// Create the name of the file
	// TODO: This should be combined into a function from ModelWorker::Write if it works
	std::filesystem::path filePath = std::filesystem::current_path();
	filePath /= "OstOutput0.txt";
	filePath = filePath.lexically_normal();

	// Create the holder for the iteration information
	std::vector<std::vector<double>> iterationValues;

	// Open the file
	std::fstream workerFile;
	workerFile.open(filePath, std::ios::in);

	// Read each line and extract the iteration information
	int lineCounter = 0;
	std::string line;

	bool parseLine = false;
	while (std::getline(workerFile, line)) {
		// Check for when the iteration block begins
		if (line.substr(0, 3) == "Run") {
			parseLine = true;
			continue;
		}

		if (parseLine && line.length() > 2) {
			// Parse the values from the line
			iterationValues.push_back(ParseLineToDoubles(line));

		} else if (parseLine) {
			// Disable the parse and continue to the next steps
			parseLine = false;
			break;

		}

		// Increment the line counter 
		lineCounter++;
	}

	// Return the iteration values to the calling function
	return iterationValues;

		

}