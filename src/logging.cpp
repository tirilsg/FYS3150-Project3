#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include "logging.h"

void saveDataToTxt(const std::string& filename, const arma::mat& data) { //a simple function that takes a data set and name of a file as arguments
    std::ofstream file(filename.c_str()); //opens file
    if (file.is_open()) { //stores data in file, and closes it
        file << data;
        file.close(); 
    }
}

