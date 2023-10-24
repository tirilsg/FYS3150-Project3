#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>

#ifndef _logging_h__  
#define __logging_h__

void saveDataToTxt(const std::string& filename, const arma::mat& data);


#endif