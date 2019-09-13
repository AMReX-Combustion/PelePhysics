#include "chemistry_file.H"
#include <sstream>
#include <iostream>
#include <fstream>

/*Poly fits for the diffusion coefficients, dim NO*KK*KK */
void egtransetCOFD(double* COFD) {
    
    std::ifstream filetoread;
    filetoread.open("COFD.dat");

    char idx[10];
    char diffcoefs[25];

    for (int i = 0; i < 99856 ; i++) {
	    filetoread>>idx;
	    filetoread>>diffcoefs;
	    //std::cout<< idx << " "<< diffcoefs << std::endl;
	    COFD[std::stoi(idx)] = atof(diffcoefs);
	    filetoread.ignore(80,'\n');
    }

}
