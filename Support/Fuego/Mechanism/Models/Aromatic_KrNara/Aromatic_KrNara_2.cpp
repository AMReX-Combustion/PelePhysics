#include "chemistry_file.H"
#include <sstream>
#include <iostream>
#include <fstream>

/*Poly fits for the viscosities, dim NO*KK */
void egtransetCOFETA(double* COFETA) {

    std::ifstream filetoread;
    filetoread.open("COFETA.dat");

    char diffcoefs[25];

    for (int i = 0 ; i < 632 ; i++) {
	    filetoread>>diffcoefs;
	    COFETA[i] = atof(diffcoefs);
	    filetoread.ignore(80,'\n');
    }
}


/*Poly fits for the conductivities, dim NO*KK */
void egtransetCOFLAM(double* COFLAM) {

    std::ifstream filetoread;
    filetoread.open("COFLAM.dat");

    char diffcoefs[25];

    for (int i = 0 ; i < 632 ; i++) {
	    filetoread>>diffcoefs;
	    COFLAM[i] = atof(diffcoefs);
	    filetoread.ignore(80,'\n');
    }
}


/*List of specs with small weight, dim NLITE */
void egtransetKTDIF(int* KTDIF) {
    KTDIF[0] = 1;
    KTDIF[1] = 5;
}


/*Poly fits for thermal diff ratios, dim NO*NLITE*KK */
void egtransetCOFTD(double* COFTD) {

    std::ifstream filetoread;
    filetoread.open("COFTD.dat");

    char diffcoefs[25];

    for (int i = 0 ; i < 1264 ; i++) {
	    filetoread>>diffcoefs;
	    COFTD[i] =atof(diffcoefs);
	    filetoread.ignore(80,'\n');
    }
}

