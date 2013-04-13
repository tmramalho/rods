/*
	Monte carlo 3D rods simulation
	2008-2009 Tiago Ramalho
*/

#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>

#include "sim.h"
#include "lab.h"
#include "mt19937ar.h"

void get_parameters(int argc, char *argv[], lab *mylab);

int main(int argc, char *argv[]) {
	init_genrand(time(NULL));

	lab *laboratory = new lab();

	get_parameters(argc, argv, laboratory);

	return 0;
}

void get_parameters(int argc, char *argv[], lab *mylab) {

	int size = 100;
	int dtime = 10000;
	double dT = 0.01;
	double rho = 0.1;
	double Tf = 0.2;
	double Ti = 0.1;
	int numSamples = 200;
	int mode = 0;

	std::cout << "                          _|            \n_|  _|_|    _|_|      _|_|_|    _|_|_|  \n_|_|      _|    _|  _|    _|  _|_|      \n_|        _|    _|  _|    _|      _|_|  \n_|          _|_|      _|_|_|  _|_|_|    \n" << std::endl;
	std::cout << "Usage : " << argv[0] << " mode Ti rho Tf dT dt L numSamples[1]" << std::endl;

	switch(argc) {
		case 9:
			numSamples = atoi(argv[8]);
		case 8:
			size = atoi(argv[7]);
		case 7:
			dtime = atoi(argv[6]);
		case 6:
			dT = atof(argv[5]);
		case 5:
			Tf = atof(argv[4]);
		case 4:
			rho = atof(argv[3]);
		case 3:
			Ti = atof(argv[2]);
		case 2:
			mode = atoi(argv[1]);
		default:
			break;
	}

	mylab->setTempRange(Ti, Tf, dT);
	mylab->setDims(size, rho);

	std::cout << "using params: Ti: " << Ti;
	std::cout << " Tf: " << Tf << " dT: " << dT;
	std::cout << " rho: " << rho << " dt: " << dtime;
	std::cout << " L: " << size << std::endl;

	switch(mode) {
		case 0:
			mylab->simSimple(dtime);
			break;
		case 1:
			mylab->simCorr(dtime);
			break;
		case 2:
			mylab->simMeas();
			break;
	}
}

