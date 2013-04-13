#ifndef LAB_H
#define LAB_H

/*
	Laboratory class
	This class will manage the sim simulation
*/

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <sys/stat.h>

#include "mt19937ar.h"
#include "sim.h"
#include "len.h"

#define VERBOSE 1

class lab {

	private:

		double Ti, Tf, dT, dens;

		int size;

		rods *curSim;

		void saveState(double T);

		void saveParameter(const char *name, double temp, int *array, int size);
		void saveParameter(const char *name, double temp, float *array, int size);
		void saveParameter(const char *name, double temp, double *array, int size);

		int waitForEquilibrium();

		double integrateAcFunction(double *xsi, int tend);

		void saveAcFunction(int mode, double *xsi, int tend);

		double calcAcFunction(int *vArray, int dt, int mode = 0);

		void saveCorrTime(double T, double res);

		double findCorrTime(double T);

		double calcAverage(int *vArray, int size);
		double calcAverage(double *vArray, int size);

		double calcSqAverage(int *vArray, int size);
		double calcSqAverage(double *vArray, int size);

		int *sumDists(int *dist1, int *dist2);

		struct stat fileinfo;

	public:

		lab();

		void setTempRange(double Ti, double Tf, double dT);

		void setDims(int size, double dens);

		void simSimple(int runtime, int startMode);

		void simCorr(int runtime);

		void simMeas();

		~lab();

};

#endif