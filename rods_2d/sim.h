#ifndef SIM_H
#define SIM_H

/*
	Simulation of our model
*/

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "mt19937ar.h"

class rods {

	private:

		double beta, d;

		double *prob;

		int N, L, numPos, totalx, totaly, energy;

		double xLen, yLen;

		int *pos;

		char *lat;

		void sweep();

		void showNeighbors(int pos);

		double calcAvxLen();

		double calcAvyLen();

	public:

		rods(int size, double betaVal, double density, char *lattice = NULL);

		~rods();

		void run();

		char *getLat() { return lat; }

		int getEnergy();

		int recalcEnergy();

		double getMRatio();

		double getAvxLen();

		double getAvyLen();

		double getRhox() { return xLen; }

		double getRhoy() { return yLen; }

		double calcXLen();

		double calcYLen();

		double getAvLen();

		int *getXDist();

		int *getYDist();

		int getPredDirection();
};

#endif