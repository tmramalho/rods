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

		int L, N, A, numPos, totalx, totaly, totalz, energy;

		double xLen, yLen, zLen;

		int *pos;

		char *lat;

		void sweep();

		void showNeighbors(int pos);

		char randDirection();

		double calcXLen();

		double calcYLen();

		double calcZLen();

	public:

		rods(int size, double betaVal, double density);

		~rods();

		void run();

		char *getLat() { return lat; }

		int getEnergy() { return energy; }

		int recalcEnergy();

		double *getMRatio();

		double getAvxLen();

		double getAvyLen();

		double getAvzLen();

		double getRhox() { return xLen; }

		double getRhoy() { return yLen; }

		double getRhoz() { return zLen; }

		double getAvLen();

		int *getXDist();

		int *getYDist();

		int getPredDirection();
};

#endif