#ifndef SIM_H
#define SIM_H

/*
	Fake ising model
	Using the nicer c++ implementation
	of rods to run the ising model
	full of hacks
	not production quality code!
*/

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "mt19937ar.h"

#define XNN 1
#define YNN L

class rods {

	private:

		double beta, d;

		double *prob;

		int L, N, A, numPos, totalx, totaly, totalz, energy, mag;

		double xLen, yLen, zLen;

		int *pos;

		short *lat;

		void sweep();

		void showNeighbors(int pos);

		char randDirection();

	public:

		rods(int size, double betaVal, double density);

		~rods();

		void run();

		short *getLat() { return lat; }

		int getEnergy() { return energy; }

		int getMag() { return mag; }

		int recalcEnergy();

};

#endif
