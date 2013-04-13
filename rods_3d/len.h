#ifndef LEN_H
#define LEN_H

/*
	Lendist
	This class manages the distribution of array lengths
*/

#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>

#include "mt19937ar.h"
#include "sim.h"

class len {

	private:
		int *master;
		int size;
		int numSamples;

	public:

		len(int maxLengthSize);

		~len();

		void addNewMeasurement(int *dist);

		double getAvLen(int pos);

};

#endif