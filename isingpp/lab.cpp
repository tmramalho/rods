#include "lab.h"

lab::lab() {

}

lab::~lab() {

}

void lab::setTempRange(double _Ti, double _Tf, double _dT) {
	Ti = _Ti;
	Tf = _Tf;
	dT = _dT;
}

void lab::setDims(int _size, double _dens) {
	size = _size;
	dens = _dens;
}

void lab::simSimple(int runtime) {
	int *eArray = new int[runtime];
	int *mArray = new int[runtime];

	for(double T = Ti; T <= Tf; T += dT) {
	
		curSim = new rods(size, 1/T, dens);
		
		std::cout << std::endl << "Running simulation at " << T << std::endl;
		std::cout << "Initial energy: " << curSim->getEnergy() << std::endl;

		for(int i=0; i < runtime; i++) {
			curSim->run();
			eArray[i] = curSim->getEnergy();
			mArray[i] = curSim->getMag();
		}

		std::cout << "Final energy: " << curSim->getEnergy() << std::endl;
		std::cout << "Final mag: " << curSim->getMag() << std::endl;

		//write the final state to a file
		//saveState(T);
		saveParameter("energy", T, eArray, runtime);
		saveParameter("mag", T, mArray, runtime);

		delete curSim;
	}
	
	delete eArray;

}

void lab::simCorr(int runtime) {

	int *eArray = new int[runtime];
	int nmax = 5;
	double *cArray = new double[nmax];

	for(double T = Ti; T <= Tf; T += dT) {
		curSim = new rods(size, (double)1/T, dens);

		std::cout << std::endl << "Running simulation at " << T << std::endl;
		std::cout << "Initial energy: " << curSim->getEnergy() << std::endl;

		std::cout << "I had to wait: " << waitForEquilibrium() << std::endl;

		std::cout << "Final energy: " << curSim->getEnergy() << std::endl;

		for(int n = 0; n < nmax; n++) {
			std::cout << std::endl << "Measuring energy" << std::endl;

			for(int i = 0; i < runtime; i++) {
				curSim->run();
				eArray[i] = curSim->getEnergy();
			}

			#ifdef VERBOSE
			saveParameter("energy", T, eArray, runtime);
			#endif

			std::cout << "Calculating correlation time: " << std::endl;
			//calcAcFunction(eArray, runtime, 0);
			cArray[n] = calcAcFunction(eArray, runtime, 1);
		}
		saveCorrTime(T, calcAverage(cArray, nmax));

		delete curSim;
	}

	delete[] eArray;
	delete[] cArray;
}

void lab::simMeas() {

	int numSamples = 200;
	//alloc arrays
	int *eArray = new int[numSamples];
	double *m1Array = new double[numSamples];
	double *m2Array = new double[numSamples];
	double *lArray = new double[numSamples];
	double *ratioVec;
	len *xLen, *yLen;
	std::fstream filestr;

	for(double T = Ti; T <= Tf; T += dT) {
		curSim = new rods(size, (double)1/T, dens);
		xLen = new len(size);
		yLen = new len(size);
		
		int n = 0;
		std::cout << "Looking up correlation time for " << T << "..." << std::endl;
		int ctime = (int) floor(findCorrTime(T));

		if(ctime == 0) {
			std::cout << "Failed! Using ";
			ctime = 100;
		} else {
			std::cout << "Got it! Using ";
		}

		std::cout << "Waiting for equilibrium..." << std::endl;
		std::cout << "I had to wait: " << waitForEquilibrium() << std::endl;
		
		/*********************************************
		* Measure parameters every 2 correlation times
		*********************************************/

		int runtime = numSamples * 2 * ctime;
		std::cout << "Simulating..." << std::endl;
		std::cout << ctime << " -> " << runtime << std::endl;

		for(int i = 0; i < runtime; i++) {
			curSim->run();
			if(i % (2*ctime) == 0) {

				eArray[n] = curSim->getEnergy();

				
				n++;
				std::cout << ".";
				
				delete ratioVec;
			}
		}
		std::cout << std::endl;
		
		/**********************************
		* Calculate averages of all parameters
		**********************************/
		
		double avE = calcAverage(eArray, numSamples);
		double avsqE = calcSqAverage(eArray, numSamples);
		double numParticles = (double) size * size * size * dens;
		double c = (avsqE - avE*avE) / ( numParticles * T * T);
		double avM1 = calcAverage(m1Array, numSamples);
		double avsqM1 = calcSqAverage(m1Array, numSamples);
		double avM2 = calcAverage(m2Array, numSamples);
		double avsqM2 = calcSqAverage(m2Array, numSamples);
		double avl = calcAverage(lArray, numSamples);
		double avsql = calcSqAverage(lArray, numSamples);

		/**********************************
		* Write averages of all parameters
		**********************************/

		filestr.open("master.dat", std::fstream::in | std::fstream::out | std::fstream::app);
		filestr.seekg (0, std::fstream::end);
		if(filestr.tellg() == 0) {
			//quick check to see if master.dat is empty -> write header
			filestr << "T\trho\t<E>\t<E^2>\tc\t<M>\t<M^2>\t<l>\t<l^2>\t";
			filestr << "<l_dominant>\t<l_minority>\t<\\Delta>\t<\\Delta^2>" << std::endl;
		}

		filestr << T << "\t" << dens << "\t" << avE << "\t" << avsqE << "\t";
		filestr << c << "\t" << avM1 << "\t" << avsqM1 << "\t";
		filestr << avM2 << "\t" << avsqM2 << "\t" << avl << "\t" << avsql << std::endl;

		filestr.close();

		/**********************************
		* Write average length distribution
		**********************************/
		//TODO:not implemented yet

		#ifdef VERBOSE
		saveState(T);
		#endif

		delete curSim;
		delete xLen;
		delete yLen;
	}

	delete[] eArray;
	delete[] m1Array;
	delete[] m2Array;
	delete[] lArray;
}

int lab::waitForEquilibrium() {
	bool eq = false;
	int counter = 0;

	while(!eq) {
		int total = 0;
		int prevEnergy = curSim->getEnergy();
		int curEnergy;
		for(int i=0; i < 10000; i++) {
			curSim->run();
			curEnergy = curSim->getEnergy();
			total += prevEnergy - curEnergy;
			prevEnergy = curEnergy;
		}
		counter += 10000;
		std::cout << "Slope( " << counter << "):" << total << std::endl;
		if(counter > 20000) eq = true;
	}

	return counter;
}

void lab::saveState(double T) {
	/*std::stringstream oss1;
	std::stringstream oss2;
	std::stringstream oss3;
	std::string filename;
	oss1 << "lattice - T" << T << "rho" << dens << "x.dat";
	FILE *fp1 = fopen(oss1.str().c_str(), "w");
	oss2 << "lattice - T" << T << "rho" << dens << "y.dat";
	FILE *fp2 = fopen(oss2.str().c_str(), "w");
	oss3 << "lattice - T" << T << "rho" << dens << "z.dat";
	FILE *fp3 = fopen(oss3.str().c_str(), "w");

	short *s = curSim->getLat();

	for(int i=0; i < size; i++) {
		for(int j=0; j < size; j++) {
			for(int k=0; k < size; k++) {
				int pos = i + size*j + size*size*k;
				if(s[pos] == 1) fprintf(fp1, "%d %d %d\n", i, j, k);
				else if(s[pos] == 2) fprintf(fp2, "%d %d %d\n", i, j, k);
				else if(s[pos] == 4) fprintf(fp3, "%d %d %d\n", i, j, k);
			}
		}
	}

	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
*/
}

void lab::saveParameter(const char* name, double temp, int *array, int size) {
	std::stringstream oss;
	std::string filename;
	oss << name << temp << ".dat";
	//filename = oss.str();

	FILE *fp = fopen(oss.str().c_str(), "w");

	for(int i=0; i < size; i++) {
		fprintf(fp, "%d\n", array[i]);
	}

	fclose(fp);
}

void lab::saveParameter(const char* name, double temp, float *array, int size) {
	std::stringstream oss;
	std::string filename;
	oss << name << temp << ".dat";
	//filename = oss.str();

	FILE *fp = fopen(oss.str().c_str(), "w");

	for(int i=0; i < size; i++) {
		fprintf(fp, "%f\n", array[i]);
	}

	fclose(fp);
}

void lab::saveParameter(const char* name, double temp, double *array, int size) {
	std::stringstream oss;
	std::string filename;
	oss << name << temp << ".dat";
	//filename = oss.str();

	FILE *fp = fopen(oss.str().c_str(), "w");

	for(int i=0; i < size; i++) {
		fprintf(fp, "%lf\n", array[i]);
	}

	fclose(fp);
}

double lab::calcAcFunction(int *vArray, int dt, int mode) {

	int tend = floor(0.9*(double)dt);
	double sum[3];
	double xsi[tend];

	if( mode == 0 ) {
		/***************
		3.21 Newman
		***************/
		for(int t = 0; t < tend; t++) {
			int tmax = dt-t;
			double ct = 1/((double)tmax);
			sum[0] = sum[1] = sum[2] = 0;
			for(int k = 0; k < tmax; k++) {
				sum[0] += vArray[k]*vArray[k+t];
				sum[1] += vArray[k];
				sum[2] += vArray[k+t];
			}
			xsi[t] = ct*sum[0]-ct*ct*sum[1]*sum[2];
		}
	} else {
		/***************
		3.17 Newman
		***************/
		double maverage = 0;
		for(int i = 0; i < dt; i++) {
			maverage += vArray[i];
		}
		maverage /= dt;

		for(int t = 0; t < tend; t++) {
			xsi[t] = 0;
			for(int k = 0; k < dt-t; k++) {
				xsi[t] += (vArray[k] - maverage) * (vArray[k+t] - maverage);
			}
		}
	}

	#ifdef VERBOSE
	saveAcFunction(mode, xsi, tend);
	#endif

	double result = integrateAcFunction(xsi, tend);

	std::cout << "Correlation time: " << result << std::endl;

	return result;

}

double lab::integrateAcFunction(double *xsi, int tend) {

	double tau = 0;

	for(int t = 0; t < tend; t++) {
		//3.20 Newman
		tau += (xsi[t]/xsi[0]);
		if(xsi[t]/xsi[0] < 0.2) break;
	}

	return tau;

}

void lab::saveAcFunction(int mode, double *xsi, int tend) {
	std::stringstream oss;
	std::string filename;
	oss << "ac" << mode << ".dat";
	filename = oss.str();
	FILE *fp = fopen(filename.c_str(), "w");

	for(int t = 0; t < tend; t++) {
		//if(xsi[t] < 0) break;
		fprintf(fp, "%f\n", xsi[t]/xsi[0]);
	}

	fclose(fp);

}

void lab::saveCorrTime(double T, double res) {
	FILE *fp = fopen("ctimes.dat", "a");

	fprintf(fp, "%f %f\n", T, res);

	fclose(fp);

}

double lab::findCorrTime(double T) {
	double Tr, ctime = 0, res = 0;

	//read ctimes file
	FILE *fp = fopen("ctimes.dat", "r");
	if(fp == NULL) {
		std::cout << "ctimes .dat does not exist!" << std::endl;
		return 0;
	}

	while(fscanf(fp, "%lf %lf", &Tr, &res) != EOF) {
		if(floor(Tr*100) == floor(T*100)) {
			ctime = res;
		}
	}

	fclose(fp);

	return ctime;
}

double lab::calcAverage(int *vArray, int size) {
	double result = 0;

	for(int i = 0; i < size; i++) {
		result += (double) vArray[i];
	}

	return (result / size);
}

double lab::calcAverage(double *vArray, int size) {
	double result = 0;

	for(int i = 0; i < size; i++) {
		result += vArray[i];
	}

	return (result / size);
}

double lab::calcSqAverage(int *vArray, int size) {
	double result = 0;

	for(int i = 0; i < size; i++) {
		result += (double) vArray[i]*vArray[i];
	}

	return (result / size);
}

double lab::calcSqAverage(double *vArray, int size) {
	double result = 0;

	for(int i = 0; i < size; i++) {
		result += vArray[i]*vArray[i];
	}

	return (result / size);
}
