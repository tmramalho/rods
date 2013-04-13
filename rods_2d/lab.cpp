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

void lab::simSimple(int runtime, int startMode) {
	int *eArray = new int[runtime];
	char *prevLattice = NULL;

	for(double T = Ti; T <= Tf; T += dT) {
		curSim = new rods(size, 1/T, dens, prevLattice);
		std::cout << std::endl << "Running simulation at " << T << std::endl;
		std::cout << "Initial energy: " << curSim->getEnergy() << std::endl;

		for(int i=0; i < runtime; i++) {
			curSim->run();
			eArray[i] = curSim->getEnergy();
		}

		std::cout << "Final energy: " << curSim->getEnergy() << std::endl;
		std::cout << "M ratio: " << curSim->getMRatio() << std::endl;
		std::cout << "Av length: " << curSim->getAvLen() << std::endl;
		std::cout << "Av x length: " << curSim->getAvxLen() << std::endl;
		std::cout << "Av y length: " << curSim->getAvyLen() << std::endl;
		std::cout << "rho_x: " << curSim->getRhox() << std::endl;
		std::cout << "rho_y: " << curSim->getRhoy() << std::endl;

		#ifdef VERBOSE
		//write the final state to a file
		saveState(T);
		saveParameter("energy", T, eArray, runtime);
		#endif

		if(startMode == 1) {
			char *temp = curSim->getLat();
			prevLattice = new char[(size*size)];
			for(int n = 0; n < (size * size); n++) {
				prevLattice[n] = temp[n];
			}
		}

		delete curSim;
	}

	if(prevLattice != NULL) delete[] prevLattice;
	delete eArray;

}

void lab::simCorr(int runtime) {

	int *eArray = new int[runtime];

	for(double T = Ti; T <= Tf; T += dT) {
		curSim = new rods(size, (double)1/T, dens);

		std::cout << std::endl << "Running simulation at " << T << std::endl;
		std::cout << "Initial energy: " << curSim->getEnergy() << std::endl;

		std::cout << "I had to wait: " << waitForEquilibrium() << std::endl;

		std::cout << "Final energy: " << curSim->getEnergy() << " vs ";
		std::cout << curSim->recalcEnergy() << std::endl;

		std::cout << std::endl << "Measuring energy..." << std::endl;

		for(int i = 0; i < runtime; i++) {
			curSim->run();
			eArray[i] = curSim->getEnergy();
		}

		saveParameter("energy", T, eArray, runtime);

		std::cout << std::endl << "Calculating acfunction..." << std::endl;
		//calcAcFunction(eArray, runtime, 0);
		double res = calcAcFunction(eArray, runtime, 1);
		saveCorrTime(T, res);

		delete curSim;
	}

	delete[] eArray;
}

void lab::simMeas() {

	int numSamples = 200;
	//alloc arrays
	int *eArray = new int[numSamples];
	double *mrArray = new double[numSamples];
	double *lArray = new double[numSamples];
	double *ldArray = new double[numSamples];
	double *lmArray = new double[numSamples];
	len *xLen, *yLen, *totLen;
	std::fstream filestr;

	for(double T = Ti; T <= Tf; T += dT) {
		curSim = new rods(size, (double)1/T, dens);
		xLen = new len(size);
		yLen = new len(size);
		totLen = new len(size);
		int n = 0;
		std::cout << "Looking up correlation time for " << T << "..." << std::endl;
		int ctime = (int) floor(findCorrTime(T));

		if(ctime == 0) {
			std::cout << "Failed! Using ";
			ctime = 100;
		} else {
			std::cout << "Got it! Using ";
		}

		int runtime = numSamples * 2 * ctime;

		std::cout << ctime << " " << runtime << std::endl;
		std::cout << "Waiting for equilibrium..." << std::endl;
		std::cout << "I had to wait: " << waitForEquilibrium() << std::endl;
		std::cout << "Simulating..." << std::endl;

		for(int i = 0; i < runtime; i++) {
			curSim->run();
			if(i % (2*ctime) == 0) {
				eArray[n] = curSim->getEnergy();
				mrArray[n] = curSim->getMRatio();
				lArray[n] = curSim->getAvLen();
				if(curSim->getPredDirection() == 1) {
					ldArray[n] = curSim->getAvxLen();
					lmArray[n] = curSim->getAvyLen();
				} else if(curSim->getPredDirection() == 2) {
					ldArray[n] = curSim->getAvxLen();
					lmArray[n] = curSim->getAvyLen();
				}
				int *xdist = curSim->getXDist();
				int *ydist = curSim->getYDist();
				totLen->addNewMeasurement(sumDists(xdist, ydist));
				xLen->addNewMeasurement(curSim->getXDist());
				yLen->addNewMeasurement(curSim->getYDist());
				std::cout << "Measured! At " << i << std::endl;
				std::cout << eArray[n] << "; " << mrArray[n] << "; " << lArray[n] << std::endl;
				//if(mrArray[n] == 1) saveState(T+i);
				n++;
			}
		}

		double avE = calcAverage(eArray, numSamples);
		double avsqE = calcSqAverage(eArray, numSamples);
		double numParticles = (double) size * size * dens;
		double c = (avsqE - avE*avE) / ( numParticles * T * T);
		double avM = calcAverage(mrArray, numSamples);
		double avsqM = calcSqAverage(mrArray, numSamples);
		double avl = calcAverage(lArray, numSamples);
		double avsql = calcSqAverage(lArray, numSamples);
		double avld = calcAverage(ldArray, numSamples);
		double avlm = calcAverage(lmArray, numSamples);

		std::cout << "<Energy> = " << avE << std::endl;
		std::cout << "<Energy^2> = " << avsqE << std::endl;
		std::cout << "Specific heat = " << c << std::endl;
		std::cout << "<Mratio> = " << avM << std::endl;
		std::cout << "<l> = " << avl << std::endl;

		/**********************************
		* Write averages of all parameters
		**********************************/

		filestr.open("master.dat", std::fstream::out | std::fstream::app);

		if(stat("master.dat", &fileinfo) == 0) {
			//quick check to see if master.dat is empty -> write header
			//filestr << "T;rho;<E>;<E^2>;c;<M>;<M^2>;<l>;<l^2>;";
			//filestr << "<l_dominant>;<l_minority>" << std::endl;
		}

		filestr << T << ";" << dens << ";" << avE << ";" << avsqE << ";";
		filestr << c << ";" << avM << ";" << avsqM << ";" << avl << ";";
		filestr << avsql << ";" << avld << ";" << avlm << std::endl;
		filestr.close();

		/**********************************
		* Write average length distribution
		**********************************/
		filestr.open("dists.dat", std::fstream::out | std::fstream::app);

		if(stat("dists.dat", &fileinfo) == 0) {
			//quick check to see if master.dat is empty -> write header
			filestr << "T;<{l_x}_i>;<{l_y}_i>" << std::endl;
		}

		filestr << T << std::endl;
		for(int i = 0; i < size; i++)
			filestr << totLen->getAvLen(i) << ";";
		filestr << std::endl;
		for(int i = 0; i < size; i++)
			filestr << xLen->getAvLen(i) << ";";
		filestr << std::endl;
		for(int i = 0; i < size; i++)
			filestr << yLen->getAvLen(i) << ";";
		filestr << std::endl;

		filestr.close();

		#ifdef VERBOSE
		saveState(T);
		#endif

		delete curSim;
		delete xLen;
		delete yLen;
	}

	delete[] eArray;
	delete[] mrArray;
	delete[] lArray;
	delete[] ldArray;
	delete[] lmArray;
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
		if(counter > 200000 && abs(total) < 10) eq = true;
	}

	return counter;
}

void lab::saveState(double T) {
	std::stringstream oss;
	std::string filename;
	oss << "lattice - T" << T << "rho" << dens << ".dat";
	filename = oss.str();
	FILE *fp = fopen(filename.c_str(), "w");

	char *s = curSim->getLat();

	for(int i=0; i < ( size*size ); i++) {
		putc((char) s[i], fp);
		if((i+1) % size == 0) putc((char)255, fp);
	}

	fclose(fp);

}

void lab::saveParameter(const char* name, double temp, int *array, int size) {
	std::stringstream oss;
	std::string filename;
	oss << name << " - T" << temp << "rho" << dens << ".dat";
	filename = oss.str();

	FILE *fp = fopen(filename.c_str(), "w");

	for(int i=0; i < size; i++) {
		fprintf(fp, "%d\n", array[i]);
	}

	fclose(fp);
}

void lab::saveParameter(const char* name, double temp, float *array, int size) {
	std::stringstream oss;
	std::string filename;
	oss << name << " - T" << temp << "rho" << dens << ".dat";
	filename = oss.str();

	FILE *fp = fopen(filename.c_str(), "w");

	for(int i=0; i < size; i++) {
		fprintf(fp, "%f\n", array[i]);
	}

	fclose(fp);
}

void lab::saveParameter(const char* name, double temp, double *array, int size) {
	std::stringstream oss;
	std::string filename;
	oss << name << " - T" << temp << "rho" << dens << ".dat";
	filename = oss.str();

	FILE *fp = fopen(filename.c_str(), "w");

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
		if(xsi[t]/xsi[0] < 0.1) break;
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
		if(floor(Tr*1000) == floor(T*1000)) {
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

int *lab::sumDists(int *xdist, int *ydist) {
	int *result = new int[size];
	for(int i = 0; i < size; i++)
		result[i] = xdist[i] + ydist[i];
	return result;
}

