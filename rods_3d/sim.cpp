#include "sim.h"

rods::rods(int size, double betaVal, double density) {

	L = size;
	d = density;
	energy = 0;
	totalx = totaly = 0;

	if(d > 1 || d < 0) {
		d = 0.2;
		std::cerr << "density is a real value in [0,1]" << std::endl;
	}
	A = L*L;
	N = L*L*L;
	beta = betaVal;
	prob = new double[3];

	//initialize probability lookup table
	for(int i = 1; i < 3; i++) {
		prob[i] = exp(-beta*i);
		//std::cout << "prob " << i << ":" << prob[i] << std::endl;
	}

	//alloc lattice
	lat = new char[N];
	memset(lat, 0, N*sizeof(char));

	//lista de posicoes com rods
	numPos = (int) (N*d);
	pos = new int[numPos];

	//initialize lattice
	for(int i = 0; i < numPos; i++) {
		int n = 0;
		do {
			n = genrand_int32() % N;
		} while(lat[n] != 0);
		lat[n] = this->randDirection();
		pos[i] = n;
		//std::cout << i << ": Tossed " << (short)lat[n] << " at " << n << std::endl;
	}

	energy = recalcEnergy();

}

rods::~rods() {
	delete[] pos;
	delete[] lat;
	delete[] prob;
}

char rods::randDirection() {
	int magicNumber = genrand_int32() % 3;
	if(magicNumber == 0) return 1;
	else if(magicNumber == 1) return 2;
	else if(magicNumber == 2) return 4;
	return 0;
}

void rods::sweep() {
	int nn;
	int max = numPos;
	int count = 0;

	for(int k = 0; k < max; k++) {
		//eij energia na posicao j coord i
		int i1, i2, ex1, ex2, ey1, ey2, ez1, ez2, ei = 0, ef = 0;
		int v1, oldv1;
		int nPos = genrand_int32() % numPos;
		i1 = pos[nPos];
		//choose an empty site
		do {
			i2 = genrand_int32() % N;
		} while (lat[i2] != 0);

		oldv1 = v1 = lat[i1];

		//calc energies at position 1
		if( (nn=i1+1) >= N ) nn-=N;
		ex1 = lat[nn] & 1;
		if( (nn=i1-1) < 0 ) nn+=N;
		ex1 += lat[nn] & 1;
		if( (nn=i1+L)>=N ) nn-=N;
		ey1 = (lat[nn] & 2) >> 1;
		if( (nn=i1-L)<0 ) nn+=N;
		ey1 += (lat[nn] & 2) >> 1;
		if( (nn=i1+A)>=N ) nn-=N;
		ez1 = (lat[nn] & 4) >> 2;
		if( (nn=i1-A)<0 ) nn+=N;
		ez1 += (lat[nn] & 4) >> 2;

		if(v1 == 1) ei = ex1;
		else if(v1 == 2) ei = ey1;
		else ei = ez1;

		//flip orientation 0.5 chance
		if(genrand_real1() > 0.5) {
			if(v1 == 1) v1 = (genrand_real1() > 0.5)?2:4;
			else if(v1 == 2) v1 = (genrand_real1() > 0.5)?1:4;
			else v1 = (genrand_real1() > 0.5)?2:1;
		}

		//switch to calc energy
		lat[i2] = v1;
		lat[i1] = 0;

		//calc energies at position 2
		if( (nn=i2+1) >= N ) nn-=N;
		ex2 = lat[nn] & 1;
		if( (nn=i2-1) <  0) nn+=N;
		ex2 += lat[nn] & 1;
		if( (nn=i2+L)>=N ) nn-=N;
		ey2 = (lat[nn] & 2) >> 1;
		if( (nn=i2-L)<0 ) nn+=N;
		ey2 += (lat[nn] & 2) >> 1;
		if( (nn=i2+A)>=N ) nn-=N;
		ez2 = (lat[nn] & 4) >> 2;
		if( (nn=i2-A)<0 ) nn+=N;
		ez2 += (lat[nn] & 4) >> 2;

		if(v1 == 1) ef = ex2;
		else if(v1 == 2) ef = ey2;
		else ef = ez2;

		//calculate energy difference
		// -1 * (ef - ei)
		int delta = ei - ef;

		/******************************************************
		Metropolis algortihm
		the logic is reversed because we already did the switch
		and we want to keep our memory accesses to a minimum
		*******************************************************/
		if(delta > 0 && genrand_real1() > prob[delta]) {
			//undo switch
			lat[i2] = 0;
			lat[i1] = oldv1;
		} else {
			//keep the switch
			pos[nPos] = i2;
			energy += 2*delta;
			count++;
		}
	}
}

void rods::run() {
	this->sweep();
}

int rods::recalcEnergy() {
	int en = 0;
	//calculate initial energy
	for(int i = 0; i < numPos; i++) {
		int nn, n = pos[i];
		int sum = 0;
		//calculate the sum of the neighbouring spins
		if(lat[n] == 1) {
			/*
			1 & 1 == 1
			1 & 2/0 == 0
			*/
			if( (nn=n+1) >= N ) nn-=N;
			sum = lat[nn] & 1;
			if( (nn=n-1) < 0 ) nn+=N;
			sum += lat[nn] & 1;
		} else if (lat[n] == 2) {
			/*
			2 & 2 == 2
			2 & 1/0 == 0
			*/
			if( (nn=n+L)>= N) nn-=N;
			sum = (lat[nn] & 2) >> 1;
			if( (nn=n-L)<0 ) nn+=N;
			sum += (lat[nn] & 2) >> 1;
		} else {
			/*
			4 & 4 == 4
			4 & 2/1/0 == 0
			*/
			if( (nn=n+A)>=N ) nn-=N;
			sum = (lat[nn] & 4) >> 2;
			if( (nn=n-A)<0 ) nn+=N;
			sum += (lat[nn] & 4) >> 2;
		}

		//calc energy
		en -= sum;
	}

	return en;
}

void rods::showNeighbors(int pos) {
	int nn;
	std::cout << " ";
	if( (nn=pos-L)>=N ) nn-=N;
	std::cout << (int)lat[nn] << std::endl;
	if( (nn=pos-1) >= N ) nn-=N;
	std::cout << (int)lat[nn] << (int)lat[pos];
	if( (nn=pos+1) < 0 ) nn+=N;
	std::cout << (int)lat[nn] << std::endl;
	if( (nn=pos+L)<0 ) nn+=N;
	std::cout << " " << (int)lat[nn] << std::endl;
}

double *rods::getMRatio() {
	int *count = new int[3];
	double *result = new double[2];
	int buffer;
	count[0] = count[1] = count[2] = 0;

	for(int i = 0; i < numPos; i++) {
		if(lat[pos[i]] == 1) count[0]++;
		else if(lat[pos[i]] == 2) count[1]++;
		else count[2]++;
	}

	for(int i = 0; i < 3; i++) {
		for(int j = i+1; j < 3; j++) {
			if(count[j] > count[i]) {
				buffer = count[i];
				count[i] = count[j];
				count[j] = buffer;
			}
		}
	}

	result[0] = fabs(count[0] - count[1])/(double)numPos;
	result[1] = fabs(count[0] - count[2])/(double)numPos;

	delete count;

	return result;
}

double rods::calcXLen() {
	int l = 0, suml = 0, total = 0;

	for(int i = 0; i < N; i++) {
		if(lat[i] == 1 && l != 0) {
			l++;
		}
		if(lat[i] == 1 && l == 0) {
			total++;
			l++;
		}
		if(lat[i] != 1 && l != 0) {
			suml += l;
			l = 0;
		}
	}
	totalx = total;
	return suml;
}

double rods::calcYLen() {
	int l = 0, suml = 0, total = 0;

	for(int i = 0; i < L; i++) {
		int zPos = A * i;
		for(int j = 0; j < L; j++) {
			int xPos = zPos + j;
			for(int n = 0; n < L; n++) {
				int curPos = xPos + n * L;
				int firstPos = xPos;
				if(lat[curPos] == 2) {
					if(l == 0) total++;
					l++;
				}
				if(lat[curPos] != 2 && l != 0) {
					suml += l;
					l = 0;
				}
				//end
				if(n == (L - 1)) {
					if(lat[curPos] == 2 && lat[firstPos] == 2) {
						total--;
					}
					suml += l;
					l = 0;
				} 
			}
		}
	}
	
	totaly = total;
	
	return suml;
}

double rods::calcZLen() {
	int l = 0, suml = 0, total = 0;

	for(int i = 0; i < L; i++) {
		int yPos = L * i;
		for(int j = 0; j < L; j++) {
			int xPos = yPos + j;
			for(int n = 0; n < L; n++) {
				int curPos = xPos + n * A;
				int firstPos = xPos;
				if(lat[curPos] == 4) {
					if(l == 0) total++;
					l++;
				}
				if(lat[curPos] != 4 && l != 0) {
					suml += l;
					l = 0;
				}
				//end
				if(n == (L - 1)) {
					if(lat[curPos] == 4 && lat[firstPos] == 4) {
						total--;
					}
					suml += l;
					l = 0;
				}
			}
		}
	}

	totalz = total;

	return suml;
}

double rods::getAvxLen() {
	if(totalx == 0) return 0;
	return xLen/(double)totalx;
}

double rods::getAvyLen() {
	if(totaly == 0) return 0;
	return yLen/(double)totaly;
}

double rods::getAvzLen() {
	if(totalz == 0) return 0;
	return zLen/(double)totalz;
}

double rods::getAvLen() {
	xLen = this->calcXLen();
	yLen = this->calcYLen();
	zLen = this->calcZLen();

	return(xLen+yLen+zLen)/(double)(totalx+totaly+totalz);
}

int *rods::getXDist() {
	//int l = 0;
	int *dist = new int[L+1];
	memset(dist, 0, (L+1)*sizeof(int));

	/*for(int i = 0; i < N; i++) {
		if((i+1) % L == 0) {
			//end of the line
			if(lat[i] == 1) {
				l++;
			}
			if(lat[i] == 1 && lat[(i+1) - L] == 1) {
				int ii = (i+1) - L;
				int templ = 0;
				//recalc first string
				while(lat[ii] == 1){
					templ++;
					ii++;
				}
				if(templ < 100) {
					//subtract first string from dist
					dist[templ]--;
					//and add it to this one
					l += templ;
				}
			}
			dist[l]++;
			l = 0;
		} else {
			if(lat[i] == 1) {
				l++;
			}
			if(lat[i] != 1 && l != 0) {
				dist[l]++;
				l = 0;
			}
		}
	}*/

	return dist;
}

int *rods::getYDist() {
	//int l = 0;
	int *dist = new int[L+1];
	memset(dist, 0, (L+1)*sizeof(int));

	/*for(int i = 0; i < L; i++) {
		for(int j = 0; j < L; j++) {
			int curPos = i + L * j;
			
			//j e o indice y, a primeira condicao testa a fronteira do y
			//i sera a primeira posicao da nossa coluna
			//(ou seja i + L*j com j=0)
			*/
			/*if((j+1) % L == 0) {
				if(lat[cpos] == 2) l++;
				if(lat[cpos] == 2 && lat[i] == 2) {
					int ii = i;
					int templ = 0;
					//recalc first string
					while(lat[ii] == 2) {
						templ++;
						ii += L;
					}
					if(templ < 100) {
						//subtract first string from dist
						dist[templ]--;
						//and add it to this one
						l += templ;
					}
				}
				dist[l]++;
				l = 0;
			} else {
				if(lat[cpos] == 2) {
					l++;
				}
				if(lat[cpos] != 2 && l != 0) {
					dist[l]++;
					l = 0;
				}
			}
		}
	}*/

	return dist;
}

int rods::getPredDirection() {
	if(totalx == totaly && totalx == totalz) return 0;
	if(totalx > totaly) {
		if(totalx > totalz) return 1;
		else return 3;
	} else {
		if(totaly > totalz) return 2;
		else return 3;
	}
	return 0;
}