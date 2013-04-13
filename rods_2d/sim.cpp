#include "sim.h"

rods::rods(int size, double betaVal, double density, char *lattice) {

	L = size;
	d = density;
	energy = 0;
	totalx = totaly = 0;

	if(d > 1 || d < 0) {
		d = 0.2;
		std::cerr << "density is a real value in [0,1]" << std::endl;
	}
	N = L*L;
	beta = betaVal;
	prob = new double[3];

	//initialize probability lookup table
	for(int i = 1; i < 3; i++) {
		prob[i] = exp(-beta*i);
		//std::cout << "prob " << i << ":" << prob[i] << std::endl;
	}

	//lista de posicoes com rods
	numPos = (int) (N*d);
	pos = new int[numPos];

	if(lattice == NULL) {
		//alloc lattice
		lat = new char[N];
		memset(lat, 0, N*sizeof(char));

		//initialize lattice
		for(int i = 0; i < numPos; i++) {
			int n = 0;
			do {
				n = genrand_int32() % N;
			} while(lat[n] != 0);
			lat[n] = genrand_real1()>0.5?1:2;
			pos[i] = n;
			//std::cout << i << ": Tossed " << (short)lat[n] << " at " << n << std::endl;
		}

	} else {
		lat = lattice;
		int n = 0;

		for(int i = 0; i < N; i++) {
			if(lat[i] != 0) {
				pos[n] = i;
				n++;
			}
			if(n > numPos) {
				std::cerr << "Wrong lattice!" << std::endl;
				break;
			}
		}
	}

	energy = recalcEnergy();

}

rods::~rods() {
	delete[] pos;
	delete[] lat;
	delete[] prob;
}

void rods::sweep() {
	int nn;
	int max = numPos;

	for(int k = 0; k < max; k++) {
		//eij energia na posicao j coord i
		int i1, i2, ex1, ex2, ey1, ey2, ei = 0;
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

		ei = (v1 == 1)?ex1:ey1;

		//flip orientation 0.5 chance
		if(genrand_real1() > 0.5)
			v1 = (v1 == 1)?2:1;

		//temporary switch to calc energy
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

		//undo switch to calc energy
		lat[i2] = 0;
		lat[i1] = oldv1;

		int delta = 0;

		//calculate energy difference
		if(v1 == 1) {
			delta -= ex2 - ei;
		} else {
			delta -= ey2 - ei;
		}

		//decide whether to flip
		if(delta <= 0) {
			lat[i2] = v1;
			lat[i1] = 0;
			pos[nPos] = i2;
			energy += 2*delta;
		} else if(genrand_real1()<prob[delta]) {
			lat[i2] = v1;
			lat[i1] = 0;
			pos[nPos] = i2;
			energy += 2*delta;
		}
	}
}

void rods::run() {
	this->sweep();
}

int rods::getEnergy() {
	return energy;
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
		} else {
			/*
			2 & 2 == 2
			2 & 1/0 == 0
			*/
			if( (nn=n+L)>= N) nn-=N;
			sum = (lat[nn] & 2) >> 1;
			if( (nn=n-L)<0 ) nn+=N;
			sum += (lat[nn] & 2) >> 1;
		}

		//calculate energy difference
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

double rods::getMRatio() {
	int xcount = 0, ycount = 0;

	for(int i = 0; i < numPos; i++) {
		if(lat[pos[i]] == 1) xcount++;
		else ycount++;
	}

	double result = fabs(ycount - xcount)/(double)numPos;

	return result;
}

 double rods::calcXLen() {
	int l = 0, suml = 0, total = 0;

	for(int i = 0; i < N; i++) {
		if((i+1) % L == 0) {
			if(lat[(i+1) - L] == 1) {
				if(lat[i] == 1) l++;
				suml += l;
			} else {
				if(lat[i] == 1 && l != 0) {
					l++;
				}
				if(lat[i] == 1 && l == 0) {
					total++;
					l++;
				}
				suml += l;
			}
			l = 0;
		} else {
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
	}
	totalx = total;
	return suml;
}

double rods::calcYLen() {
	int l = 0, suml = 0, total = 0;

	for(int i = 0; i < L; i++) {
		for(int j = 0; j < L; j++) {
			int cpos = i + L * j;
			/*
			j e o indice y, a primeira condicao testa a fronteira do y
			i sera a primeira posicao da nossa coluna
			(ou seja i + L*j com j=0)
			*/
			if((j+1) % L == 0) {
				if(lat[i] == 2) {
					if(lat[cpos] == 2) l++;
					suml += l;
				} else {
					if(lat[cpos] == 2 && l != 0) {
						l++;
					}
					if(lat[cpos] == 2 && l == 0) {
						total++;
						l++;
					}
					suml += l;
				}
				l = 0;
			} else {
				if(lat[cpos] == 2 && l != 0) {
					l++;
				}
				if(lat[cpos] == 2 && l == 0) {
					total++;
					l++;
				}
				if(lat[cpos] != 2 && l != 0) {
					suml += l;
					l = 0;
				}
			}
		}
	}
	totaly = total;
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

double rods::getAvLen() {
	xLen = this->calcXLen();
	yLen = this->calcYLen();

	return(xLen+yLen)/(double)(totalx+totaly);
}

int *rods::getXDist() {
	int l = 0;
	int *dist = new int[L+1];
	memset(dist, 0, (L+1)*sizeof(int));

	for(int i = 0; i < N; i++) {
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
	}

	return dist;
}

int *rods::getYDist() {
	int l = 0;
	int *dist = new int[L+1];
	memset(dist, 0, (L+1)*sizeof(int));

	for(int i = 0; i < L; i++) {
		for(int j = 0; j < L; j++) {
			int cpos = i + L * j;
			/*
			j e o indice y, a primeira condicao testa a fronteira do y
			i sera a primeira posicao da nossa coluna
			(ou seja i + L*j com j=0)
			*/
			if((j+1) % L == 0) {
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
	}

	return dist;
}

int rods::getPredDirection() {
	if(totalx == totaly) return 0;
	return ((totalx > totaly) ? 1 : 2);
}