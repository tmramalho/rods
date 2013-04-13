#include "sim.h"

rods::rods(int size, double betaVal, double density) {

	L = size;
	d = density;
	energy = 0;
	mag = 0;

	if(d > 1 || d < 0) {
		d = 0.2;
		std::cerr << "density is a real value in [0,1]" << std::endl;
	}
	N = L*L;
	beta = betaVal;
	prob = new double[6];

	//initialize probability lookup table
	for(int i = 2; i < 5; i += 2) {
		prob[i] = (double)exp(-2*beta*i);
		//std::cout << "prob " << i << ":" << prob[i] << std::endl;
	}

	//alloc lattice
	lat = new short[N];
	memset(lat, 0, N*sizeof(short));

	//initialize lattice
	for(int i=0; i < N; i++) {
		mag += lat[i] = genrand_real1() > 0.5?1:-1;
	}

	energy = recalcEnergy();

}

rods::~rods() {
	delete[] lat;
	delete[] prob;
}

void rods::sweep() {
	int i;
	int nn, sum, delta;
	short *s = lat; //quick hack

	for(int k=0; k<N; k++) {

		//choose a site
		i = N*genrand_real1();

		//calculate the sum of the neighbouring spins
		if ((nn=i+XNN)>=N) nn-=N;
		sum = s[nn];
		if ((nn=i-XNN)<0) nn+=N;
		sum += s[nn];
		if ((nn=i+YNN)>=N) nn-=N;
		sum += s[nn];
		if ((nn=i-YNN)<0) nn+=N;
		sum += s[nn];

		//calculate energy difference
		delta = sum*s[i];

		//decide whether to flip spin
		if(delta <= 0) {
			s[i] = -s[i];
			mag += 2*s[i];
			energy += 4*delta;
		} else if (genrand_real1()<prob[delta]) {
			s[i] = -s[i];
			mag += 2*s[i];
			energy += 4*delta;
		}
	}
}

void rods::run() {
	this->sweep();
}

int rods::recalcEnergy() {
	int en = 0;
	int nn;
	short *s = lat; //quick hack
	//calculate initial energy
	for(int i=0; i < N; i++) {
		int sum = 0;
		//calculate the sum of the neighbouring spins
		if ((nn=i+XNN)>=N) nn-=N;
		sum = s[nn];
		if ((nn=i-XNN)<0) nn+=N;
		sum += s[nn];
		if ((nn=i+YNN)>=N) nn-=N;
		sum += s[nn];
		if ((nn=i-YNN)<0) nn+=N;
		sum += s[nn];

		//calculate energy difference
		en -= sum*s[i];
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
