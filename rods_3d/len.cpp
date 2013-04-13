#include "len.h"

len::len(int maxLengthSize) {
	numSamples = 0;
	size = maxLengthSize;
	master = new int[size];
	memset(master, 0, size*sizeof(int));
}

len::~len() {
	delete[] master;
}

void len::addNewMeasurement(int *dist) {
	numSamples++;
	for(int i = 0; i < size; i++) {
		master[i] += dist[i];
	}
	delete dist;
}

double len::getAvLen(int pos) {
	if(pos < size) return (double)master[pos]/numSamples;
	else return 0;
}

