#! /usr/bin/python
# -*- coding: utf-8 -*-
import os

rho = 0.2
boxSize = 20
numSamples = 30

#for k in range(0,30):
i = 0.1
while i <= 0.3:
	if(i<0.18):
		eqTime = 500000
	else:
		eqTime = 200000
	command = "qsub ~/msimple ~/r3dout ~/r3data ~/r3d/rods 1 "+ str(i) + " " + str(i) + " 0.001 " + str(rho) + " " + str(eqTime) + " " + str(boxSize) + " " + str(eqTime) + " " + str(numSamples)
	print command
	os.system(command)
	if(i<0.2):
		i+=0.001
	else:
		i+=0.01
