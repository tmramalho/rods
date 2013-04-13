#! /usr/bin/python
# -*- coding: utf-8 -*-

import os
import re
import math

def calcAverage(tList):
	acc = 0
	total = 0
	for i in tList:
		acc += float(i)
		total += 1
	return acc/total

def calcStdDev(tList):
	acc = 0
	total = 0
	mean = calcAverage(tList)
	for i in tList:
		acc += (float(i) - mean)*(float(i) - mean)
		total += 1
	return math.sqrt(acc/(total-1))

def calcError(tList):
	acc = 0
	total = 0
	mean = calcAverage(tList)
	for i in tList:
		acc += (float(i) - mean)*(float(i) - mean)
		total += 1
	return math.sqrt(acc/(total-1)/total)

if __name__ == '__main__':
	dataStr = []
	data = {}
	avData = {}
	total = 0

	#read data from individual directories
	for root, dirs, files in os.walk('./r3dout'):
		for name in files:
			if name == "master.dat":
				filename = os.path.join(root, name)
				with open(filename) as f:
					f.readline()
					dataStr.append(f.readline())
					total+=1

	#split strings and create arrays with individual
	#values for each temperature
	for i in dataStr:
		tSet = i.split()
		if not tSet[0] in data:
			data[tSet[0]] = [[],[],[],[],[],[],[]]
			avData[tSet[0]] = []
		#c
		data[tSet[0]][0].append(float(tSet[4]))
		#m1
		data[tSet[0]][1].append(float(tSet[5]))
		#\xsim1
		data[tSet[0]][2].append(float(tSet[6])-float(tSet[5])*float(tSet[5]))
		#m2
		data[tSet[0]][3].append(float(tSet[7]))
		#\xsim2
		data[tSet[0]][4].append(float(tSet[8])-float(tSet[7])*float(tSet[7]))
		#l
		data[tSet[0]][5].append(float(tSet[9]))
		#\xsil
		data[tSet[0]][6].append(float(tSet[10])-float(tSet[9])*float(tSet[9]))

	#calc averages of each parameter for each temperature
	#also errors
	for key in data:
		avData[key].append(calcAverage(data[key][0]))
		avData[key].append(calcError(data[key][0]))
		avData[key].append(calcAverage(data[key][1]))
		avData[key].append(calcError(data[key][1]))
		avData[key].append(calcAverage(data[key][2]))
		avData[key].append(calcError(data[key][2]))
		avData[key].append(calcAverage(data[key][3]))
		avData[key].append(calcError(data[key][3]))
		avData[key].append(calcAverage(data[key][4]))
		avData[key].append(calcError(data[key][4]))
		avData[key].append(calcAverage(data[key][5]))
		avData[key].append(calcError(data[key][5]))
		avData[key].append(calcAverage(data[key][6]))
		avData[key].append(calcError(data[key][6]))

	with open("mineresult.dat", "w" ) as f:
		for key in sorted(avData.iterkeys()):
			f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (tuple([key]) + tuple(avData[key])))

	os.system("gnuplot < plotcena.p")
	
	print "Read " + str(total) + " files"