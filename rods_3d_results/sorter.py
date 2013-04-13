#! /usr/bin/python
# -*- coding: utf-8 -*-

import os
import re

def calcAverage(tList):
	acc = 0
	total = 0
	for i in tList:
		acc += float(i)
		total += 1
	return acc/total

temp = re.compile("Running simulation at (.+)?",re.I | re.U)
cTime = re.compile("Correlation time: (.+)?",re.I | re.U)
cenas = []
#os.getcwd()
for root, dirs, files in os.walk("./logs"):
	for name in files:
		filename = os.path.join(root, name)
		#print filename
		with open(filename) as f:
			data = f.read()
			try:
				temperature = temp.search(data).group(1)
			except AttributeError:
				print "fail"
				continue
			ctimes = cTime.findall(data)
			if len(ctimes) == 0:
				continue
			try:
				val = calcAverage(ctimes)
			except ValueError:
				print "fail"
			else:
				cenas.append(str(temperature) + " " + str(val))
cenas.sort()
with open("ctimes.dat", 'a') as f:
	for i in cenas:
		f.write(i + "\n")
