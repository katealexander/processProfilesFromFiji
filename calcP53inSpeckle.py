#!/usr/bin/python
# Binarization method
# Will organize profiles taken from Fiji for making heatmaps in four channels
# Will return one file for every channel

import sys, os, re
import numpy as np
from scipy.signal import argrelextrema

def main(args):
	if not len(args) == 2: sys.exit("USAGE: python calcP53inSpeckle.py profiles")
	
	# Read in data; list of lists; [[SC35], [exons], [p53], [introns]]
	f = open(args[1])
	# read first line
	line = f.readline()[:-1]
	sc35 = line.split("\t")[1]
	exons = line.split("\t")[3]
	p53 = line.split("\t")[5]
	introns = line.split("\t")[7]
	cell = [[sc35], [exons], [p53], [introns]]
	line = f.readline()[:-1]
	cellList = []
	while line != "":
		sc35 = line.split("\t")[1]
		exons = line.split("\t")[3]
		p53 = line.split("\t")[5]
		introns = line.split("\t")[7]
		if line.split("\t")[0] == "0":
			#append previous cell to cellList
			cellList.append(cell)
			#start a new cell
			cell = [[sc35], [exons], [p53], [introns]]
		else:
			#append to current cell
			cell[0].append(sc35)
			cell[1].append(exons)
			cell[2].append(p53)
			cell[3].append(introns)
		line = f.readline()[:-1]
	f.close()

	#append final cell to cellList
	cellList.append(cell)
	
	#### Now each list in cellList is one cell ####

	# MaxMin normalize and then binarize speckle channel
	indexI = 0
	for i in cellList:
		indexJ = 0
		for j in i:
			j = [float(x) for x in j]
			jmin, jmax = min(j), max(j)
			for k, val in enumerate(j):
				j[k] = (val - jmin) / (jmax - jmin)
				if indexJ == 0: 
					if j[k] > 0.5: #theshold to be called a speckle
						j[k] = 1.0
					else:
						j[k] = 0.0
				i[indexJ] = j
			indexJ += 1
		cellList[indexI] = i
		indexI += 1
		
	#### Now speckle list will be binarized ####

	# calculate amount of p53 in speckle
	for cell in cellList:
		speckle = cell[0]
		notInSpeckle = []
		for i in speckle:
			if i == 1.0:
				notInSpeckle.append(0.0)
			else:
				notInSpeckle.append(1.0)
		p53 = cell[2]
		amountSpeckle = sum(speckle)
		amountNotSpeckle = sum(notInSpeckle)
		amountP53inSpeckle = sum([a*b for a, b in zip(speckle, p53)])
		amountP53notInSpeckle = sum([a*b for a, b in zip(notInSpeckle, p53)])
		relP53inSpeckle = amountP53inSpeckle/float(amountSpeckle)
		relP53notInSpeckle = amountP53notInSpeckle/float(amountNotSpeckle)
		print str(relP53inSpeckle) + "\t" + str(relP53notInSpeckle) + "\t" + str(relP53inSpeckle/relP53notInSpeckle)



if __name__ == "__main__": main(sys.argv)
