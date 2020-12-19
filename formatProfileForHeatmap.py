#!/usr/bin/python
# Will organize profiles taken from Fiji for making heatmaps in four channels
# Will return one file for every channel

import sys, os, re
import numpy as np
from scipy.signal import argrelextrema

def main(args):
	if not len(args) == 3: sys.exit("USAGE: python formatProfileForHeatmap.py profiles filePrefix")

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

	# MaxMin normalize
	indexI = 0
	for i in cellList:
		indexJ = 0
		for j in i:
			j = [float(x) for x in j]
			jmin, jmax = min(j), max(j)
			for k, val in enumerate(j):
				j[k] = (val - jmin) / (jmax - jmin)
			i[indexJ] = j
			indexJ += 1
		cellList[indexI] = i
		indexI += 1
		
	#### Now each list will be normalized between 0 and 1 ####

	# Flip so that brightest speckle is on the right
	indexCellList = 0
	for cell in cellList:
		# get index of brightest speckle pixel
		indexMax = cell[0].index(1.0)
		if indexMax < len(cell[0])/2:
			indexChannel = 0
			for channel in cell:
				channel.reverse()
				cell[indexChannel] = channel
				indexChannel += 1
		cellList[indexCellList] = cell
		indexCellList += 1

	# Sort by closest p53 local max to center and calculate the percentage of cells that have p53 local maximum within 2 pixels of SC35 max
	distToCenterList = []
	cellsWithin2 = 0
	numAtSpeckle = 0
	withinNotSpeckle = 0
	nearCenter = 0
	for cell in cellList:
		x = np.array(cell[2])
		y = argrelextrema(x, np.greater)
		#take only top four local maximum
		if len(y[0]) > 4:
			top3 = []
			for i in y[0]:
				# is cell[2][i] in top 3
				if len(top3) < 4:
					top3.append(i)
				else:
					value1 = cell[2][top3[0]]
					value2 = cell[2][top3[1]]
					value3 = cell[2][top3[2]]
					value4 = cell[2][top3[3]]
					values = [value1, value2, value3, value4]
					indexMinValue = values.index(min(values))
					if cell[2][i] > value1 or cell[2][i] > value2 or cell[2][i] > value3 or cell[2][i] > value4:
						top3[indexMinValue] = i
			y = top3
		else:
			y = y[0].tolist()
			
		center = len(cell[2])/2
		distToCenter = 30
		for i in y:
			dist = abs(i - center)
			if dist < distToCenter:
				distToCenter = dist
		distToCenterList.append(distToCenter)
		if distToCenter <= 2:
			nearCenter += 1
		
		indexMax = cell[0].index(max(cell[0]))
		isWithin = False
		for i in y:
			if abs(i-indexMax) <= 2:
				isWithin = True
		if isWithin:
			cellsWithin2 += 1
		
		
	print distToCenterList	
		
	sortedCellList = [x for (y,x) in sorted(zip(distToCenterList, cellList), key=lambda pair: pair[0])]
	
	print cellsWithin2
	print nearCenter
	print len(cellList)
#	print withinNotSpeckle
#	print len(cellList) - numAtSpeckle
#	print numAtSpeckle

	
	# Write output
	# outfile names
	outSC35 = args[2] + "_SC35.txt"
	outExon = args[2] + "_exon.txt"
	outP53 = args[2] + "_p53.txt"
	outIntron = args[2] + "_introns.txt"
	sc35 = open(outSC35, 'w')
	exon = open(outExon, 'w')
	p53 = open(outP53, 'w')
	intron = open(outIntron, 'w')
	for cell in sortedCellList:
		print >>sc35, "\t".join(str(x) for x in cell[0])
		print >>exon, "\t".join(str(x) for x in cell[1])
		print >>p53, "\t".join(str(x) for x in cell[2])
		print >>intron, "\t".join(str(x) for x in cell[3])
		
		
	
	
	

		



if __name__ == "__main__": main(sys.argv)
