'''
Program to calculate the SJ state between perfectly reflecting mirrors
'''

import causalsets as cs
import scipy as sci
import numpy as np
np.set_printoptions(threshold=np.nan)
import sys
import matplotlib.pyplot as plt
import time
import scipy.linalg as linalg


def G_ret(coords,p):
	'''
	Calculates G-ret for p with no lightcone intersection with the mirrors.
	'''
	#p,q are labels of points
	#Only need points in the past of p
	
	xp,tp = p
	coords[0].append(xp)
	coords[1].append(tp)
	s_list = cs.sort(coords)
	#Initiate G matrix
	size = len(coords[0])
	G = np.zeros((size,size))
	position = np.where(s_list == (xp,tp))
	for i, xt in enumerate(s_list):
		xq = xt[0]
		tq = xt[1]
		if tq  == 0:
			break
		c = causal(xq,tq,xp,tp)
		if c:
			G[position][i] += -0.5
		else:
			None
			
	return G


def G_mirror(coords,p,l1,l2):
	'''

	Function to calculate G function with mirrors.

	'''
	xp,tp = p
	coords[0].append(xp)
	coords[1].append(tp)
	#Time order the co-ordinates but keep track of p.
	a = list(enumerate(coords[0]))
	b = list(enumerate(coords[1]))
	sorted_list = sorts(a,b)
	position = np.where(sorted_list == (xp,tp))
	t_intersection1 = tp - xp
	#Intersection of lightcone with left hand mirror
	mirror_intL = l1 + t_intersection1
	t_intersection2 = tp + xp
	#Intersection of lightcone with right hand mirror
	mirror_intR = t_intersection2 - l2

	#If neither mirror is in the past lightcone of p
	if mirror_intR <= 0 and mirror_intL <=0:
		Gret = G_ret(coords, p)
		return Gret

	else:
		size = len(coords[0])
		G = np.zeros((size,size))	
		for i, xt in enumerate(sorted_list):
			xq, tq = xt[0], xt[1]
			causalcheck = causal(xq,tq,xp,tp)
			if tp - tq == 0:
				break
			print causalcheck
			if causalcheck:
				G[position][i] += -0.5
				#Check if also in lightcone of image
				if mirror_intR <=0 and mirror_intL >=0:
					#Only left mirror is in the past lightcone of p
					imx = (2*l1) - xp
					imagecheck = causal(xq,tq,imx,tp)
				elif mirror_intR >=0 and mirror_intL <=0:
					#Only right mirror is in the past lightcone of p
					imx = (2*l2) - xp
					imagecheck = causal(xq,tq,imx,tp)
				if imagecheck:
					G[position][i] += 0.5
					print 'Zero'
				else:
					print 'Nonzero'
			else:
				None
					 
	return G

		
def causal(xq,tq,xp,tp):
	'''
	
	Function to test causal relation between two points
	
	'''
	
	if xp - xq == 0:
		return True
	elif abs((tp-tq)/(xp-xq)) >= 1:
		return True
	else:
		return False
	
def getkey(item):
	'''
	For sorting by time values.
	'''
	return(item[1])

def sorts(ex,et):
	'''
	Time orders the co-ordinates.
	Returns a list of co-ordinates in [x,t] form i.e. list of co-ordinate pairs.
	'''
	
	coords = []
        print ex
	for i in sprinkle[0]:
		x = sprinkle[0][i]
		t = sprinkle[1][i]
		coords.append([x,t])
	time_ordered = sorted(coords, key=getkey)
	return time_ordered


if __name__ == '__main__':
	tstart = time.time()
	#Ensure that X is much larger than L2-L1 and that L2-L1 is much larger than T
	sprinkling,X,T,L1,L2,xp,tp = int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]),int(sys.argv[5]),int(sys.argv[6]),int(sys.argv[7])

	#Make rectangle of Minkowski space
	Cs = cs.square(X,T, sprinkling)
	p = xp, tp
	Green = G_mirror(Cs,p,L1,L2)
	GT = Green.transpose()
	delta = 1j*(Green - GT)
	eigenvalues, eigenvectors = linalg.eig(delta)
	tend = time.time()
	print np.where(Green != 0)
	print 'Non-zero eigenvalues:'
	for i,e in enumerate(eigenvalues):
		#Account for 10^-15 machine error
		if abs(e.real) >= 0.000000000000001:
			print e.real
			#print eigenvectors[i]
	'''	
	ts = Cs[1]
	data = open('Casdata.dat','w')
	for i, x in enumerate(Cs[0]):
		eggy = abs(eigenvectors[0][i])
		dat = '%s %s %s \n' % (x,ts[i],eggy)
		data.write(dat)
	data.close()
	'''
	
	print 'Computation time:', tend-tstart, 'seconds'
	
	
