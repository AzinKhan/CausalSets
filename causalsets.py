'''
Sprinkling

'''
import numpy as np
import numpy.ma as ma
#np.set_printoptions(threshold='nan')
import scipy
import scipy.linalg as linalg
#from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import gaussian_kde
import sys
import time
thresh = 0.000000000000001
identify = True

def square(xmax,tmax,size):
	'''
	size = Number of points sprinkled into initial rectangle
	'''
	#Make rectangle:
	tlist = []
	xlist = []
	for p in range(0,size): #Generate random x,t values
		x = np.random.randint(-xmax, xmax) #x,t are integer valued
		t = np.random.randint(0,tmax)
		xlist.append(x)
		tlist.append(t)	
	return xlist, tlist
	
def split(square,G):
    #G is boost parameter (gradient)
	xcoords = []
	tcoords = []
	tmax = max(square[1])
	for i,s in enumerate(square[0]):
		x = s
		t = square[1][i]
		if x >= 0 and t > (G*x) and t < (tmax - x):#Top boundary of kite is null line
			xcoords.append(x)
			tcoords.append(t)
		elif x <= 0 and t > (G*-1*x) and t < (tmax + x):
			xcoords.append(x)
			tcoords.append(t)
		else:
			None
	return xcoords, tcoords

def makekite(x,t,s,g):

	'''
	Function to create randomly distributed points within a 2-D flat space.
	x = largest x co-ordinate value
	t = largest t co-ordinate value
	s = number of points sprinkled on
	g = gradient of lower edge of kite (boost parameter)
	
	The upper edges of the shape are null lines.

	If g = 1, and x = t obtain a Minkowski space diamond.

	For 'Kite': g > 1 and t > x
	
	Returns X and T co-ordinates as two lists.

	'''

	init = square(x,t,s) # Make initial square.
	xs,ts = split(init,g)
	return xs,ts

def getkey(item):
	'''
	For sorting by time values.
	'''
	return(item[1])

def sort(sprinkle):
	'''
	Time orders the co-ordinates.
	Returns a list of co-ordinates in [x,t] form i.e. list of co-ordinate pairs.
	'''
	coords = []
	for i, s in enumerate(sprinkle[0]):
		x = sprinkle[0][i]
		t = sprinkle[1][i]
		coords.append([x,t])
	time_ordered = sorted(coords, key=getkey)
	return time_ordered
	
def causalmatrix(kite, G):
	ordered = sort(kite)
	#orderedcopy = sort(kite)
	size = len(ordered)
	C = np.zeros((size,size)) #Initialize causal matrix
	
	for i,f in enumerate(ordered):
		x1 = f[0]
		t1 = f[1]
		for k in range(i,len(ordered)): #Iterate only over future lying points
			x2 = ordered[k][0]
			t2 = ordered[k][1]
			# Find intersection of lightcone and edge of kite
			if x1 < 0:
				tint = (G/(G+1))*(t1 + x1)
			else:
				tint = (G/(G+1))*(t1 - x1)
			xint = tint*(1/G)
			if i == k:
				C[i][k] = 1
			elif t1 == t2:
				None
			elif x2 != x1 and abs((t2-t1)/(x2-x1)) >= 1:
				C[i][k] = 1
			elif x2 == x1:
				C[i][k] = 1
			#Extra part for identification of edges.
			if identify:
				if x2 != (-1*xint) and (t2-t1)/(x2-(-1*xint)) >= 1:
					C[i][k] = 1
					
	return C

def eigenvalue(C):

	'''
	Takes causal matrix and determines iSJ function and its eigenvalues and eigenvectors
	'''

	# Massless scalar field
	R = 0.5*C
	Rt = np.transpose(R)
	SJ = 1j*(R - Rt)
	evalues, rv = linalg.eig(SJ) # Determine eigenvectors and eignenvalues of iSJ function

        evectors = []
        # DO THIS FUCKING THING
        for i, value in enumerate(evalues):
                evectors.append(rv[:,i])
	
	return evalues, evectors

def SJWhight(VV):

	'''
	Computes SJ Whightman function from given eigenvalues and eigenvectors.
	Use output from function(eigenvalue) as input. The eigenvectors are normalized within the function.
	'''
	evectors = VV[1]
	values = VV[0].real
	#Iterate over vectors to create spectral projection
	length = len(evectors)
	#Initialize Q
	Q = np.zeros((length,length), dtype=complex)
	vectors = []
	#Normalize the eigenvectors
	for i, v in enumerate(evectors):
		norm = linalg.norm(v)
		normvect = v/norm
		vectors.append(normvect)

	for i,v in enumerate(vectors):
		eigenvalue = values[i]
		#Only require eigenvectors with positive eigenvalues
		if eigenvalue >= 0:
			#Convert to matrix object
			vec = np.matrix(v)
			vecH = vec.H
			#Calculate outer product of the two vectors
			outer = np.outer(vec,vecH)
			proj = outer*eigenvalue
			Q += proj
		else:
			None

	return Q 

def image(X,Y,Z):

	'''
	Creates an image array from with XY co-ordinates and associated Z value
	'''

	xsize = 2*max(X)
	ysize = max(Y)
	a = np.zeros((xsize,ysize))
	for i, x in enumerate(X):
		xval = int(x + (0.5*xsize))
		yval = int(Y[i])
		value = Z[i]
		a[(xval-1)][(yval-1)] = value

	return a

def Entropy(C_matrix):
	'''
	Take causal matrix and return entropy of space-time as per Sorkin's paper.
	'''
	EV = eigenvalue(C_matrix)
	W = SJWhight(EV)
	R = 0.5*C_matrix
	Rt = np.transpose(R)
	sj = 1j*(R - Rt)

	sjmask, Wmask = zeroprojection(sj,W)

	#Initialize entropy
	S = 0
	#Now calculate eigenvectors of new

	evals,evectors = linalg.eig(Wmask,sjmask, overwrite_a=False, overwrite_b=False)


	for p, l in enumerate(evals):
		if abs(l) < thresh:
			print 'Problem.'
			
		elif np.isnan(l):
			#Remove any NaN values
			print 'nan'
			None
		else:
			L = abs(l)
			ent = l*np.log(L)
			S += ent 
		
	return S
		
	
def zeroprojection(sj,W):
	'''

	Function to constrain image of SJ i.e. projects out any vectors with zero eigenvalue.
	Changes representation basis to that of eigenvectors of SJ.

	'''
	evalues, evectors = linalg.eig((sj))
	SJeigs = sorted(evalues,reverse=True)
	Weigs = sorted(linalg.eigvals(W),reverse=True)
	shape = sj.shape
	sj_init = np.zeros((shape), dtype='cfloat')
	P_init = np.zeros((shape), dtype='cfloat')
	#Change basis to diagonalize SJ matrix.
	Scounter = 0 
	for i, e in enumerate(SJeigs):
		if e > thresh:
			sj_init[i][i] = e.real
		else:
			Scounter = i
			break

	#Transform W to this basis by constructing transformation matrices
	evectlist = []
	#Make list of tuples with eigenvalues and associated eigenvectors
	for i, e in enumerate(evalues):
                tup = (evectors[:,i],e)
		evectlist.append(tup)
	#Sort the eigenvalue/vector list to be in the same order as those in SJ 
	finallist = sorted(evectlist,reverse=True, key=getkey)
	
	for i, x in enumerate(finallist):
		P_init[:,i] = x[0]
	
	P = P_init.transpose()
	#Find inverse of P
	shapep = P.shape
	Pinv = linalg.solve(P, np.identity(shapep[0]))
	#Pinv = linalg.inv(P)
	Wtransformed = P*W*Pinv
        
        '''
	#Transformation check
	check = sj_init - (P*sj*Pinv)
	for i in range(shape[0]):
		#These should all be zero
		print check[i][i]
	'''

	Wtransformed = np.zeros((shape), dtype='cfloat')
	for h, w in enumerate(Weigs):
		Wtransformed[h][h] = w.real
	sj_diag = sj_init[:Scounter,:Scounter]
	W_diag = Wtransformed[:Scounter,:Scounter]
	
	return sj_diag, W_diag
		
	
if __name__ == '__main__':
	# User input from command line
	tstart = time.time()
	xsize, tsize, sprinklesize, gradient, vectornum = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
	xs,ts = makekite(xsize,tsize,sprinklesize,gradient)
	kite = xs, ts
	CM = causalmatrix(kite, gradient)
	vals,vectors = eigenvalue(CM)
	sj = SJWhight((vals,vectors))
	eigenvals = vals[0].real
	s = Entropy(CM)
	print 'Entropy:', s
	eigenvectors = []
	#Normalize the eigenvectors
	for i, v in enumerate(vectors):
		norm = linalg.norm(v)
		normvect = v/norm
		eigenvectors.append(normvect)
	'''	
	data = open('yarmulkemode.dat','w')
	for i, x in enumerate(xs):
		eigenvectorval = 1*abs(eigenvectors[vectornum][i])
		dat = '%s %s %s \n' % (x,ts[i],eigenvectorval)
		data.write(dat)
	data.close()

	data1 = open('Causalmode10.dat','w')
	for i, x in enumerate(xs):
		eigenvectorval = 100000*abs(vectors[10][i])
		dat1 = '%s %s %s \n' % (x,ts[i],eigenvectorval)
		data1.write(dat1)
	data1.close()

	data2 = open('Causalmode20.dat','w')
	for i, x in enumerate(xs):
		eigenvectorval = 100000*abs(vectors[20][i])
		dat2 = '%s %s %s \n' % (x,ts[i],eigenvectorval)
		data2.write(dat2)
	data2.close()
	'''
	
	data = open('test.dat', 'w')
	for i, x in enumerate(xs):
		t = ts[i]
		s = abs(sj[100][i])
		dat = '%s %s %s \n' % (x,t,s)
		data.write(dat)
	data.close()
	
	#print 'First 3 eigenvalues:', vals[0],vals[2], vals[4]
	tend = time.time() - tstart
	print 'Calculation time:', tend, 'seconds'
