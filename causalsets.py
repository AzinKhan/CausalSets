'''
Sprinkling

'''
import numpy as np
import scipy.linalg as linalg
from scipy.stats import gaussian_kde
import sys
import time
thresh = 0.000000000000001
identify = True

def getkey(item):
    '''
    For sorting by time values.
    '''
    return(item[1])

class CausalSet(object):

    def __init__(self, xmax, tmax, sprinklenum, gradient):
        self.xmax = xmax
        self.tmax = tmax
        self.sprinklenum = sprinklenum
        self.gradient = gradient
        self.xcoords = []
        self.tcoords = []
        self.ordered = []
        self.entropy = 0
        self.causalmatrix = 0
        self.SJ = 0
        self.Q = 0


    def makesquare(self):
        '''
        size = Number of points sprinkled into initial rectangle
        '''
        #Make rectangle:
        tlist = []
        xlist = []
        for p in range(0,self.sprinklenum): #Generate random x,t values
            x = np.random.randint(-self.xmax, self.xmax) #x,t are integer valued
            t = np.random.randint(0,self.tmax)
            xlist.append(x)
            tlist.append(t) 
        self.square = xlist, tlist
        
    def split(self):
        #G is boost parameter (gradient)
        G = self.gradient
        xcoords = []
        tcoords = []
        square = self.square
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
        self.xcoords = xcoords
        self.tcoords = tcoords

    def makekite(self):

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

        self.makesquare() # Make initial square.
        self.split()


    def sort(self):
        '''
        Time orders the co-ordinates.
        Returns a list of co-ordinates in [x,t] form i.e. list of co-ordinate pairs.
        '''
        coords = []
        for i, s in enumerate(self.xcoords):
            x = s
            t = self.tcoords[i]
            coords.append([x,t])
        self.ordered = sorted(coords, key=getkey)
        
    def makecausalmatrix(self):
        self.sort()
        ordered = self.ordered 
        G = self.gradient
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
                        
        self.causalmatrix = C

    def iSJ(self):
        '''
        Get iSJ function from Causal matrix
        '''
        # Massless scalar field
        R = 0.5*self.causalmatrix
        Rt = np.transpose(R)
        self.SJ = 1j*(R - Rt)

    def eigenvalue(self):

        '''
        Get list of eigenvectors and eigenvalues of SJ function
        
        '''

        SJ = self.SJ

        evalues, rv = linalg.eig(SJ) # Determine eigenvectors and eignenvalues of iSJ function
        evectors = []
        # DO THIS FUCKING THING
        for i, value in enumerate(evalues):
            evectors.append(rv[:,i])

        self.eigenvalues = evalues
        self.eigenvectors = evectors
        
    def SJWhight(self):

        '''
        Computes SJ Whightman function from given eigenvalues and eigenvectors.
        Use output from function(eigenvalue) as input. The eigenvectors are normalized within the function.
        '''
        evectors = self.eigenvectors
        values = self.eigenvalues.real
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

        self.Q = Q

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

    def Entropy(self):
        '''
        Take causal matrix and return entropy of space-time as per Sorkin's paper.
        '''
        W = self.Q
        self.iSJ()
        sjmask, Wmask = self.zeroprojection()

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
                self.entropy += ent 
            
        
    def zeroprojection(self):
        '''

        Function to constrain image of SJ i.e. projects out any vectors with zero eigenvalue.
        Changes representation basis to that of eigenvectors of SJ.

        '''
        sj = self.SJ
        W = self.Q
        evalues = self.eigenvalues
        evectors = self.eigenvectors
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
            tup = (evectors[i],e)
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
    xsize, tsize, sprinklesize, gradient = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    spacetime = CausalSet(xsize, tsize, sprinklesize, gradient)
    spacetime.makekite()
    spacetime.makecausalmatrix()
    spacetime.iSJ()
    spacetime.eigenvalue()
    spacetime.SJWhight()
    spacetime.Entropy()
    print 'Entropy:', spacetime.entropy
    tend = time.time() - tstart
    print 'Calculation time:', tend, 'seconds'
