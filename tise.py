import numpy as np
from scipy.integrate import odeint
import scipy.optimize as opt
import matplotlib.pyplot as mpl
from mpl_toolkits.mplot3d import Axes3D



y0 = [0.0, 1.0] # inital psi and psi'

# time-independent schroedinger's equation with semi-infinite V
# given y = [psi, psi'] returns [psi', psi'']
def tise(y, x, E):
	psi, u = y # u is d/dx*psi
	dydx = [u, 2*(x*psi - E*psi)]
	return dydx


# takes array of psi(l), returns all indices just before psi(l) crosses zero
def findSignChange(f):
	changes = [];
	for i in range(0, len(f)-1):
		if(f[i+1]*f[i]<=0):
			changes.append(i)
	return changes

#get function psi(E) with the max length l
def getPsiFcn(l):
	#returns psi at length l
	def psi_l(E):
		x = np.linspace(0,l,101)
		psi = odeint(tise, y0, x, args=(E,))
		return psi[-1,0]
	return psi_l

def findEigenVals(lengths, Es, getPsiFcn):
	rootConvergence = []
	for length in lengths:

		psiFcn = getPsiFcn(length)

		#sample values of psi_l(E)
		psi_l = []
		for E in Es:
			psi_l.append(psiFcn(E))
		guesses = findSignChange(psi_l)
		roots = []
		for guess in guesses:
		    #root = opt.newton(f, E[guesses[k]])
			upper = Es[guess+1]
			lower = Es[guess]
			root = opt.bisect(psiFcn ,lower, upper)
			roots.append(root)
		rootConvergence.append(roots)
	return rootConvergence

def normalize(f,dx):
	integral= 0
	for i in f:
		integral += i**2
		
	integral = integral*(dx)
	return integral


numE = 20
numL = 50

E0 = 1.5
E1 = 10.5
L0 = .1
L1 = 30

Es = np.linspace(E0,E1,numE)
lengths = np.linspace(L0,L1,numL)

eigenMat = findEigenVals(lengths,Es,getPsiFcn)
mpl.clf()
for i in range(0,len(lengths)):
	for r in eigenMat[i]:
		mpl.scatter(lengths[i], r)
mpl.show()
"""
mpl.rcParams['legend.fontsize'] = 10
fig = mpl.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(guesses, roots, lengths, label='parametric curve')
mpl.plot(rootConvergence.slice())
mpl.show()
"""
"""
x = np.linspace(L0,L1,numL)
s = odeint(tise, y0, x, args=(1.85575704051,))

norm = normalize(s[:,0]	,(L1-L0)/numL)
mpl.clf()
mpl.plot(x,s[:,0]*norm)
mpl.show()
"""
