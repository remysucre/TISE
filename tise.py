import numpy as np
from scipy.integrate import odeint
import scipy.optimize as opt
import matplotlib.pyplot as mpl

# horizontal axis in the plots of psi
x0 = 0
x1 = 10
numX = 100
dx = (x1 - x0)/numX
x = np.linspace(x0, x1, numX)

# guesses of E
E0 = 1.5
E1 = 10.5
numE = 10
Es = np.linspace(E0,E1,numE)

# increasing L in psi(L) = 0
L0 = 7
L1 = 30
numL = 10
lengths = np.linspace(L0,L1,numL)

# inital psi and psi'
y0 = [0.0, 1.0]

# time-independent Schroedinger's equation with semi-infinite V
# given y = [psi, psi'] returns [psi', psi'']
def tise(y, x, E):
    psi, u = y # u = d/dx*psi
    dydx = [u, 2*(x*psi - E*psi)] # right-hand side of ODE
    return dydx


# takes array of psi(l), returns all indices just before psi(l) crosses zero
def findSignChange(f):
    changes = [];
    for i in range(0, len(f)-1):
        if(f[i+1]*f[i]<=0):
            changes.append(i)
    return changes

#returns function psi(E) at x = l
def getPsiFcn(l):
    #psi at length l
    def psi_l(E):
        x = np.linspace(x0,l,numX)
        psi = odeint(tise, y0, x, args=(E,))
        return psi[-1,0]
    return psi_l

# find the sets of eigenvalues for lengths ls
# with guesses Es
def findEigenVals(ls, Es, getPsiFcn):
    rootConvergence = [] # sets of eigenvalues

    for l in ls:
        psiFcn = getPsiFcn(l) # psi(l)
        psi_l = [] #sample values of psi_l(E)

        for E in Es:
            psi_l.append(psiFcn(E)) # shooting E

        guesses = findSignChange(psi_l) # brackets for root finder
        roots = [] # roots for current l

        # find roots for each bracket
        for guess in guesses:
            lower = Es[guess]
            upper = Es[guess+1]
            root = opt.bisect(psiFcn ,lower, upper)
            roots.append(root)

        rootConvergence.append(roots)
    return rootConvergence

# normalize the eigenfunction f stored in an array
def normalize(f,dx):
    integral= 0 # integrate f^2 over x

    for i in f:
        integral += i**2

    integral = integral*(dx)
    return f/np.sqrt(integral)




# sets of eigenvalues for different L in psi(L) = 0
eigenMat = findEigenVals(lengths,Es,getPsiFcn)

mpl.clf()

# plot the sets of eigenvalues obtained for
# each L in psi(L) = 0
for i in range(0,len(lengths)):
    for r in eigenMat[i]:
        mpl.scatter(lengths[i], r)

mpl.show()

# converging eigenvalues
eigenVals = eigenMat[-1]
print(eigenVals)

# number of eigenfunctions to show
eigenPlotNum = 5

# plot the potential V = x
mpl.plot(x, x, 'k', label='V')

# plot the first several converging eigenfunctions
for i in range(eigenPlotNum):
        psi = odeint(tise, y0, x, args=(eigenVals[i],))
        psi_norm = normalize(psi, dx)
        mpl.plot(x, psi_norm[:, 0], label='$\Psi_'+str(i)+'$')

mpl.xlim(0, x1)
mpl.xlabel('x')

mpl.ylim(-1, x1)
mpl.ylabel(r'$\Psi$')

mpl.gca().set_aspect('equal', adjustable='box')
mpl.grid()
mpl.legend()

mpl.savefig('eigen-functions.pdf')
# mpl.show()
