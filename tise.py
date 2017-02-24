import numpy as np
from scipy.integrate import odeint
import scipy.optimize as opt
import matplotlib.pyplot as mpl
from mpl_toolkits.mplot3d import Axes3D

E0 = 1.5 # constant energy
E1 = 10.5
iterations = 20

E = np.linspace(E0,E1,iterations)
y0 = [0.0, 1.0] # inital psi and psi'
l = 10
x = np.linspace(0, l, 101) # position

# time-independent schroedinger's equation with semi-infinite V
def tise(y, x, E):
    psi, u = y # u is d/dx*psi
    dydx = [u, 2*(x*psi - E*psi)]
    return dydx

# numerical solution of TISE
psi_l = np.zeros(iterations)

"""
for i in range(iterations):
    sol = odeint(tise, y0, x, args=(E[i],))
    # plot psi over psi'
    mpl.plot(x, sol[:, 0], 'b', label='$\Psi$')
    #mpl.plot(x, sol[:, 1], 'g', label="$\dot{\Psi}$")
    psi_l[i] = sol[-1,0]
"""


"""
#mpl.legend(loc='best')
mpl.xlabel('x')
mpl.ylabel('psi(x)')

#mpl.show()

mpl.clf()
mpl.plot(E,psi_l)
mpl.show()
"""


def findSignChange(f):
    changes = [];
    for i in range(0, len(f)-1):
        if(f[i+1]*f[i]<=0):
            changes.append(i)
    return changes

def getFcn(l):
    def fcn(E):
        x = np.linspace(0,l,101)
        psi = odeint(tise, y0, x, args=(E,))
        return psi[-1,0]
    return fcn


#for each length, i
lengths = np.linspace(5,30,101)
print('sign changes:')
guesses = findSignChange(psi_l)
print(E[guesses])
mpl.clf()
rootConvergence = []


mpl.rcParams['legend.fontsize'] = 10

fig = mpl.figure()
ax = fig.gca(projection='3d')


for j in range(0,101):

    f = getFcn(lengths[j])

    psi_l = np.zeros(len(E))
    for k in range(0,len(E)):
        psi_l[k] = f(E[k])
    guesses = findSignChange(psi_l)

    roots = []
    for k in range(0,len(guesses)):
        root = opt.newton(f, E[guesses[k]])
        # mpl.scatter(k,root)
        ax.scatter(k, root, lengths[j], label='parametric curve')
        roots.append(root)
    #rootConvergence.append(roots)

# ax.plot_surface(guesses, roots, lengths, label='parametric curve')

#mpl.plot(rootConvergence.slice())
mpl.show()

"""
num_l = 1
roots = np.zeros(num_l)
for i in range(num_l):
    roots[i] = bisect(getFcn(l*i/num_l), E0,E1)

mpl.clf()
mpl.plot(np.linspace(0,l,num_l), roots)
mpl.show()
"""
