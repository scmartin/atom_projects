"""
1d MMC program.

A Metropolis MC program to calculate averages and probability
distributions of an arbitrary potential. To add a new potential
define it as a function and add it to the Else-If list in
the PotCalc function.
"""

import numpy as np
import random as ran


def harmOsc(x):
    global k
    k = 10.  # kcal/mol/m^2
    v = k*x**2/2.
    return v


def Morse(x):
    De = 1.0  # depth of the potential well (dissociation E)
    a = 1.0   # molecular parameter for potential minimum position
    v = De*(1.0-np.exp(-a*x))**2
    return v


def Pot3(x):
    a = np.sqrt(48.)  # kcal/mol/A^2
    b = 1.  # kcal/mol/A^4
    v = -0.5*a*x**2 + 0.25*b*x**4 + 0.25*a**2/b
    return v


def UmPot(x):
    c = -12  # gaussian height
    d = 0.35  # width scaling (2*sigma^2)^-1
    v = c*np.exp(-d*x**2)
    return v


def PIMCpot(x, ind):
    """Calculate the effective potential including due to the harmonic
    constraint between PIMC beads
    """
    omega = 1.
    mass = 1.
    pes = PotCalc(pot, x[ind])
    if ind < nbeads-1:
        energy = 0.5*mass*omega**2*((x[ind-1]-x[ind])**2+(x[ind+1]-x[ind])**2)+pes
    else:
        energy = 0.5*mass*omega**2*((x[ind-1]-x[ind])**2+(x[0]-x[ind])**2)+pes
    return energy


def PotCalc(pot, x):
    if pot == 'HO':  # Harmonic oscillator potential
        v = harmOsc(x)
    elif pot == 'Morse':  # Morse potential
        v = Morse(x)
    elif pot == 'Pot3':  # quartic potential
        v = Pot3(x)
    elif pot == 'UmPot':  # quartic potential + umbrella potential
        v = UmPot(x) + Pot3(x)
    else:
        v = False
    return v


def Tot_Pot():
    v = 0
    delx = 0
    for bead in range(nbeads):
        if bead == 0:
            delx = xo[bead]-xo[nbeads]
        else:
            delx = xo[bead]-xo[bead+1]
        v+=HarmOsc(delx)
    return v


def Generate_Init():
    for bead in range(nbeads):
        if bead == 0:
            xo.append(x+move*(ran.random()-0.5))
        else:
            xo.append(xo[bead-1]+move*(ran.random()-0.5))
    for i in range(nbeads):
        vo.append(PotCalc(pot, xo[i]))
    return xo, vo


def Center_of_Mass(x):
    com=0
    for bead in range(nbeads):
        com+=x[bead]
    com = com/float(nbeads)
    return com


def Metropolis_Engine(xo,vo):
    xn=xo
    vn=vo
    beadind=int(ran.random()*(nbeads))
    print(beadind)
    xn[beadind] = xo[beadind] + move*(ran.random()-0.5)
    if pimc == 'y':
        vn[beadind] = PIMCpot(xn, beadind)
    else:
        vn[beadind] = PotCalc(pot, xn[beadind])
    print("old potentials")
    print(vo)
    print("new potentials")
    print(vn)
    print()
    if vn[beadind] > vo[beadind]:
        prob = np.exp(-beta*(vn[beadind]-vo[beadind]))
        if prob < ran.random():
            xn[beadind] = xo[beadind]
            vn[beadind] = vo[beadind]

    x = Center_of_Mass(xn)
    return x, xn, vn


hist = np.zeros([10000])
binsize = 0.001  # A
move = .5  # A
beta = 1.0/0.6  # mol/kcal
nbeads = 1
v = False
xo = []
vo = []
x = ran.random()
pimc = input("PIMC (y or n)?")
if pimc != 'y':
    pimc == 'n'
else:
    nbeads = int(input("enter number of PIMC beads"))
while v == False:
    print('enter a potential')
    pot = input('HO,Morse,Pot3,UmPot  ')
    v = PotCalc(pot, x)

# Generate initial coords
xo, vo = Generate_Init()

count = 0.0
samp = int(input('number of MC steps:  '))
sqrs = x**2
sqsum = sqrs
vsum = v
xsum = 0.0

# xfile = open('xavg.txt','w')
# sqrfile = open('sqravg.txt','w')
histfilename = str(pot)+'.txt'
histfile = open(histfilename, 'w')

for i in range(samp):
    x, xn, vn = Metropolis_Engine(xo, vo)
    count += 1.0
#    print(x)
    bin = int((x+5.)/binsize)
#    print(bin)
    hist[bin] += 1
    xsum += x
    sqsum += sqrs
    vsum += v
    xavg = xsum/count
    sqrav = sqsum/count
    vavg = vsum/count
"""    if count % 100 == 0:
        xfile.write(str(count) + ' ' + str(xavg) + '\n')
        sqrfile.write(str(count) + ' ' + str(sqrav) + '\n')
"""
xavg = xsum/count
sqrav = sqsum/count
vavg = vsum/count
print(xavg)  #  Angstroms
print(sqrav)  #  A^2
print(vavg)  #  kcal/mol
for i,p in enumerate(hist):
    x = float(i)*binsize-5.+binsize*0.5
    if pot=='UmPot':
        vum = UmPot(x)
        histfile.write(str(x) + ' ' + str(-p/(vum*559.303)) + '\n' )
        # normalized based on XMGRACE integration
    elif pot=='HO':
        p_x = np.sqrt(beta*k/np.pi)*np.exp(-beta*k*x**2)  # analytical P(x)
        histfile.write(str(x) + ' ' + str(p/(count*binsize)) +
                   ' ' + str(p_x) + '\n')
    else:
        histfile.write(str(x) + ' ' + str(p/(count*binsize)) + '\n')
