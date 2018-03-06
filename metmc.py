"""
1d MMC program.

A Metropolis MC program to calculate averages and probability
distributions of an arbitrary potential. To add a new potential
define it as a function and add it to the Else-If list in
the PotCalc function.
"""

import numpy as np
import random as ran


def harmOsc(x, ind):
    k = 10.  # kcal/mol/m^2
    v = k*x[ind]**2/2.
    return v


def Morse(x, ind):
    De = 1.0  # depth of the potential well (dissociation E)
    a = 1.0   # molecular parameter for potential minimum position
    v = De*(1.0-np.exp(-a*x[ind]))**2
    return v


def Pot3(x, ind):
    a = np.sqrt(48.)  # kcal/mol/A^2
    b = 1.  # kcal/mol/A^4
    v = -0.5*a*x[ind]**2 + 0.25*b*x[ind]**4 + 0.25*a**2/b
    return v


def UmPot(x, ind):
    c = -12  # gaussian height
    d = 0.35  # width scaling (2*sigma^2)^-1
    v = c*np.exp(-d*x[ind]**2)
    return v


def PIMCpot(x, ind):
    """Calculate the effective potential including due to the harmonic
    constraint between PIMC beads
    """
    omega = 1.
    mass = 1.
    pes = PotCalc(Pot3, x, ind)
    energy = 0.5*mass*omega**2*((x[ind-1]-x[ind])**2+(x[ind+1]-x[ind])**2)+pes
    return energy


def PotCalc(pot, x, ind):
    if pot == 'HO':  # Harmonic oscillator potential
        v = harmOsc(x, ind)
    elif pot == 'Morse':  # Morse potential
        v = Morse(x, ind)
    elif pot == 'Pot3':  # quartic potential
        v = Pot3(x, ind)
    elif pot == 'UmPot':  # quartic potential + umbrella potential
        v = UmPot(x, ind) + Pot3(x, ind)
    elif pot == 'PIMC':
        v = PIMCpot(x, ind)
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
        vo.append(PotCalc(pot, xo, i))
    return xo


def Center_of_Mass(x):
    com=0
    for bead in range(nbeads):
        com+=x[bead]
    com = com/float(nbeads)
    return com


def Metropolis_Engine(xo,vo):
    xn=xo
    vn=vo
    beadind=int(ran.random()*(nbeads+1))
    xn[beadind] = xo[beadind] + move*(ran.random()-0.5)
    vn[beadind]=PotCalc(pot, xn, beadind)
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
nbeads = 10
v = False
xo = []
vo = []

while v == False:
    print('enter a potential')
    pot = input('HO,Morse,Pot3,UmPot,PIMC  ')
    v = PotCalc(pot, xo, 0)

# Generate initial coords
x = ran.random()
xo, vo = Generate_Init()

count = 0.0
samp = int(input('number of MC steps:  '))
sqrs = x**2
sqsum = sqrs
vsum = v

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
