"""
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
    d = 0.35  #  width scaling (2*sigma^2)^-1 
    v = c*np.exp(-d*x**2)
    return v

def PotCalc(pot):
    if pot == 'HO':  #  Harmonic oscillator potential
        v = harmOsc(x)
    elif pot == 'Morse':  #  Morse potential
        v = Morse(x)
    elif pot == 'Pot3':  #  quartic potential
        v = Pot3(x)
    elif pot == 'UmPot':  #  quartic potential + umbrella potential
        v = UmPot(x) + Pot3(x)
    else:
        v = False
    return v

x = np.sqrt(np.sqrt(48))
#x = ran.random()  # Angstroms
hist = np.zeros([10000])
binsize = 0.001  # A
move = .5  # A
xsum = x
beta = 1.0/0.6  #  mol/kcal
v = False
while v == False:
    print('enter a potential')
    pot = input('HO,Morse,Pot3,UmPot  ')
    v = PotCalc(pot)
    
count = 0.0
samp = int(input('number of MC steps:  '))
sqrs = x**2
sqsum = sqrs
vsum = v

#xfile = open('xavg.txt','w')
#sqrfile = open('sqravg.txt','w')
histfile = open('umb_hist.txt','w')

for i in range(samp):
    xold = x
    vold = v
    x += move*(ran.random()-0.5)
    v = PotCalc(pot)
    if v > vold:
        prob = np.exp(-beta*(v-vold))
        if prob < ran.random():
            x = xold
            v = vold
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
        histfile.write(str(x) + ' ' + str(p/(count*binsize)) + '\n' )
