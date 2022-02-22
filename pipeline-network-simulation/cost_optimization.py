import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

# function to calculate pressure drop with Churchill eqn
def pdrop(q,l,d):
    re = (4 * q * density) / (np.pi * viscosity * d)
    a = (2.457 * np.log(1/((7/re)**0.9 + 0.27 * eps / d))) ** 16
    b = (37530/re)**16
    f = 2 * ((8/re)**12 + 1/((a+b)**1.5))**(1/12)
    return 2 * f * density * ((4*q)/(np.pi*d**2))**2 * l / d

# function with mass balances and energy balances
def equations(q,d):
    return [q[0]-q[1]-q[2],
            q[1]-q[3]-q[4],
            q[3]-q[5]-endflow,
            q[2]+q[4]-q[6]-q[7],
            q[5]+q[6]-q[8]-endflow,
            q[7]-q[9]-endflow,
            q[8]+q[9]-endflow,
            pdrop(q[1],lengths[1],d) + pdrop(q[4],lengths[4],d) - pdrop(q[2],lengths[2],d),
            pdrop(q[3],lengths[3],d) + pdrop(q[5],lengths[5],d) - pdrop(q[4],lengths[4],d) - pdrop(q[6],lengths[6],d),
            pdrop(q[6],lengths[6],d) + pdrop(q[8],lengths[8],d) - pdrop(q[7],lengths[7],d) - pdrop(q[9],lengths[9],d)]

# initialize variables
density = 1113 * 0.65 + 998 * 0.35 # kg/m^3
endflow = 72 / (60 * 1000) # kg/s at 3,5,6,7
viscosity = 0.0161 * 0.65 + 0.0010 * 0.35 # Pa s
eps = 5E-5 # metres
lengths = np.array([230,460,585,460,380,380,460,380,380,460]) # metres

# initial guesses
q0 = [endflow * 4]
for i in range(1,10):
    if i == 1 or i == 2:
        q0.append(q0[0]/2)
    elif i == 3 or i == 4:
        q0.append(q0[1]/2)
    elif i == 5:
        q0.append(q0[3])
    elif i == 6 or i == 7:
        q0.append((q0[2]+q0[4])/2)
    elif i == 8:
        q0.append(q0[5]+q0[6])
    else:
        q0.append(q0[7])

# initialize variables for cost analysis
diameters = np.linspace(1,4,13) * 0.0254
tac = []
aoc = []
acc = []
# function to calculate aoc
f_acc = lambda d : (1+1)*(5.92/0.3048)*((d/0.0254)**1.25)*0.24*sum(lengths)


# calculate acc, aoc, and tac for different diameters 
for i in range(len(diameters)):
    sum_term = 0
    q = fsolve(equations,q0,args=(diameters[i]))
    for j in range(10):
        sum_term += (q[j] * pdrop(q[j],lengths[j],diameters[i]))
    aoc.append(((sum_term / 1000) / (0.65 * 0.80)) * 8400 * 0.105)
    acc.append(f_acc(diameters[i]))
    tac.append(aoc[i] + acc[i])


print(f'Minimum TAC is ${round(min(tac),2)}/year. Optimum diameter is {diameters[tac.index(min(tac))]} metres')

# plotting
fig, axs = plt.subplots(1, 3)
axs[0].plot(diameters, aoc,'-')
axs[0].set_title('AOC vs. Diameter')
axs[0].set(xlabel='diameter (m)', ylabel='$/year')
axs[1].plot(diameters, acc,'-')
axs[1].set_title('ACC vs. Diameter')
axs[1].set(xlabel='diameter (m)', ylabel='$/year')
axs[2].plot(diameters, tac,'-')
axs[2].set_title('TAC vs. Diameter')
axs[2].set(xlabel='diameter (m)', ylabel='$/year')

plt.show()