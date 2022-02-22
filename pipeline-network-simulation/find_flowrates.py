import numpy as np
from scipy.optimize import fsolve
import pandas as pd
from tabulate import tabulate

# function to calculate pressure drop with Churchill eqn
def pdrop(q,l):
    re = (4 * q * density) / (np.pi * viscosity * diameter)
    a = (2.457 * np.log(1/((7/re)**0.9 + 0.27 * eps / diameter))) ** 16
    b = (37530/re)**16
    f = 2 * ((8/re)**12 + 1/((a+b)**1.5))**(1/12)
    return 2 * f * density * ((4*q)/(np.pi*diameter**2))**2 * l / diameter

# function with mass balances and energy balances
def equations(q):
    return [q[0]-q[1]-q[2],
            q[1]-q[3]-q[4],
            q[3]-q[5]-endflow,
            q[2]+q[4]-q[6]-q[7],
            q[5]+q[6]-q[8]-endflow,
            q[7]-q[9]-endflow,
            q[8]+q[9]-endflow,
            pdrop(q[1],lengths[1]) + pdrop(q[4],lengths[4]) - pdrop(q[2],lengths[2]),
            pdrop(q[3],lengths[3]) + pdrop(q[5],lengths[5]) - pdrop(q[4],lengths[4]) - pdrop(q[6],lengths[6]),
            pdrop(q[6],lengths[6]) + pdrop(q[8],lengths[8]) - pdrop(q[7],lengths[7]) - pdrop(q[9],lengths[9])]

# initialize variables
diameter = 0.0254 # metres
density = 1113 * 0.65 + 998 * 0.35 # kg/m^3
endflow = 72 / (60 * 1000) # m^3/s at 3,5,6,7
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

# solve for q
x = fsolve(equations,q0)

# append results into a pandas dataframe
if __name__ =="__main__":
    pressures = [pdrop(x[i],lengths[i]) for i in range(10)]
    labels = ['Pipe 01', 'Pipe 12', 'Pipe 14', 'Pipe 23','Pipe 24', 'Pipe 35', 'Pipe 45', 'Pipe 46','Pipe 57', 'Pipe 67']
    df = pd.DataFrame(labels,columns=['Pipe Section'])
    df['Flowrate (m^3/s)'] = x
    df['Pressure Drop (Pa)'] = pressures
    print(tabulate(df, headers='keys', tablefmt='fancy_grid'))