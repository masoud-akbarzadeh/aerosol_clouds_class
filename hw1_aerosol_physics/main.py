# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import matplotlib.pyplot as plt
import numpy as np

# all diameters are in nm
# Dp = np.linspace(1e-9, 1e-8, 100)

S_g = 1.7 # standard deviation of the lognormal distribution
N = 1000 # number of particles in the distribution (cm^-3)
D_pg = 200 # geometric mean diameter (nm)

Dp = np.geomspace(1e-8, 1e-6, 100)*1e9 # convert to nm
dN_dDp = N/(np.sqrt(2*np.pi)*np.log(S_g)*Dp)*np.exp(-(np.log(Dp/D_pg))**2/(2*(np.log(S_g))**2))
print (Dp, dN_dDp)
print((np.sqrt(2*np.pi)*np.log(S_g)*Dp))
print(np.exp(-np.log(Dp/D_pg)**2/(2*(np.log(S_g))**2)))
print(-np.log(Dp/D_pg)**2)
plt.plot(Dp, dN_dDp)
plt.xscale('log')
plt.legend(['dN/dDp'])
plt.xlabel('Dp (nm)')
plt.ylabel('dN/dDp (cm^-3)')
plt.show()


dN_dlnDp = dN_dDp*Dp
plt.plot(Dp, dN_dlnDp)
plt.xscale('log')
plt.legend(['dN/dlnDp'])
plt.xlabel('Dp (nm)')
plt.ylabel('dN/dlnDp')
plt.show()


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
