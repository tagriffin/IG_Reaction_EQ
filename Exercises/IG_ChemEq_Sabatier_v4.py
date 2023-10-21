# example Sabatier reaction
import cantera as ct
import numpy as np
from scipy.optimize import fsolve
from src.main.utils import Properties_NASA_v5 as Props
import matplotlib.pyplot as plt

# Create the ideal gas mixture with the species of interest
# Get first all of the Species objects defined in the GRI 3.0 mechanism
species = {S.name: S for S in ct.Species.listFromFile('gri30.yaml')}

# Create an IdealGas object with species representing the species being considered
considered_species = [species[S] for S in ('CO2','H2O','H2','CH4')]
gas1 = ct.Solution(thermo='IdealGas', species=considered_species)

#gas1 = ct.Solution('gri30.xml')   this would be the whole gri mech species

# Set the initial composition of the mixture
gas1.X = np.zeros(4)
gas1.X = 'CO2:1, H2:4'

gas1.TP = 300, 1e5 #K and Pa 
gas1.equilibrate('TP')

# Setting up the species for use in the properties determination

CO2 = Props.NASA_Species('CO2')
H2O = Props.NASA_Species('H2O')
H2 = Props.NASA_Species('H2')
CH4 = Props.NASA_Species('CH4')

print(gas1())

T_eq = np.zeros(100)
P_eq = np.zeros(5)
gas1_CO2 = np.zeros([100,5])
gas1_H2O = np.zeros([100,5])
gas1_H2 = np.zeros([100,5])
gas1_CH4 = np.zeros([100,5])

X_rxn = np.zeros([100,5])
X_rxn_Solve = np.zeros([100,5])
Heat_rxn = np.zeros([100,5])
Gibbs_rxn = np.zeros([100,5])

for i in range (100):
    T_eq[i] = 300 + 1500/99*i
    for j in range(5):
        if j==0:
            P_eq[0] = 1e5
        else:
            P_eq[j] = 1e5*10*j
        gas1.TP = T_eq[i],P_eq[j] #K und Pa
        gas1.equilibrate('TP')
        gas1_CO2[i,j] = gas1.X[0]
        gas1_H2O[i,j] = gas1.X[1]
        gas1_H2[i,j] = gas1.X[2]
        gas1_CH4[i,j] = gas1.X[3]
        X_rxn[i,j] = 5 * gas1.X[3] /(1 + 2 * gas1.X[3])
        Heat_rxn[i,j] = CH4.h_mol(T_eq[i],P_eq[0]) + 2*H2O.h_mol(T_eq[i],P_eq[0]) - \
        CO2.h_mol(T_eq[i],P_eq[0]) - 4*H2.h_mol(T_eq[i],P_eq[0])
        Gibbs_rxn[i,j] = CH4.g_mol(T_eq[i],P_eq[0]) + 2*H2O.g_mol(T_eq[i],P_eq[0]) - \
        CO2.g_mol(T_eq[i],P_eq[0]) - 4*H2.g_mol(T_eq[i],P_eq[0])
        eq_fcn = lambda X: np.exp(-Gibbs_rxn[i,j]*1e3/(8.314 * T_eq[i])) - 4*X**3 * (5-2*X)/((1-X) * (4-4*X)**4)
        X_guess = 0.9
        X_rxn_Solve[i,j] = fsolve(eq_fcn, X_guess)

for j in range (5):
    plt.plot(T_eq, gas1_CH4[:,j], label='CH4 @ %g bar' %(P_eq[j]/1e5))
    #plt.plot(T_eq, gas1_CO2[:,j], ':', label='CO2 @ %g bar' %(P_eq[j]/1e5))
    plt.plot(T_eq, X_rxn[:,j], ':', label='X_rxn @ %g bar' %(P_eq[j]/1e5))
    plt.plot(T_eq, X_rxn_Solve[:,j], ':', label='X_rxn_solve @ %g bar' %(P_eq[j]/1e5))
    
plt.grid(True)
plt.legend()
plt.xlabel('Temperature [K]')
plt.ylabel('Mole fraction')
plt.show(block = 'false')

for j in range (5):
    plt.plot(T_eq, Heat_rxn[:,j], label='Heat of reaction' %(P_eq[j]/1e5))
    plt.plot(T_eq, Gibbs_rxn[:,j], label='Gibbs Energy of reaction' %(P_eq[j]/1e5))
    
plt.grid(True)
plt.legend()
plt.xlabel('Temperature [K]')
plt.ylabel('Heat of Reaction kJ/mole')
plt.show(block = 'false')
