# example water-gas-shift-adv thermo week 7
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# Create the ideal gas mixture with the species of interest
# Get first all of the Species objects defined in the GRI 3.0 mechanism
species = {S.name: S for S in ct.Species.listFromFile('gri30.yaml')}

# Create an IdealGas object with species representing complete combustion
considered_species = [species[S] for S in ('CO','O2','CO2')]
gas1 = ct.Solution(thermo='IdealGas', species=considered_species)

#gas1 = ct.Solution('gri30.xml')   this would be the whole gri mech species

# Set the initial composition of the mixture
gas1.X = np.zeros(3)
gas1.X = 'CO:1, O2:0.5'


gas1.TP = 2000, 1e5 #K and Pa 
gas1.equilibrate('TP')

print(gas1())


#plt.plot(phi, T_complete, label='complete combustion', lw=2)
#plt.plot(phi, T_incomplete, label='incomplete combustion', lw=2)
#plt.grid(True)
#plt.xlabel('Equivalence ratio, $\phi$')
#plt.ylabel('Temperature [K]');
#plt.show()
