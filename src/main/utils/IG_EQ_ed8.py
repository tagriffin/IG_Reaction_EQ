'''This file creates an object for the calculation of reaction equilibrium
with an ideal gas mixture'''

import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

import cantera as ct
from src.main.utils import Properties_NASA_v5 as NASA
import CoolProp.CoolProp as CP

# Class for pure or mixture fuels, type 1
class IG_Reaction(object):
    def __init__(self, name, Components, StoichCoeff): # NASA_name, Hu MJ per kg dry fuel
        """
        Pass parameters describing molecules
        """
        self.name = name
        self.comps = Components
        self.coeffs = StoichCoeff
        self.num_Comps = len(Components)
        self.coeff_sum = np.sum(self.coeffs)
        
    def DeltaH_Rxn(self, Temp): # returns Heat of Reaction in kJ/mole
        DelH = 0
        Press = 1e5
        for i in range(self.num_Comps):
            DelH += NASA.NASA_Species(self.comps[i]).h_mol(Temp, Press) * self.coeffs[i]
        return DelH

    def DeltaG_Rxn(self, Temp): # returns Heat of Reaction in kJ/mole
        DelG = 0
        Press = 1e5
        for i in range(self.num_Comps):
            DelG += NASA.NASA_Species(self.comps[i]).g_mol(Temp, Press) * self.coeffs[i]
        return DelG

    def Kp_Rxn_298(self): # returns Equilibrium Constant Kp at 25 °C
        DelG = 0
        Temp = 298.15
        Press = 1e5
        Rm = 8.314 # J/(kg mole)
        for i in range(self.num_Comps):
            DelG += NASA.NASA_Species(self.comps[i]).g_mol(Temp, Press) * self.coeffs[i]
        return np.exp(-(DelG * 1000)/(Rm * Temp))

    def Kp_Rxn(self, Temp): # returns Equilibrium Constant Kp at Temperature
        DelG = 0
        Press = 1e5
        Rm = 8.314 # J/(kg mole)
        for i in range(self.num_Comps):
            DelG += NASA.NASA_Species(self.comps[i]).g_mol(Temp, Press) * self.coeffs[i]
        return np.exp(-(DelG * 1000)/(Rm * Temp))

    def z(self, n0): # returns z_min, z_max, given starting composition in moles
        z_min = -np.max(n0)
        z_max = np.max(n0)
        n = np.zeros(self.num_Comps)
        molfcn = np.zeros(self.num_Comps)
        for i in range(self.num_Comps):
            if self.coeffs[i] < 0:
                if -n0[i]/self.coeffs[i] < z_max:
                    z_max = -n0[i]/self.coeffs[i]
            if self.coeffs[i] > 0:
                if -n0[i]/self.coeffs[i] > z_min:
                    z_min = -n0[i]/self.coeffs[i]
        return z_min, z_max

    def x(self, n0, epsilon): # returns mole fractions as function of epsilon, given starting composition in moles
        z_min = self.z(n0)[0]
        z_max = self.z(n0)[1]
        n = np.zeros(self.num_Comps)
        molfcn = np.zeros(self.num_Comps)
        z = z_min + epsilon * (z_max - z_min)
        
        for i in range(self.num_Comps):
            n[i] = n0[i] + z * self.coeffs[i]

        molsum = np.sum(n)
        for i in range(self.num_Comps):
            molfcn[i] = (n0[i] + z * self.coeffs[i])/molsum
        return molfcn

    def epsilon_EQ(self, n0, Temp, press):
        def EQ_fcn(epsilon):
            x_term = np.zeros(self.num_Comps)
            for i in range(self.num_Comps):
                x_term[i] = self.x(n0, epsilon)[i]**self.coeffs[i]
            return self.Kp_Rxn(Temp) - np.prod(x_term) * (press/1e5)**self.coeff_sum
        epsilon = 0.9 # initial guess
        return optimize.fsolve(EQ_fcn, epsilon)

    def x_epsilon_plot(self, n0):
        epsilon = np.linspace(0,1,100)
        x_values = np.zeros([self.num_Comps, 100])
        
        for i in range(0,len(epsilon)):
            for j in range(self.num_Comps):
                x_values[j,i] = self.x(n0, epsilon[i])[j]

        for j in range(self.num_Comps):
            plt.plot(epsilon,x_values[j, :],label = self.comps[j])
    
        plt.title('Mole Fractions as Function of epsilon for ' + self.name)
        plt.xlabel('epsilon')
        plt.legend(loc='best')
        plt.grid()
        plt.show()

    def affinity_plot(self, n0, T, p):  # plot the affinity of the mixture as a function of epsilon for T, p
        epsilon = np.linspace(0,1,100)
        x_values = np.zeros([self.num_Comps, 100])
        affinity = np.zeros(100)
        Rm = 8.314 # J/(mol K)
        for i in range(0,len(epsilon)):
            for j in range(self.num_Comps):
                x_values[j,i] = self.x(n0, epsilon[i])[j]
                affinity[i] += - self.coeffs[j]*(NASA.NASA_Species(self.comps[j]).g_mol(T, p)*1000 + Rm*T*np.log(x_values[j,i]))
                                                       
        for j in range(self.num_Comps):
            plt.plot(epsilon,affinity)
    
        plt.title('Affinity as Function of epsilon for ' + self.name)
        plt.xlabel('epsilon')
        plt.grid()
        plt.show()

    def gibbs_plot(self, n0, T, p):  # plot the affinity of the mixture as a function of epsilon for T, p
        epsilon = np.linspace(0.001,0.999,100)
        n_values = np.zeros([self.num_Comps, 100])
        x_values = np.zeros([self.num_Comps, 100])
        z_min = self.z(n0)[0]
        z_max = self.z(n0)[1]     
        gibbs = np.zeros(100)
        Rm = 8.314 # J/(mol K)
        for i in range(0,len(epsilon)):
            z = z_min + epsilon[i] * (z_max - z_min)
            for j in range(self.num_Comps):
                n_values[j,i] = n0[j] + z * self.coeffs[j]
                x_values[j,i] = self.x(n0, epsilon[i])[j]
                gibbs[i] += n_values[j,i]*(NASA.NASA_Species(self.comps[j]).g_mol(T, p)*1000 + Rm*T*np.log(x_values[j,i]))
                                                       
        for j in range(self.num_Comps):
            plt.plot(epsilon,gibbs)
           
        plt.title('Gibbs Energy (J/mol) as Function of epsilon for ' + self.name)
        plt.xlabel('epsilon')
        plt.grid()
        plt.show()  
        
                                          
                                          

