import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy import optimize
import seaborn as sns

T = np.array([298.15, 600., 800., 1000., 1200.])
G_CaO = np.array([-646.3, -663.2, -679.0,-697.4, -718.0])           #kJ/mol
G_CO2 = np.array([-457.3, -526.6, -576.7, -629.4, -684.3])          #kJ/mol
G_CaCO3 = np.array([-1234.9, -1273.7, -1309.2, -1350.6, -1397.6])   #kJ/mol

R = 8.314e-3 # kJ/(mol K)

def mu_CO2(T_value, x_CO2):
    f = interp1d(T,G_CO2)
    G_CO2_value = f(T_value)
    return G_CO2_value + R*T_value*np.log(x_CO2)

def G_CO2_fit(T_value):
    g1 = interp1d(T,G_CO2)
    return g1(T_value)

def G_CaO_fit(T_value):
    g2 = interp1d(T,G_CaO)
    return g2(T_value)

def G_CaCO3_fit(T_value):
    g3 = interp1d(T,G_CaCO3)
    return g3(T_value)

def DeltaG_rxn(T_value):
    return G_CaO_fit(T_value) + G_CO2_fit(T_value) - G_CaCO3_fit(T_value)

def A_rxn(T_value,x_value):
    return -(DeltaG_rxn(T_value) + R*T_value*np.log(x_value))

def T_Gl(x_value):
    def A_rxn_T(T_Gl):
        x = x_value
        return A_rxn(T_Gl,x)
    return optimize.fsolve(A_rxn_T,300)

def x_CO2_Gl(T_value):
    return np.exp(-DeltaG_rxn(T_value)/(R*T_value))

print('Gleichgewichtsmolanteil CO2 bei 300 K ist %.4f' %x_CO2_Gl(300))
print('Temperatur für Gleichgewicht bei x_CO2 = 0.06 ist %.4f K' %T_Gl(0.06))
print('Temperatur für Gleichgewicht bei x_CO2 = 0.10 ist %.4f K' %T_Gl(0.10))
print('Temperatur für Gleichgewicht bei x_CO2 = 0.99 ist %.4f K' %T_Gl(0.99))

x_plot = np.linspace(0.01,1,100)
sns.set_context('paper')
total = plt.figure()

diag1 = total.add_subplot(1,2,1)
diag1.plot(T,G_CO2_fit(T), label = 'G_CO2')
diag1.plot(T,G_CaO_fit(T), label = 'G_CaO')
diag1.plot(T,G_CaCO3_fit(T), label = 'G_CaCO3')
diag1.plot(T,G_CO2_fit(T)+G_CaO_fit(T)-G_CaCO3_fit(T), label = 'G_rxn')
diag1.scatter(T,G_CO2)
diag1.scatter(T,G_CaO)
diag1.scatter(T,G_CaCO3)
diag1.set_xlabel('Temperatur, K')
diag1.set_ylabel('Molar Gibbs Energie, kJ/mol')
diag1.legend(loc='upper left', bbox_to_anchor=(0.15,0.85))
diag1.grid('on')

diag2 = total.add_subplot(1,2,2)
diag2.plot(x_plot,A_rxn(900,x_plot), label = 'A_rxn (900 K)')
diag2.plot(x_plot,A_rxn(1000,x_plot), label = 'A_rxn (1000 K)')
diag2.plot(x_plot,A_rxn(1100,x_plot), label = 'A_rxn (1100 K)')
diag2.plot(x_plot,A_rxn(1150,x_plot), label = 'A_rxn (1150 K)')
diag2.plot(x_plot,np.zeros(100), 'r:')

diag2.set_xlabel('Molanteil CO2, x')
diag2.set_ylabel('Affinität, kJ/mol')
diag2.legend(loc='best')
diag2.grid('off')

total.tight_layout()
plt.savefig('Gibbs_Aufgabe_13-7')

total.show()
