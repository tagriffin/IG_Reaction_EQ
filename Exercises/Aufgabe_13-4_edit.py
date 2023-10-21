import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate

from src.main.utils import IG_Props_Comb as Props

H2 = Props.Molecule('H2')
H2O = Props.Molecule('H2O')
CO = Props.Molecule('CO')
CO2 = Props.Molecule('CO2high')
CH4 = Props.Molecule('CH4')

# Absolute ("konventionelle Enthalpien") bei 500 und 1000 K
T1 = 500
T2 = 1100

table = [[T1, H2.h_abs(T1),H2O.h_abs(T1),CO.h_abs(T1), CO2.h_abs(T1), CH4.h_abs(T1)], \
         [T2, H2.h_abs(T2),H2O.h_abs(T2),CO.h_abs(T2), CO2.h_abs(T2), CH4.h_abs(T2)]]

print(tabulate(table,headers = ['T(K)', 'H2','H2O', 'CO', 'CO2', 'CH4'])) \
                             # tabulate only works on table with two rows


T = np.linspace(298,1100,100)
plt.plot(T, H2.h_abs(T), label = 'h_H2')
plt.plot(T, H2O.h_abs(T), label = 'h_H2O')
plt.plot(T, CO.h_abs(T), label = 'h_CO')
plt.plot(T, CO2.h_abs(T), label = 'h_CO2')
plt.plot(T, CH4.h_abs(T), label = 'h_CH4')
plt.xlabel('Temperatur, K')
plt.ylabel('Absolute Enthalpie, kJ/mol')
plt.legend(loc='best')
plt.show()
