from src.main.utils import IG_EQ_ed8 as EQ
import matplotlib.pyplot as plt
import numpy as np

CO_ox = EQ.IG_Reaction('CO2 Oxidation',['CO2', 'CO', 'O2'], [1, -1, -0.5])
print(CO_ox.DeltaH_Rxn(298.15))
print(CO_ox.DeltaG_Rxn(298.15))
print(CO_ox.Kp_Rxn_298())
print(CO_ox.Kp_Rxn(2000))
print(CO_ox.Kp_Rxn(2500))

epsilon_EQ = CO_ox.epsilon_EQ([0, 1, 0.5], 2500, 1e5)
print(CO_ox.epsilon_EQ([0, 1, 0.5], 2500, 1e5))
print(CO_ox.x([0, 1, 0.5], epsilon_EQ))
print('......................')
print()

epsilon_EQ = CO_ox.epsilon_EQ([0, 1, 0.5], 2500, 10e5)
print(CO_ox.epsilon_EQ([0, 1, 0.5], 2500, 10e5))
print(CO_ox.x([0, 1, 0.5], epsilon_EQ))
print('......................')
print()
CO_ox.affinity_plot([0, 1, 0.5], 2500, 10e5)
