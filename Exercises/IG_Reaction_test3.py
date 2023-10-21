from src.main.utils import IG_EQ_ed8 as EQ
import matplotlib.pyplot as plt
import numpy as np

N2O4Rxn = EQ.IG_Reaction('N2O4 Dissociation',['N2O4', 'NO2'], [-1, 2])
print(N2O4Rxn.DeltaH_Rxn(298.15))
print(N2O4Rxn.DeltaG_Rxn(298.15))
print(N2O4Rxn.Kp_Rxn_298())
print(N2O4Rxn.Kp_Rxn(500))

epsilon_EQ = N2O4Rxn.epsilon_EQ([1, 0.0], 298.15, 1e5)
print(N2O4Rxn.epsilon_EQ([1, 0.0], 298.15, 1e5))
print(N2O4Rxn.x([1, 0.0], epsilon_EQ))
print('......................')
print()

N2O4Rxn2 = EQ.IG_Reaction('N2O4 Dissociation',['N2O4', 'NO2', 'N2'], [-1, 2, 0])
print(N2O4Rxn2.DeltaH_Rxn(298.15))
print(N2O4Rxn2.DeltaG_Rxn(298.15))
print(N2O4Rxn2.Kp_Rxn_298())
print(N2O4Rxn2.Kp_Rxn(500))


epsilon_EQ = N2O4Rxn2.epsilon_EQ([0.2, 0.0, 0.8], 298.15, 1e5)
print(N2O4Rxn2.epsilon_EQ([0.2, 0.0, 0.8], 298.15, 1e5))
print(N2O4Rxn2.x([0.2, 0.0, 0.8], epsilon_EQ))
N2O4Rxn.affinity_plot([0.2, 0.0, 0.8], 1000,1e5)
N2O4Rxn.gibbs_plot([0.2, 0.0, 0.8], 1000,1e5)
