from src.main.utils import IG_EQ_ed8 as EQ

N2O4Rxn = EQ.IG_Reaction('N2O4 Dissociation',['N2O4', 'NO2'], [-1, 2])
T = 298.15
p = 1e5
print(N2O4Rxn.DeltaH_Rxn(T))
print(N2O4Rxn.DeltaG_Rxn(T))
print(N2O4Rxn.Kp_Rxn_298())
print(N2O4Rxn.Kp_Rxn(500))
n0 = [1, 0.0]
N2O4Rxn.gibbs_plot(n0, T, p)
epsilon_EQ = N2O4Rxn.epsilon_EQ(n0, T, p)
print(N2O4Rxn.epsilon_EQ(n0, T, p))
print(N2O4Rxn.x([1, 0.0], epsilon_EQ))
print('......................')
print()

N2O4Rxn2 = EQ.IG_Reaction('N2O4 Dissociation',['N2O4', 'NO2', 'N2'], [-1, 2, 0])
print(N2O4Rxn2.DeltaH_Rxn(T))
print(N2O4Rxn2.DeltaG_Rxn(T))
print(N2O4Rxn2.Kp_Rxn_298())
print(N2O4Rxn2.Kp_Rxn(500))
n0 = [0.2, 0.0, 0.8]
N2O4Rxn2.gibbs_plot(n0, T, 1e5)
epsilon_EQ = N2O4Rxn2.epsilon_EQ(n0, T, p)
print(N2O4Rxn2.epsilon_EQ(n0, T, p))
print(N2O4Rxn2.x(n0, epsilon_EQ))

N2O4Rxn2.gibbs_plot(n0, 500, 1e5)
N2O4Rxn2.affinity_plot(n0, 500, 1e5)
