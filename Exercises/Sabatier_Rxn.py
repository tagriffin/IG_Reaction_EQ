from src.main.utils import IG_EQ_ed8 as EQ

Sabatier1 = EQ.IG_Reaction('CO2 + 4H2 = CH4 + 2H2O',['CO2', 'H2', 'CH4', 'H2O'], [-1, -4, 1, 2])
Sabatier2 = EQ.IG_Reaction('CO + 3H2 = CH4 + H2O',['CO', 'H2', 'CH4', 'H2O'], [-1, -3, 1, 1])

print(Sabatier1.DeltaH_Rxn(298.15))
print(Sabatier1.DeltaG_Rxn(298.15))
print(Sabatier1.Kp_Rxn_298())

print(Sabatier2.DeltaH_Rxn(298.15))
print(Sabatier2.DeltaG_Rxn(298.15))
print(Sabatier2.Kp_Rxn_298())

n0 = [1, 4, 0, 0 ]

Sabatier1.x_epsilon_plot(n0)


epsilon1EQ_500 = Sabatier1.epsilon_EQ(n0, 500 + 273.15, 1e5)
molfcn = Sabatier1.x(n0, epsilon1EQ_500)
print(molfcn[2], molfcn[3])
