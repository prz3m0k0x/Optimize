#Model Ludwika naprezen sigma_pl = sigma_y + K * epsilon ** n
#Nalezy znalezc parametry K oraz n tak, aby najlepiej oddawaly pomiary
#funckja celu ma postac funkcji dwoch zmiennych (K oraz n):
#f(n,K) = [sigma_pl - (K * epsilon**n )] -> min

from scipy import optimize
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
filename = 'tension_test_data'
data = pd.read_csv(
    'tension_test_data.csv',
    sep= ';'
)
print(data.columns)
eps = np.array(data.loc[:]['true_strain']) #NIE ZGADZALA SIE NAZWA KOLUMNY
sigma = np.array(data.loc[:]['plastic_stress_MPa']) 

eps_mod = np.linspace(0, max(eps), 100)

def ludwik(eps, b0, b1, b2):
    return b0 + b1 * eps**b2
#defining residual
def res(beta):
    return sigma - ludwik(eps, *beta)

beta0 = (1.0, 1.0, .5)

results = optimize.least_squares(res, beta0) #Tutaj dziwnie rozpakowane byly elementy rozwiazania optymalnego, nieobslugiwane chyba przez Python 3.13?

beta_opt = results.x

print("Beta =")
print(beta_opt)

corr_matrix = np.corrcoef(sigma, ludwik(eps, *beta_opt))
R = corr_matrix[0,1]
print('Coefficient of determination, R2 = ')
print(R**2)

fig = plt.figure()
plt.plot(eps, sigma, 'ro', label = 'Tension test data')
plt.plot(eps_mod, ludwik(eps_mod, *beta_opt),label= "Ludwik: $\sigma_{pl} = \\beta_0 + \\beta_1 \cdot \epsilon^{\\beta_2}$")
plt.xlabel("$\epsilon$")
plt.ylabel("$\sigma_{pl}$")
fig.legend(bbox_to_anchor=(.13,.87), loc = 2)
plt.show()