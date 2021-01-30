from pandas import read_excel, DataFrame
from math import sqrt, lgamma, log, pi, exp
from numpy import diag, ones, dot, array, asmatrix, shape, zeros, hstack, delete
from numpy.linalg import inv, det
from numpy import transpose as t
from numpy import sum as Sum
from scipy.stats import t as student

from Support.gamrnd_Koop import gamrnd_Koop
from Support.my_NLRM import my_NLRM
from Support.norm_rnd import norm_rnd

# mozna bude nutny doplnujici balicek xlrd (>>: pip install xldr)

# nacteni dat do DataFrame
df = read_excel('capm2_data.xlsx')

# definice zavisle promenne
y = df['GM'] - df['RKFREE']

# definice matice X, kde prvni sloupce obsahuje pouze 1
X = DataFrame(
    {
        'const' : ones(len(df)),
        'x_1' : df['MKT'] - df['RKFREE']
    }
)

# nastaveni apriornich hyperparametru
beta_0 = [0, 1]
cov_beta_0 = diag([0.05**2, 0.5**2])
s2_0 = 0.2**2
h_0 = 1/s2_0
nu_0 = 10
V_0 = cov_beta_0 * (nu_0 - 2)/nu_0 * h_0

# odhad NLRM s prirozene konjugovanym priorem
res_GM = my_NLRM(y, X, beta_0, V_0, h_0, nu_0)

# +------------------------+
# | 1) prezentace vysledku |
# +------------------------+
# PREZENTACE VYSLEDKU ZDE

# +-------------------------------+
# | 2) test hypotezy, ze beta = 1 |
# +-------------------------------+
# odhad omezeneho modelu za predpokladu, ze beta = 1
# y = alpha + 1 * x + epsilon ---> y - x = alpha + epsilon
y_a = (df['GM'] - df['RKFREE']) - (df['MKT'] - df['RKFREE'])
X_a = array(ones((len(df),1)))
V_0_a = asmatrix(V_0[1][1])
beta_0_a = asmatrix(beta_0[1])

# nutna uprava apriornich hustot - vstupuji zde hyperparametry jen pro prvni
# parametr (neovlivni to odhady, ale ovlivnilo by to vypocet marginalni verohodnosti)
res_GM_rest = my_NLRM(y_a,X_a,beta_0_a,V_0_a,h_0,nu_0)

# PREZENTACE VYSLEDKU ZDE

# logaritmus Bayesova faktoru porovnavajici model omezeny a neomezeny
log_BF = res_GM_rest['log_ML'] - res_GM['log_ML']
BF = exp(log_BF)
 
# PREZENTACE VYSLEDKU ZDE

# +-------------------------------+
# | 3) test hypotezy, ze beta > 1 |
# +-------------------------------+
# lze analyticky (z posteriorni marginalni hustoty pro beta ~ t-rozdeleni nebo
# simulacne h|y ~ G(h_1, nu_1) a beta|h,y ~ N(beta_1, 1/h * V_1))

# pocet simulaci
S = 1_000

# zasobnik pro simulovane hodnoty beta (po sloupcich)
# ulozen v transponovane forme narozdil od Matlabu, protoze lze jednoduseji
# zapisovat hodnoty v cyklu

beta_sim = array([[0],[0]])

for i in range(S):
    h_sim = gamrnd_Koop(res_GM['h_1'],res_GM['nu_1'],1,1)

    # kde je pouzita transpozice na vektor beta_1, protoze python nedokaze scitat
    # jinak orientovane vektory
    rand_pick = norm_rnd(inv(h_sim) * res_GM['V_1']) + t(res_GM['beta_1'])
    beta_sim = hstack((beta_sim, rand_pick))

# zde se maze prvni sloupec matice beta_sim, ktery byl pouzit pouze k nastaveni
# sirky matice pro jednodussi zapis dat
beta_sim = delete(beta_sim,0,1)

# Vypocet pravdepodobnosti, ze beta > 1
pr_beta = Sum(beta_sim[1] > 1)/S

# PREZENTACE VYSLEDKU ZDE

# analyticky vypocet pravdepodobnosti
# a) standardizace skalovaneho t-rozdeleni (p(beta|y)) pro druhy prvek
#    vektoru parametru beta 
zscore = (1 - res_GM['beta_1'][0,1])/(res_GM['b1_std'][1])

# b) vypocet odpovidajiciho kvantilu ze standardizovaneho centrovaneho
pr_beta_analyticky = 1 - student.pdf(zscore, res_GM['nu_1'])

# PREZENTACE VYSLEDKU ZDE
