# Matematicke metody:
# Import jednotlivych matematickych metod abych nebylo nutne pouzivat
# plnou cestu funkce v dane knihovne.
print("Import matematickych funkci ... ", end="")
from math import exp
from math import log
from math import pi
from math import sqrt
from numpy import array
from numpy import diag
from numpy import eye
from numpy import mean
from numpy import ones
from numpy import std
from numpy import zeros
from numpy.linalg import det
from numpy.linalg import inv
from numpy.random import gamma
from numpy.random import normal
from numpy.random import randn
from numpy.random import uniform
from scipy.stats import multivariate_normal
from scipy.stats import norm
from scipy.stats import t as student
from scipy.stats.distributions import chi2
from pandas import read_excel
from pandas import read_csv
from pandas import DataFrame
print("OK")

# Podpurne a vizualizacni metody:
# Predevsim pro prehlednejsi vypis vysledku a praci se skripty.
print("Import podpurnych funkci ... ", end="")
from pdb import set_trace
from tabulate import tabulate
from termcolor import colored
from time import time
from warnings import filterwarnings
print("OK")

# Knihovny:
# Import celych knihoven kvuli zachovani dalsich rozsirujicich funkci,
# ktere tyto knihovny poskytuji.
print("Import knihoven ... ", end="")
import matplotlib.pyplot as plt
import numpy as np
import pandas
print("OK")

# Vlastni metody:
# Pomocne metody prevzate od doc. Daniela Nemce nebo z LeSageho
# ekonometrickeho toolboxu.
print("Import funkci ze slozky 'support' ... ", end="")
from support.a_post import a_post
from support.gamm_rnd_koop import gamm_rnd_koop
from support.gamm_rnd_koop2 import gamm_rnd_koop2
from support.geweke import geweke
from support.lik_ces import lik_ces
from support.log_post_ces import log_post_ces
from support.momentg import momentg
from support.my_boot import my_boot
from support.my_nlrm import my_nlrm
from support.new_geweke import new_geweke
from support.norm_rnd import norm_rnd
from support.nu_lambda_post import nu_lambda_post
from support.prior_ces import prior_ces
from support.progress_info import progress_bar
from support.wish_rnd import wish_rnd
print("OK")
