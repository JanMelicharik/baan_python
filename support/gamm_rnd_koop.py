# Funkce gamma, zapis dle Koopa (2003)
# Python funkce numpy.random.gamma() pouziva jinou notaci
# Vysvetleni prevodu Koopovy notace na notaci numpy:
# https://stats.stackexchange.com/questions/30655/a-question-about-parameters-of-gamma-distribution-in-bayesian-econometrics
# Analogicke reseni k funkci gamm_rnd_Koop_2.m

from numpy.random import gamma

def gamm_rnd_koop(mean, degrees_of_freedom):
    """
        mean ................. h_1 (posteriorni odhad presnostni)
        degrees_of_freedom ... nu_1
        size ................. cele cislo nebo tuple rozmeru
    """
    shape = degrees_of_freedom / 2
    # scale = degrees_of_freedom / (2 * mean)       # not sure about this setup
    scale = 2 * mean / degrees_of_freedom           # pro prvni cviceni vyhovuje tento
    return gamma(shape=shape, scale=scale)
