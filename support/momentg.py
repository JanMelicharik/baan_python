# Prepis funkce momentg z LeSageho ekonometrickeho toolboxu, pro detailni informace viz popis funkce na
# https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/54382/versions/1/previews/BSPSM_beta/auxiliary_files/momentg.m/index.html

# Pozn: Originalni script neni psan uplne stastnym zpusobem.

import sys
from termcolor import colored
from math import floor, sqrt
import numpy as np

def momentg(draws: np.ndarray):
    [ndraw, nvar] = draws.shape
    ng = 100
    if ndraw < ng:
        print(colored("momentg: Je potreba vice vzorku.", "red"))
        sys.exit(1)

    ntaper = [4, 8, 15]
    ns = floor(ndraw / ng)
    nuse = ns * ng

    pmean = [0] * nvar
    pstd = [0] * nvar
    nse = [0] * nvar
    rne = [0] * nvar
    nse1 = [0] * nvar
    rne1 = [0] * nvar
    nse2 = [0] * nvar
    rne2 = [0] * nvar
    nse3 = [0] * nvar
    rne3 = [0] * nvar

    for jf in range(nvar):
        cnt = 0
        cn = np.zeros(ng)
        cd = np.zeros(ng)
        cdn = np.zeros(ng)
        cdd = np.zeros(ng)
        cnn = np.zeros(ng)
        cvar = np.zeros(ng)
        td, tn, tdd, tnn, tdn, tvar = [0] * 6
        for ig in range(ns):
            gd, gn, gdd, gdn, gnn, gvar = [0] * 6
            for i_s in range(ns):
                g = draws[cnt][jf]
                ad = 1
                an = ad * g
                gd += ad
                gn += an
                gdn += ad * an
                gdd += ad ** 2
                gnn += an ** 2
                gvar += an * g
                cnt += 1
            
            td += gd
            tn += gn
            tdn += gdn
            tdd += gdd
            tnn += gnn
            tvar += gvar
            cn[ig] = gn / ns
            cn[ig] = gd / ns
            cn[ig] = gdn / ns
            cn[ig] = gdd / ns
            cn[ig] = gnn / ns
            cn[ig] = gvar / ns

        eg = tn / td
        varg = tvar / td - eg ** 2
        sdg = sqrt(varg) if varg > 0 else -1
        varnum = (tnn - 2 * eg * tdn + tdd * eg ** 2) / (td ** 2)
        sdnum = sqrt(varnum) if varnum > 0 else -1

        pmean[jf] = eg
        pstd[jf] = sdg
        nse[jf] = sdnum
        rne[jf] = varg / (nuse * varnum)

        barn = tn / nuse
        bard = td / nuse

        for ig in range(ng):
            cn[ig] -= barn
            cd[ig] -= bard

        rnn = [0] * (ng)
        rdd = [0] * (ng)
        rnd = [0] * (ng)
        rdn = [0] * (ng)
        for lag in range(ng):
            ann, add, an_d, adn = [0] * 4
            for ig in range(lag, ng):
                ann += cn[ig] * cn[ig - lag]
                add += cd[ig] * cd[ig - lag]
                an_d += cn[ig] * cd[ig - lag]
                adn += cd[ig] * cd[ig - lag]

            rnn[lag] = ann / ng
            rdd[lag] = add / ng
            rnd[lag] = an_d / ng
            rdn[lag] = adn / ng

        for mm in range(3):
            m = ntaper[mm]
            snn = rnn[0]
            sdd = rdd[0]
            snd = rnd[0]
            for lag in range(m - 1):
                att = 1 - lag / m
                snn += 2 * att * rnn[lag + 1]               # Zde se LeSage chtel asi vyjadrit dvema ruznymi zpusoby:
                sdd += 2 * att * rdd[lag + 1]               # 2 * a * b
                snd += att * (rnd[lag + 1] + rnd[lag + 1])  # a * (b + b) ... :) Zamysleni nad motivy necham na ctenari.

            varnum = ns * nuse * (snn - 2 * eg * snd + sdd * eg ** 2) / (td ** 2)
            sdnum = sqrt(varnum) if varnum > 0 else -1

            if mm == 0:
                nse1[jf] = sdnum
                rne1[jf] = varg / (nuse * varnum)
            if mm == 1:
                nse2[jf] = sdnum
                rne2[jf] = varg / (nuse * varnum)
            if mm == 2:
                nse3[jf] = sdnum
                rne3[jf] = varg / (nuse * varnum)

    return {
        "ndraw": ndraw,
        "nvar": nvar,
        "meth": "momentg",
        "pmean": np.array(pmean),
        "pstd ": np.array(pstd),
        "nse": np.array(nse),
        "rne": np.array(rne),
        "nse1": np.array(nse1),
        "rne1": np.array(rne1),
        "nse2": np.array(nse2),
        "rne2": np.array(rne2),
        "nse3": np.array(nse3),
        "rne3": np.array(rne3),
    }
