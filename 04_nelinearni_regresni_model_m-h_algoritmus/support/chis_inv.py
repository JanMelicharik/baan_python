# Vraci inverzni kvantil z chi-kvadrat rozdeleni pro danou hodnotu x
# Pro uplny popis viz. tuto funkci pro MATLAB

def chis_inv(p, a):
    x = gamm_inv(p, a * 0.5) * 2