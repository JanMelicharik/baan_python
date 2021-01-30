function y = gamm_rnd_Koop_2(mu,nu,m)
%Generator nahodnych cisel z G(mu,nu) rozdeleni
%dle Koop (2003) - vyuziva starsi gamm_rnd z Econ. toolboxu
%mu ... stredni hodnota
%nu ... pocet stupnu volnosti
%m x 1 ... rozmer nahodneho vektoru y

A = nu/2;
B = nu/(2*mu);
y = gamm_rnd(m,1,A,B);
end
