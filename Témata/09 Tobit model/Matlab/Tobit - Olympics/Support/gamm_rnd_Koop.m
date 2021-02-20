function y = gamm_rnd_Koop(mu,nu,m)
%Generator nahodnych cisel z G(mu,nu) rozdeleni
%dle Koop (2003) - vyuziva gamm_rnd z Econ. toolboxu
%mu ... stredni hodnota
%nu ... pocet stupnu volnosti
%m x 1 ... rozmer nahodneho vektoru y

A = nu/2;
B = 2*mu/nu;
y = B*gamm_rnd3(m,A);
end

