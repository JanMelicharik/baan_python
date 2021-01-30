function results = gamrnd_Koop(mu,nu,M,N)
%Vraci matici M x N nezavislych vyberu z Gamma rozdeleni dle znaceni Koop(2003)
%mu ... stredni hodnota, nu ... pocet stupnu volnosti
A = nu/2;
B = 2*mu/nu;
results = gamrnd(A,B,M,N);
end

