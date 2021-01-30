function y = gamrnd_Koop(mu,nu,M,N)
%Generovani nahodnych cisel z G(mu,nu)
%dle Koop(2003), vyuzivajici gammrnd
%statistickeho toolboxu
%M ... pocet radku
%N ... pocet sloupcu

A = nu/2;
B = 2*mu/nu;
y = gamrnd(A,B,M,N);
end
