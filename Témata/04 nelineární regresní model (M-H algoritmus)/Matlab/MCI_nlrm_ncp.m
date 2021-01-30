function MCI=MCI_nlrm_ncp(results,S,num)
%Monte Carlo integrace pro NLRM s NCP
%   results ... strukturni promenna (vysledek funkce nlrm_ncp
%   S ... pocet replikaci
%   num ... cislo parametru pro ktery chceme MC integraci provest

%posteriorni stredni hodnota a rozptyl
theta_S = results.b1(num)+results.bstd1(num)*trnd(results.nu1,S,1);
E_g_S = mean(theta_S);
var_g_S = mean(theta_S.^2)-E_g_S.^2; %lze pouzit rovnez primo funkci var

MCI.theta_S = theta_S;
MCI.mean = E_g_S;
MCI.var = var_g_S;

%NSE;
NSE=sqrt(var_g_S/S); %numerical standard error

MCI.NSE=NSE;