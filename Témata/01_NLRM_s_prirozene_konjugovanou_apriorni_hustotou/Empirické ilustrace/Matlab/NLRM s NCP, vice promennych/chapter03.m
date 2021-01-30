clear all
close all
clc

%Empiricka ilustrace - kapitola 3 - Koop 2003

%nastaveni generatoru nahodnych cisel dle casu spusteni
randn('seed',sum(100*clock));
rand('seed',sum(100*clock));

%% Nacteni dat
load hprice.txt

price = hprice(:,1); %prodejni cena domu
lot_size = hprice(:,2); %rozloha ve stopach ctverecnich
n_bed = hprice(:,3); %poèet ložnic
n_bath = hprice(:,4); %poèet koupelen
n_storey = hprice(:,5); %poèet pater

y = price; %zavisla promenna
X = [ones(size(y)) lot_size n_bed n_bath n_storey]; %matice planu
Xstar = [1 5000 2 2 1];  %hodnoty regresoru pro vypocet predikce

[n,k] = size(X);

%% apriorni hyperparametry 
b0 = [0 10 5000 10000 10000]';

s02 = 5000^2;
h0 = 1/s02;
varb0 = [10000^2 25 2500^2 5000^2 5000^2]';
nu0 = 5;
V0 = diag(varb0)*(nu0-2)/(nu0*s02);

results_inf = nlrm_ncp(y,X,b0,V0,s02,nu0,Xstar); %vyocet normalniho LRM s NCP
prt_nlrm_ncp(results_inf) %funkce pro vygenerovani prehledne tabulky s posteriornimi vysledky
fprintf('   ------------------------------------\n\n');


%% Posteriorni podil sanci -- M_1:beta_j=0 M_2:beta_j\neq 0
BF=zeros(k,1); %inicializace promenne pro ukladani Bayesovho faktoru
for i=1:k
X1 = X;
X1(:,i)=[];
b1 = b0;
b1(i)=[];
varb1 = varb0;
varb1(i) = [];
V1 = diag(varb1)*(nu0-2)/(nu0*s02);
results1=nlrm_ncp(y,X1,b1,V1,s02,nu0); %nested model beta_j=0
    BF(i)=exp(results1.lmarglik-results_inf.lmarglik);
end

HPDI95=HPDI_nlrm_ncp(results_inf,0.95); %vypocet 95% HPDI
HPDI99=HPDI_nlrm_ncp(results_inf,0.99); %vypocet 99% HPDI

probpos = zeros(k,1); %inicializace promenne pro ukladani posteriorni pravdepodobnosti beta_j>0
for i=1:k
  % pravdepodobnost beta kladne
  tnorm = -results_inf.b1(i,1)/sqrt(results_inf.s12*results_inf.V1(i,i)); %normovani pro normalizovane studentovo rozdeleni
  probpos(i,1) = 1 - tcdf(tnorm,results_inf.nu1); 
end

%Vypis vysledku
fprintf('Porovnani modelu zahrnujici Beta\n');
fprintf('        p(beta_j>0)           95%% HPDI                   99%% HPDI       Posterior Odds pro Beta_j=0\n');
for i=1:k
fprintf('Beta %2u   %5.4f   [%11.3f %11.3f]  [%11.3f  %11.3f]      %8.4f         \n',[i probpos(i) HPDI95(i,:) HPDI99(i,:) BF(i)]);
end
fprintf('   ------------------------------------\n\n');

%% Monte Carlo integrace pro Beta_2
S=[10 100 1000 10000 100000];
fprintf('Posteriorni vysledky pro Beta_2 pocitane ruznymi zpusoby\n');
fprintf('                         Mean      Std. Deviation     NSE \n');
fprintf('Analyticky         %12.4f  %12.4f         --- \n',[results_inf.b1(2) results_inf.bstd1(2)]);
fprintf('---------------        \n');
fprintf('Pocet replikaci \n');
fprintf('---------------        \n');
for i=1:length(S)
    MCI = MCI_nlrm_ncp(results_inf,S(i),2);
fprintf('S = %8u       %12.4f  %12.4f  %12.4f \n',[S(i) MCI.mean sqrt(MCI.var) MCI.NSE]);
end
fprintf('   ------------------------------------\n\n');

%% Predikce
for ii=1:size(Xstar,1)
    XXstar = Xstar(ii,:);
    yplot=0:1:140000;
    y_plotpred = my_tpdf(yplot,XXstar*results_inf.b1,results_inf.s12*(ones(size(XXstar,1))+XXstar*results_inf.V1*XXstar'),results_inf.nu1);
    figure
    plot(yplot,y_plotpred./sum(y_plotpred))
    title(['Predikcni hustota'])
    xlabel(['y^* pro [',num2str(XXstar),']'])
    ylabel('Hustota pravdepodobnosti')
end