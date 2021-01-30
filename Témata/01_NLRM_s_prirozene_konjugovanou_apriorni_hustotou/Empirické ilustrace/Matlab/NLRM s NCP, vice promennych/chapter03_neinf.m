clear all
close all
clc

%Empiricka ilustrace - kapitola 3 - Koop 2003 - temer neinformativni
%apriorni hustota


%nastaveni generatoru nahodnych cisel dle casu spusteni
randn('seed',sum(100*clock));
rand('seed',sum(100*clock));


%% Nacteni dat
load hprice.txt

price = hprice(:,1); %prodejni cena domu
lot_size = hprice(:,2); %rozloha ve stopach ctverecnich
n_bed = hprice(:,3); %po�et lo�nic
n_bath = hprice(:,4); %po�et koupelen
n_storey = hprice(:,5); %po�et pater

y = price;
X = [ones(size(y)) lot_size n_bed n_bath n_storey];
Xstar = [1 5000 2 2 1];  %predikce

[n,k] = size(X);

%% apriorni hyperparametry 
b0 = [0 10 5000 10000 10000]';

s02 = 5000^2;
h0 = 1/s02;
varb0 = [1000000^2 1000000^2 1000000^2 1000000^2 1000000^2]';
nu0 = 1;
V0 = diag(varb0)*(nu0-2)/(nu0*s02);

results_neinf = nlrm_ncp(y,X,b0,V0,s02,nu0,Xstar);
prt_nlrm_ncp(results_neinf)
fprintf('   ------------------------------------\n\n');


%% Posteriorni podil sanci -- M_1:beta_j=0 M_2:beta_j\neq 0
%Problem pro informativni apriorni hustotu
BF=zeros(k,1);
for i=1:k
X1 = X;
X1(:,i)=[];
b1 = b0;
b1(i)=[];
varb1 = varb0;
varb1(i) = [];
V1 = diag(varb1)*(nu0-2)/(nu0*s02);
results1=nlrm_ncp(y,X1,b1,V1,s02,nu0); %nested model beta_j=0
    BF(i)=exp(results1.lmarglik-results_neinf.lmarglik);
end

HPDI95=HPDI_nlrm_ncp(results_neinf,0.95);
HPDI99=HPDI_nlrm_ncp(results_neinf,0.99);

probpos = zeros(k,1);
for i=1:k
  % pravdepodobnost beta kladne
  tnorm = -results_neinf.b1(i,1)/sqrt(results_neinf.s12*results_neinf.V1(i,i)); %normovani pro normalizovane studentovo rozdeleni
  probpos(i,1) = 1 - tcdf(tnorm,results_neinf.nu1); 
end


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
fprintf('Analyticky         %12.4f  %12.4f         --- \n',[results_neinf.b1(2) results_neinf.bstd1(2)]);
fprintf('---------------        \n');
fprintf('Pocet replikaci \n');
fprintf('---------------        \n');
for i=1:length(S)
    MCI = MCI_nlrm_ncp(results_neinf,S(i),2);
fprintf('S = %8u       %12.4f  %12.4f  %12.4f \n',[S(i) MCI.mean sqrt(MCI.var) MCI.NSE]);
end
fprintf('   ------------------------------------\n\n');

%% Predikce
for ii=1:size(Xstar,1)
    XXstar = Xstar(ii,:);
    yplot=0:1:140000;
    y_plotpred = my_tpdf(yplot,XXstar*results_neinf.b1,results_neinf.s12*(ones(size(XXstar,1))+XXstar*results_neinf.V1*XXstar'),results_neinf.nu1);
    figure
    plot(yplot,y_plotpred./sum(y_plotpred))
    title(['Predikcni hustota'])
    xlabel(['y^* pro [', num2str(XXstar), ']'])
    ylabel('Hustota pravdepodobnosti')
end