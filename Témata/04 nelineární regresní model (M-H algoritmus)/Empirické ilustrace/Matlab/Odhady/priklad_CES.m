%Ilustrace Koop (2000) - kapitola 5 - prvni cast
%Independence Chain Metropolis-Hastings algoritmus
%Random Walk Chain Metropolis-Hastings algoritmus


clear all
close all
clc

%nacteni dat
load ch5data.out;
output = ch5data(:,1);
lab=ch5data(:,2);
cap=ch5data(:,3);
[n,k]=size(ch5data);

%% Vysvetlujici a vysvetlovana promenna pro CES funkci
X = [ones(n,1) lab cap]; 
y=output;

tic;

%*predpokladame neinformativni apriorni hustotu pro gamma a h ve tvaru 1/h
%*zamerime se na marginalni posteriorni hustotu pro vektor gamma
%pro sestrojeni kandidatske hustoty -- hledame modus posteriorni hustoty +
%Hessian nebo alternativne maximum verohodnostni funkce
%*lze vyuzit ruzne optimalizacni funkce (algoritmy)-- napr. Matlabovskou 'fminunc' nebo
%'pow_min' z ekonometrickeho toolboxu; 
%*pro tyto ucely je treba definovat funkci marginalni posteriorni hustoty
%nebo verohodnostni funkci (-1) nasobek teto funkce, nebot oba algoritmy hledaji minimum
%pro omezeni numerickych nepresnosti je funkce definovana jako log hustoty

gamma0 = [1;1;1;1];

[x,fval,exitflag,output,grad,hessian]=fminunc(@(x)CES_post(x,y,X),gamma0);

%result = pow_min('CES_post',gamma0,[],y,X); %muze nastat problem, ze hessian neni pozitivne definitni

postmode = x;
postvar = inv(hessian);


%pocet replikace
S = 25000;
S0 = 5000;
S1 = S-S0;

%kandidatska hustota - t-rozdeleni se stredni hodnotou a matici sigma danou
%vysledky posteriorni maximalizece
%experimentujeme s ruznymi hodnotami c a dof (stupne volnosti)
%*v pripade Independent Chain spise experimentujeme s dof
%*navic muzeme i prvotni odhady var(gamma|y) pouzit pro novou definici
%postvar, pripadne dale kombinovat s volbou c
c=1;
d=0.1;
dof=5;
vscale_ic= c*postvar;
vchol=chol(vscale_ic);
vchol=vchol';
vscale_rw= d*postvar;

gamma_ic = zeros(4,S); %vzorky pro ICMH algoritmus 
gamma_rw = zeros(4,S); %vzorky pro RWMH algoritmus 


count_ic = 0; %èítaè akceptovaných kandidátù
count_rw = 0;

gamma0_ic = [-10;-10;-10;-10];  
gamma0_rw = postmode;

%% prvni beh M-H
fprintf('Posteriorni simulace\n');
fprintf('0%%       50%%      100%% \n');
fprintf('|');
%Independent Chain
    gamma_can_ic = postmode+ vchol*tdis_rnd(4,dof); %kandidat
    log_accept_ic = min(-CES_post(gamma_can_ic,y,X)+CES_post(gamma0_ic,y,X)+log_mvtpdf(gamma0_ic,postmode,vscale_ic,dof) - log_mvtpdf(gamma_can_ic,postmode,vscale_ic,dof),0);
    if log_accept_ic > log(rand)
        gamma_ic(:,1)=gamma_can_ic;
        count_ic=count_ic+1;
    else
        gamma_ic(:,1)=gamma0_ic;
    end
   
%Random Walk M-H
    gamma_can_rw = gamma0_rw + norm_rnd(vscale_rw); %kandidat
    log_accept_rw = min(-CES_post(gamma_can_rw,y,X)+CES_post(gamma0_rw,y,X),0);
    if log_accept_rw > log(rand)
        gamma_rw(:,1)=gamma_can_rw;
        count_rw=count_rw+1;
    else
        gamma_rw(:,1)=gamma0_rw;
    end

    
    
for i=2:S
    gamma_can_ic = postmode + vchol*tdis_rnd(4,dof); %kandidat
    log_accept_ic = min(-CES_post(gamma_can_ic,y,X)+CES_post(gamma_ic(:,i-1),y,X)+log_mvtpdf(gamma_ic(:,i-1),postmode,vscale_ic,dof) - log_mvtpdf(gamma_can_ic,postmode,vscale_ic,dof),0);
    if log_accept_ic > log(rand)
        gamma_ic(:,i)=gamma_can_ic;
        count_ic=count_ic+1;
    else
        gamma_ic(:,i)=gamma_ic(:,i-1);
    end

    gamma_can_rw = gamma_rw(:,i-1)+ norm_rnd(vscale_rw); %kandidat
    log_accept_rw = min(-CES_post(gamma_can_rw,y,X)+CES_post(gamma_rw(:,1),y,X),0);
    if log_accept_rw > log(rand)
        gamma_rw(:,i)=gamma_can_rw;
        count_rw=count_rw+1;
    else
        gamma_rw(:,i)=gamma_rw(:,i-1);
    end
    
    if mod(i/S*100,5) == 0;    
      fprintf('|');
    end
end

fprintf('\n');
fprintf('Hotovo!\n');

%vyhozeni prvnich S0 replikaci
gamma_ic=gamma_ic(:,S0+1:end);
gamma_rw=gamma_rw(:,S0+1:end);

mean_ic = mean(gamma_ic,2); %posteriorni stredni hodnoty
mean_rw = mean(gamma_rw,2); 

var_ic = mean(gamma_ic.^2,2)-mean_ic.^2; %posteriorni rozptyly
var_rw = mean(gamma_rw.^2,2)-mean_rw.^2;

accept_ratio_ic = count_ic/S;
accept_ratio_rw = count_rw/S;



%% KONVERGENCNI DIAGNOSTIKY - Gewekova CD pomoci funkce momentg
alldraws = gamma_ic';
%The function momentg is taken from LeSage's toolbox
%it inputs all Gibbs draws and produces posterior
%mean, standard deviation, nse and rne
%it calculates what the book calls S(0) in various ways
%see momentg.m for more details
%calculate Geweke convergence diagnostic based on first .1
%and last .4 of draws
idraw1= round(.1*S1); %prvnich 10%
result = momentg(alldraws(1:idraw1,:));
meansa=[result.pmean]';
nsea=[result.nse1]';
idraw2= round(.6*S1)+1; %poslednich 40%
result = momentg(alldraws(idraw2:S1,:));
meansb=[result.pmean]';
nseb=[result.nse1]';

CD_ic = (meansa - meansb)./(nsea+nseb);
result = momentg(alldraws);
nse_ic=[result.nse1]';

alldraws = gamma_rw';
%The function momentg is taken from LeSage's toolbox
%it inputs all Gibbs draws and produces posterior
%mean, standard deviation, nse and rne
%it calculates what the book calls S(0) in various ways
%see momentg.m for more details
%calculate Geweke convergence diagnostic based on first .1
%and last .4 of draws
idraw1= round(.1*S1); %prvnich 10%
result = momentg(alldraws(1:idraw1,:));
meansa=[result.pmean]';
nsea=[result.nse1]';
idraw2= round(.6*S1)+1; %poslednich 40%
result = momentg(alldraws(idraw2:S1,:));
meansb=[result.pmean]';
nseb=[result.nse1]';

CD_rw = (meansa - meansb)./(nsea+nseb);
result = momentg(alldraws);
nse_rw=[result.nse1]';


t = toc;

fprintf('\n\n------------------------------------------------------\n')
fprintf('Posteriorni charakteristiky Independence Chain a Random Walk Chain M-H algoritmu\n')
fprintf('------------------------------------------------------\n')
fprintf('                 Independence chain M-H\n')
fprintf('                 Post.mean          Post.Var           NSE       Geweke CD\n')
fprintf('------------------------------------------------------\n')
for i = 1:4
fprintf('gamma_%1.0f        % 8.4f           % 8.4f          % 8.4f     % 8.4f \n',[i mean_ic(i) var_ic(i) nse_ic(i) CD_ic(i)]);
end
fprintf('   ------------------------------------\n');
fprintf('Accept.Ratio:  %8.4f\n', accept_ratio_ic)

fprintf('------------------------------------------------------\n')
fprintf('                 Random Walk Chain M-H\n')
fprintf('                 Post.mean          Post.Var           NSE       Geweke CD\n')
fprintf('------------------------------------------------------\n')
for i = 1:4
fprintf('gamma_%1.0f        % 8.4f           % 8.4f          % 8.4f     % 8.4f \n',[i mean_rw(i) var_rw(i) nse_rw(i) CD_rw(i)]);
end
fprintf('   ------------------------------------\n');
fprintf('Accept.Ratio:  %8.4f\n\n', accept_ratio_rw)
fprintf('   ------------------------------------\n');
fprintf('Celkovy pocet replikaci:       %12.4f\n',S)
fprintf('Vyhozenych replikaci:          %12.4f\n',S0)
fprintf('   ------------------------------------\n');
fprintf('Dosavadni cas (s): %12.2f \n',t);
fprintf('***************************************************************************\n\n\n');







%% Posteriorni p-hodnota
%pocitame f(X,gamma)
%gamma vezmeme napr. ze vzorku RWMH
fprintf('   ------------------------------------\n');
fprintf('   Vypocet posteriorni predikcni p-hodnoty pro sikmost a spicatost\n');
fprintf('   ------------------------------------\n');


ystar = zeros(n,S1);
fgamma_rw=zeros(n,S1);
eps = zeros(n,S1);
eps_star = zeros(n,S1);

skew = zeros(1,S1);
kurt = zeros(1,S1);
skew_star = zeros(1,S1);
kurt_star = zeros(1,S1);



for ii=1:S1
  fgamma_rw(:,ii) = gamma_rw(1,ii)*X(:,1)+(gamma_rw(2,ii)*X(:,2).^gamma_rw(4,ii)+gamma_rw(3,ii)*X(:,3).^gamma_rw(4,ii)).^(1/gamma_rw(4,ii));
  s12 = (y-fgamma_rw(:,ii))'*(y-fgamma_rw(:,ii))/n;
  ystar(:,ii) = fgamma_rw(:,ii) + sqrt(s12)*tdis_rnd(n,n); %simulace umelych dat
  eps(:,ii) = y-fgamma_rw(:,ii);
  eps_star(:,ii) = ystar(:,ii)-fgamma_rw(:,ii);
  
  skew(1,ii)= sqrt(n)*sum(eps(:,ii).^3)/(sum(eps(:,ii).^2)^1.5);
  kurt(1,ii)= n*sum(eps(:,ii).^4)/(sum(eps(:,ii).^2))^2 -3;
  skew_star(1,ii)= sqrt(n)*sum(eps_star(:,ii).^3)/(sum(eps_star(:,ii).^2)^1.5);
  kurt_star(1,ii)= n*sum(eps_star(:,ii).^4)/(sum(eps_star(:,ii).^2))^2 -3;
end

%Stredni hodnota sikmosti v pozorovanych datech
mskew=mean(skew);
%Stredni hodnota spicatosti v pozorovanych datech
mkurt=mean(kurt);

%posteriorni predikcni p-hodnota (oboustranna)
p_skew = length(find(abs(skew_star)<abs(mskew)))/S1;
p_kurt = length(find(abs(kurt_star)<abs(mkurt)))/S1;

t=toc;

fprintf('p-value sikmosti:   %8.4f\n', 1-p_skew)
fprintf('p-value spicatosti: %8.4f\n', 1-p_kurt)
fprintf('   ------------------------------------\n');
fprintf('Celkovy cas (s): %12.2f \n',t);
fprintf('***************************************************************************\n\n\n');



%% Vykresleni obrazku
fig = figure;
hist(skew_star,25)
title('Figure 5.1: Posterior Predictive Density for Skewness')
xlabel('Sikmost')
hold on
plot(mskew,0,'y*')
annotation(fig,'textarrow',[0.7822 0.5946],[0.2822 0.1347],...
    'TextEdgeColor','none',...
    'String',{'Pozorovana','sikmost'});


fig = figure;
hist(kurt_star,25)
title('Figure 5.2: Posterior Predictive Density for Kurtosis')
xlabel('Spicatost')
hold on
plot(mkurt,0,'y*')
annotation(fig,'textarrow',[0.5428 0.3723],[0.2777 0.1364],...
    'TextEdgeColor','none',...
    'String',{'Pozorovana','spicatost'});



