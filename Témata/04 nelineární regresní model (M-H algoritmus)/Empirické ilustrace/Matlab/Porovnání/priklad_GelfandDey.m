%Ilustrace Koop (2003) - kapitola 5 - druha cast
%Metoda Gelfanda a Deye

clear all
close all
clc

%% nacteni dat
load ch5data.out;
output = ch5data(:,1);
lab=ch5data(:,2);
cap=ch5data(:,3);
n = length(output);



%% Vysvetlujici a vysvetlovana promenna pro CES funkci
X = [ones(n,1) lab cap]; 
y=output;

tic;

%*predpokladame informativni nezavislou normalni-gama apriorni hustotu
%pro gamma a h
%*predchozi cast doplnime o generovani podminenych hustot z h
%*modifikuje se drobne jadrova posteriorni podminena hustota pro gamu (viz
%Koop (2003) - (5.24)
%stejny postup pro sestrojeni kandidatske hustoty -- hledame modus posteriorni hustoty +
%Hessian nebo alternativne maximum verohodnostni funkce
%*lze vyuzit ruzne optimalizacni funkce (algoritmy)-- napr. Matlabovskou 'fminunc' nebo
%'pow_min' z ekonometrickeho toolboxu; 
%*pro tyto ucely je treba definovat funkci marginalni posteriorni hustoty
%nebo verohodnostni funkci (-1) nasobek teto funkce, nebot oba algoritmy hledaji minimum
%pro omezeni numerickych nepresnosti je funkce definovana jako log hustoty



%% Apriorni hyperparametry pro M1 (M2 - nasledne zkraceny o ctvrty parametr)
gamma0 = [1;1;1;1];
V0 = 0.25*eye(4);
nu0 = 12;
s02 = 1/10;

h0 = 1/s02;

[x,fval,exitflag,output,grad,hessian]=fminunc(@(x)CES_post(x,y,X,h0,gamma0,V0),gamma0);

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
d=0.8;
vscale_rw= d*postvar;

gamma_rw = zeros(4,S); %vzorky pro RWMH algoritmus 
h_rw = zeros(1,S);

%pocitani akceptovaných kandidátù
count_rw = 0;


gamma0_rw = postmode;


%% prvni beh M-H
fprintf('Posteriorni simulace MH \n');
fprintf('0%%       50%%      100%% \n');
fprintf('|');
   
%generovani h
nu1 = n+nu0;
f = gamma0_rw(1)*X(:,1)+(gamma0_rw(2)*X(:,2).^gamma0_rw(4)+gamma0_rw(3)*X(:,3).^gamma0_rw(4)).^(1/gamma0_rw(4));
s12 = ((y-f)'*(y-f)+nu0*s02)/nu1;
h_rw(:,1)=gamm_rnd(1,1,.5*nu1,.5*nu1*s12); %Koop 5.23


%Random Walk M-H
    gamma_can_rw = gamma0_rw + norm_rnd(vscale_rw); %kandidat
    log_accept_rw = min(-CES_post(gamma_can_rw,y,X,h_rw(:,1),gamma0,V0)+CES_post(gamma0_rw,y,X,h_rw(:,1),gamma0,V0),0);
    if log_accept_rw > log(rand)
        gamma_rw(:,1)=gamma_can_rw;
        count_rw=count_rw+1;
    else
        gamma_rw(:,1)=gamma0_rw;
    end

   
for i=2:S

%generovani h   
f = gamma_rw(1,i-1)*X(:,1)+(gamma_rw(2,i-1)*X(:,2).^gamma_rw(4,i-1)+gamma_rw(3,i-1)*X(:,3).^gamma_rw(4,i-1)).^(1/gamma_rw(4,i-1));
s12 = ((y-f)'*(y-f)+nu0*s02)/nu1;
h_rw(:,i)=gamm_rnd(1,1,.5*nu1,.5*nu1*s12); %Koop 5.23
        
    gamma_can_rw = gamma_rw(:,i-1)+ norm_rnd(vscale_rw); %kandidat
    log_accept_rw = min(-CES_post(gamma_can_rw,y,X,h_rw(:,i),gamma0,V0)+CES_post(gamma_rw(:,i-1),y,X,h_rw(:,i),gamma0,V0),0);
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
gamma_rw=gamma_rw(:,S0+1:end);
h_rw=h_rw(:,S0+1:end);


%posteriorni stredni hodnoty
mean_rw = mean(gamma_rw,2); 
mean_h = mean(h_rw); 

%posteriorni rozptyly
var_rw = mean(gamma_rw.^2,2)-mean_rw.^2;
var_h = mean(h_rw.^2)-mean_h.^2;

accept_ratio_rw = count_rw/S;



%% KONVERGENCNI DIAGNOSTIKY - Gewekova CD pomoci funkce momentg
alldraws = [gamma_rw' h_rw'];
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

%% Vysledky

fprintf('\n\n------------------------------------------------------\n')
fprintf('Posteriorni charakteristiky Random Walk Chain M-H algoritmu\n')
fprintf('------------------------------------------------------\n')
fprintf('                 Random Walk Chain M-H\n')
fprintf('                 Post.mean          Post.Var           NSE       Geweke CD\n')
fprintf('------------------------------------------------------\n')
for i = 1:4
fprintf('gamma_%1.0f        % 8.4f           % 8.4f          % 8.4f     % 8.4f \n',[i mean_rw(i) var_rw(i) nse_rw(i) CD_rw(i)]);
end

fprintf('------------------------------------------------------\n')
fprintf('                 Parametr h\n')
fprintf('                 Post.mean          Post.Var           NSE       Geweke CD\n')
fprintf('------------------------------------------------------\n')
fprintf('h                % 8.4f           % 8.4f          % 8.4f     % 8.4f \n',[mean_h var_h nse_rw(5) CD_rw(5)]);
fprintf('   ------------------------------------\n');
fprintf('Accept.Ratio:  %8.4f\n', accept_ratio_rw)
fprintf('   ------------------------------------\n');
fprintf('Celkovy pocet replikaci:       %12.4f\n',S)
fprintf('Vyhozenych replikaci:          %12.4f\n',S0)
fprintf('   ------------------------------------\n');
fprintf('Dosavadni cas (s): %12.2f \n',t);
fprintf('***************************************************************************\n\n\n');





%% Metoda Gelfanda-Deye pro model M2
k = 5;
theta = [gamma_rw; h_rw];
mean_theta = mean(theta,2);
cov_theta = cov(theta');

pp = 0.01;

chi = chis_inv(1-pp,k); %(1-pp)% kvantil chi-kvadrat

g_theta = 0;

for ii=1:S1
    Theta = (mean_theta-theta(:,ii))'*inv(cov_theta)*(mean_theta-theta(:,ii));
    ff = 1/(pp*(2*pi)^(k/2))*det(cov_theta)^(-1/2)*exp(-1/2*Theta)*(Theta<=chi);
    prior = CES_prior(gamma_rw(:,ii),h_rw(:,ii),gamma0,V0,h0,nu0);
    like = CES_like(gamma_rw(:,ii),h_rw(:,ii),y,X);
    g_theta = g_theta+ff/(prior*like);
end

%inverze marginalni verohodnosti
inv_marg = g_theta/S1;
marg_M2 = 1/inv_marg;








%% Odhad modelu pro M1: gamma4 = 1
k=3;
gamma0 = [1;1;1];
V0 = 0.25*eye(3);
nu0 = 12;
s02 = 1/10;
h0 = 1/s02;

S = 25000;
S0 = 5000;
S1 = S-S0;


gamma_Gibbs = zeros(k,S);
h_Gibbs = zeros(1,S);


%% Gibbsuv vzorkovac
%prvni beh - je mozno dat vse do jednoho cyklu, kdy budeme mit S+1
%replikaci, kdy prvni vektor parametru odpovida nasim pocatecnim podminkam
h_Gibbs0 = h0;

V1 = inv(inv(V0)+h_Gibbs0*X'*X);
b1 = V1*(inv(V0)*gamma0+h_Gibbs0*X'*y);

gamma_Gibbs(:,1)=b1+norm_rnd(V1);

nu1 =n+nu0;
s12 = 1/nu1*((y-X*gamma_Gibbs(:,1))'*(y-X*gamma_Gibbs(:,1))+nu0*s02);
h_Gibbs(:,1)=gamm_rnd(1,1,.5*nu1,.5*nu1*s12);

fprintf('Posteriorni simulace - Gibbs \n');
fprintf('0%%       50%%      100%% \n');
fprintf('|');
for i=2:S
   V1 = inv(inv(V0)+h_Gibbs(1,i-1)*X'*X);
   b1 = V1*(inv(V0)*gamma0+h_Gibbs(1,i-1)*X'*y);
   gamma_Gibbs(:,i)=b1+norm_rnd(V1); %spolehlivejsi generator!!!
    %nu1 =n+nu0; %v ramci cyklu se jiz nemeni
   s12 = 1/nu1*((y-X*gamma_Gibbs(:,i))'*(y-X*gamma_Gibbs(:,i))+nu0*s02);
   h_Gibbs(:,i)=gamm_rnd(1,1,.5*nu1,.5*nu1*s12);
        
        %graficky detail pro zobrazeni prubehu simulace
    if mod(i/S*100,5) == 0;    
        fprintf('|');
    end
end
fprintf('\n');
fprintf('Hotovo!\n');

%% vyhozeni prvnich S0 vzorku
gamma_Gibbs = gamma_Gibbs(:,S0+1:end);
h_Gibbs = h_Gibbs(:,S0+1:end);

mean_Gibbs = mean(gamma_Gibbs,2); %posteriorni stredni hodnota 
mean_h_Gibbs = mean(h_Gibbs);
var_Gibbs = mean(gamma_Gibbs.^2,2)-mean_Gibbs.^2; %posteriorni rozptyl
var_h_Gibbs = mean(h_Gibbs.^2)-mean_h_Gibbs.^2;


%% KONVERGENCNI DIAGNOSTIKY - Gewekova CD pomoci funkce momentg
alldraws = [gamma_Gibbs' h_Gibbs'];
%The function momentg is taken from LeSage's toolbox
%it inputs all Gibbs draws and produces posterior
%mean, standard deviation, nse and rne
%it calculates what the book calls S(0) in various ways
%see momentg.m for more details
result = momentg(alldraws);
means=[result.pmean]';
stdevs=[result.pstd]';
nse=[result.nse]';
nse1=[result.nse1]';
nse2=[result.nse2]';
nse3=[result.nse3]';
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

CD_Gibbs = (meansa - meansb)./(nsea+nseb);
result = momentg(alldraws);
nse_Gibbs=[result.nse1]';


t = toc;

%% Prezentace vysledku
fprintf('\n\n------------------------------------------------------\n')
fprintf('Posteriorni charakteristiky pro linearni model\n')
fprintf('------------------------------------------------------\n')
fprintf('                 Gibbsuv vzorkovac\n')
fprintf('                 Post.mean          Post.Var           NSE       Geweke CD\n')
fprintf('------------------------------------------------------\n')
for i = 1:3
fprintf('gamma_%1.0f        % 8.4f           % 8.4f          % 8.4f     % 8.4f \n',[i mean_Gibbs(i) var_Gibbs(i) nse_Gibbs(i) CD_Gibbs(i)]);
end

fprintf('------------------------------------------------------\n')
fprintf('                 Parametr h\n')
fprintf('                 Post.mean          Post.Var           NSE       Geweke CD\n')
fprintf('------------------------------------------------------\n')
fprintf('h                % 8.4f           % 8.4f          % 8.4f     % 8.4f \n',[mean_h_Gibbs var_h_Gibbs nse_Gibbs(4) CD_Gibbs(4)]);
fprintf('   ------------------------------------\n');
fprintf('Celkovy pocet replikaci:       %12.4f\n',S)
fprintf('Vyhozenych replikaci:          %12.4f\n',S0)
fprintf('   ------------------------------------\n');
fprintf('Dosavadni cas (s): %12.2f \n',t);
fprintf('***************************************************************************\n\n\n');


%% Metoda Gelfanda-Deye pro model M1

k = 4;
theta = [gamma_Gibbs; h_Gibbs];
mean_theta = mean(theta,2);
cov_theta = cov(theta');

pp = 0.01;

chi = chis_inv(1-pp,k); %(1-pp)% kvantil chi-kvadrat

g_theta = 0;

for ii=1:S1
    Theta = (mean_theta-theta(:,ii))'*inv(cov_theta)*(mean_theta-theta(:,ii));
    ff = 1/(pp*(2*pi)^(k/2))*det(cov_theta)^(-1/2)*exp(-1/2*Theta)*(Theta<=chi);
    prior = CES_prior(gamma_Gibbs(:,ii),h_Gibbs(:,ii),gamma0,V0,h0,nu0);
    like = CES_like_lin(gamma_Gibbs(:,ii),h_Gibbs(:,ii),y,X);
    g_theta = g_theta+ff/(prior*like);
end

%inverze marginalni verohodnosti
inv_marg = g_theta/S1;
marg_M1 = 1/inv_marg;




%% Porovnani modelu
fprintf('\n\n------------------------------------------------------\n')
fprintf('Porovnani modelu\n')
fprintf('------------------------------------------------------\n')
fprintf('Model              Marginalni verohodnost\n')
fprintf('------------------------------------------------------\n')
fprintf('M1:gamma_4=1       % 8.4f\n',marg_M1);
fprintf('M2:gamma_4<>1      % 8.4f\n',marg_M2);
fprintf('------------------------------------------------------\n')
fprintf('Bayesuv faktor     %12.4f\n',marg_M1/marg_M2)
fprintf('***************************************************************************\n\n\n');
