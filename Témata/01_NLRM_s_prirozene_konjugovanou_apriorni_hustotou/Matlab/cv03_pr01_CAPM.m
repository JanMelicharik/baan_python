%% Bayesovsky odhad NLRM s prirozene konjugovanou apriorni hustotou (naturally conjugate prior = NCP)
clear           %vymaze vsechny promenne z pameti
close all       %zavre vsechna okna (s obrazky)
clc             %vymaze command window (prikazove okno)

%% Nastaveni cest k podpurnym funkcim a datum
addpath('..\Support') % = nastavi cestu relativne k aktualnimu adresari

%% Nacteni dat a priprava promennych
load capm2.mat

y = data.GM-data.RKFREE; % zavisla promenna
X = [ones(size(y)) data.MKT-data.RKFREE]; % matice planu (konstanta a 1 regresor)

%% Nastaveni apriornich hyperparametru
% beta|h ~ N(beta_0,V_0)
% h ~ G(h_0,nu_0)
beta_0 = [0
          1];      
cov_beta_0 = diag([0.05^2,0.5^2]);
s2_0 = 0.2^2;
h_0 = 1/s2_0;
nu_0 = 10;
V_0 = cov_beta_0*(nu_0-2)/nu_0*h_0;

%% Odhad NLRM s NCP pomoci funkce my_NLRM
res_GM = my_NLRM(y,X,beta_0,V_0,h_0,nu_0);

%% 1. prezentace vysledku
fprintf('Parametr    Prior   Prior std.    Posterior   Posterior std.\n');
fprintf('============================================================\n');
fprintf('Alpha       %6.4f   %6.4f     %6.4f  %6.4f\n',res_GM.beta_0(1),...
    res_GM.b0_std(1),res_GM.beta_1(1),res_GM.b1_std(1));
fprintf('Beta        %6.4f   %6.4f     %6.4f  %6.4f\n',res_GM.beta_0(2),...
    res_GM.b0_std(2),res_GM.beta_1(2),res_GM.b1_std(2));
fprintf('h           %6.4f   %6.4f     %6.4f  %6.4f\n',res_GM.h_0,...
    res_GM.h0_std,res_GM.h_1,res_GM.h1_std);

%% 2. test hypotezy beta = 1
%odhad omezeneho modelu za predpokladu beta = 1:
% y = alpha+ 1*x + epsilon -> y-x = alpha +epsilon
y = (data.GM-data.RKFREE)-(data.MKT-data.RKFREE);
X = ones(size(y)); % matice planu jiz jenom konstanta
res_GM_rest = my_NLRM(y,X,beta_0(1),V_0(1,1),h_0,nu_0); 
        %nutna uprava apriornich hustot ... vstupuji hyperparametry
        % jen pro prvni parametr (neovlivni to sice odhady jako takove,
        %ale ovlivnilo by to vypocet marginalni verohodnosti)

fprintf('\n\n Omezeny model pro beta = 1 \n');
fprintf('Parametr    Prior   Prior std.    Posterior   Posterior std.\n');
fprintf('============================================================\n');
fprintf('Alpha       %6.4f   %6.4f     %6.4f  %6.4f\n',res_GM_rest.beta_0(1),...
    res_GM_rest.b0_std(1),res_GM_rest.beta_1(1),res_GM_rest.b1_std(1));
fprintf('h           %6.4f   %6.4f     %6.4f  %6.4f\n',res_GM_rest.h_0,...
    res_GM_rest.h0_std,res_GM_rest.h_1,res_GM_rest.h1_std);

%logaritmus Bayesova faktoru porovnavajici model omezeny a neomezeny
log_BF = res_GM_rest.log_ML-res_GM.log_ML;
%Bayesuv faktor (odlogaritmujeme)
BF = exp(log_BF);

fprintf('\n Bayesuv faktor porovnavajici omezeny a neomezeny model\n');
fprintf('BF = %6.4f\n',BF);

%% 3. Hypoteza beta > 1
% Lze analyticky (z posteriorni marginalni hustoty pro beta = t-rozdeleni
%nebo simulacne h|y~G(h_1,nu_1) a beta|h,y ~ N(beta_1,h^-1*V_1)
MC = 100000; %pocet simulaci
beta_sim = zeros(2,MC); %predpriprava simulovanych vektoru parametru beta (po sloupcich)

for ii=1:MC
   h_sim = gamm_rnd_Koop(res_GM.h_1,res_GM.nu_1,1); % generujeme z Gamma rozdeleni
   beta_sim(:,ii) = norm_rnd(h_sim^-1*res_GM.V_1)+res_GM.beta_1; %generujeme z normalniho rozdeleni s vyuzitim funkce norm_rnd (LeSageho ekonometricky toolbox)
end

%Vypocet pravdepodobnosti beta>1
pr_beta = sum(beta_sim(2,:)>1)/MC; % zjistime v jakem podilu vzorku je simulovana hodnota beta > 1
fprintf('\n Pravdepodobnost beta > 1 \n');
fprintf('Pr = %6.4f \n',pr_beta);

%analyticky vypocet pravdepodobnosti
%funkce tcdf statistickeho toolboxu (alternativne lze vyuzit LeSageho
%econ. toolbox
%a) standardizace skalovaneho t-rozdeleni (p(beta|y)) pro druhy prvek
%vektoru parametru beta
zscore = (1-res_GM.beta_1(2))/res_GM.b1_std(2);
%b)vypocet odpovidajiciho kvantilu ze standardizovaneho centrovaneho
%t-rozdeleni
pr_beta_analyticky = 1-tcdf(zscore,res_GM.nu_1);
fprintf('Pr = %6.4f (analyticky) \n',pr_beta_analyticky);


