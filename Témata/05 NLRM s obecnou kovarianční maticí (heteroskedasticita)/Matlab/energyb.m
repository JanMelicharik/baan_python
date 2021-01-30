%% Pindyck,Rubinfeld (1997) - example 6.4 - heteroskedasticita ve zname podobe

clear
close all
clc

%% Nastaveni cest a dat
addpath('..\Support'); %slozka je v nadrazenem adresari, proto ..\

tic %meric casu - pro zjisteni rychlosti programu

load EX64.dat

date = EX64(:,1); %rok 1960-1985
Q = EX64(:,2); %log mnoztvi distribuovane energie
P = EX64(:,3); %log ceny (za jednotku, v dolarech 1975)
INC = EX64(:,4); %duchod na jednu domacnost (v dolarech 1975)

n = length(Q);

YEAR = (1:n)';


%% Odhady Gibbsovym vzorkovacem za predpokladu zname heteroskedasticity

y = Q;
X = [ones(size(y)) YEAR P log(INC)];

[n,k] = size(X);

%prior hyperparameters
beta_0 = [3.5 0.1 -0.5 0.5]';
h_0 = 1/0.1^2;
nu_0 = 5;
var_beta_0 = [18^2 0.1^2 0.1^2 0.25^2]';
var_h_0 = 2*h_0^2/nu_0;
V_0 = diag(var_beta_0);

S0 = 15000+1;
S1 = 15000;
S = S0+S1;

%promenne "zpusobujici" heteroskedasticitu
z = [YEAR P log(INC)];
p = size(z,2);

beta = zeros(k,S);
h = zeros(1,S);
alpha_rw = zeros(p,S);
count_rw = 0; %èítaè akcepotvaných kandidátù

%po prvnich 50000 replikacich se vzal vysledny vektor rozptylu (a kovarianci) jako zaklad
%kovariancni matice; po te se dale ladi "d" pro ziskani zadouci
%akceptacni pravdepodobnosti;
%nebo - najde se modus posteriorni podminene hustoty pro alfa a prislusny
%hessian (podmineno apriornimi strednimi hodnotami, nebo jeste lepe
%posteriornimi strednimi hodnotami z prvnich 5000 replikaci

%postvar = eye(3);
%postvar = diag([0.132; 15.6423; 4.6199]); 
postvar = [   0.0775   -0.7646    0.3380
   -0.7646   16.5090   -8.1914
    0.3380   -8.1914    4.1386
          ];
d=0.1;
vscale_rw=d*postvar;
alpha_0 = [0.1278
    6.4524
   -3.5437];
%alpha_0 = [-0.1707; -6.2956; 3.1309];
%alpha_0= [1.4918; -7.4384; 4.3193];
Omega = diag((ones(n,1)+z*alpha_0).^2); %kov. matice Omega
invOmega = inv(Omega);

h(1) = h_0;
beta(:,1) = beta_0;
alpha_rw(:,1) = alpha_0;

error_check = 0; %pocitadlo chyb pri Gibbsove vzorkovaci (uprava kovariancni matice)

%% Gibbsuv vzorkovac
 fprintf('Posteriorni simulace\n')
 fprintf('0%%     50%%      100%%\n')
 fprintf('|')
 %pomocna promenna pro graficke znazorneni prubehu
 pom_step = 0.05; %0.05 = posun po 5% 
 pom_graph = pom_step;   
 
for s=2:S
 Omega = diag((ones(n,1)+z*alpha_rw(:,s-1)).^2); %kov. matice Omega
 invOmega = inv(Omega);
 b_Omega =inv(X'*invOmega*X)*X'*invOmega*y;
 V_1 = inv(inv(V_0)+h(:,s-1)*X'*invOmega*X);
 beta_1 = V_1*(inv(V_0)*beta_0+h(:,s-1)*X'*invOmega*X*b_Omega);
 %(nesystematicky zasah -- nekdy nastava numericky problem, ze Matlab nevezme V1 jako pozitivne
 %definitni matici => ignorujeme chybu a "udelame z V1 pozitivne definitni
 %prictenim maleho cisla 
  try
   beta(:,s)=beta_1+norm_rnd(V_1);
  catch
   beta(:,s)=beta_1+norm_rnd(V_1+0.000000001);
   error_check = error_check+1;
  end
 
  nu_1 =n+nu_0;
  h_1 = (1/nu_1*((y-X*beta(:,s))'*invOmega*(y-X*beta(:,s))+nu_0*h_0^-1))^-1;
  h(:,s)=gamm_rnd_Koop(h_1,nu_1,1);

%R-W M-H pro alpha
    a_can_rw = alpha_rw(:,s-1) + norm_rnd(vscale_rw); %kandidat
    log_accept_rw = min(a_post(a_can_rw,beta(:,s),h(:,s),y,X,z)...
        -a_post(alpha_rw(:,s-1),beta(:,s),h(:,s),y,X,z),0);
    if log_accept_rw > log(rand)
        alpha_rw(:,s)=a_can_rw;
        count_rw=count_rw+1;
    else
        alpha_rw(:,s)=alpha_rw(:,s-1);
    end   

    %graficke zobrazeni prubehu po 5 %
    if s/S>=pom_graph    
     fprintf('|');
     pom_graph = pom_graph+pom_step;
    end
    
 end  
 fprintf('\n')
 
%% vyhozeni prvnich S0 vzorku a pocitani posteriornich momentu
beta(:,1:S0) = [];
h(:,1:S0) = [];
alpha_rw(:,1:S0) = [];

beta_mean = mean(beta,2); %posteriorni stredni hodnota 
h_mean = mean(h);
alpha_mean = mean(alpha_rw,2);
var_beta_1 = mean(beta.^2,2)-beta_mean.^2; %posteriorni rozptyl
var_h_1 = mean(h.^2)-h_mean.^2;
var_alpha_1 = mean(alpha_rw.^2,2)-alpha_mean.^2;

accept_ratio_rw = count_rw/S; %prumerna akceptace


%KONVERGENCNI DIAGNOSTIKY - Gewekova CD pomoci funkce momentg
theta = [beta' h' alpha_rw'];
res_converg = Geweke(theta);

%HPDI pro beta a dlasi parametry
 hperc=0.95; %HPDI
 HPDI_beta = (quantile(beta',[0.05,0.95]))';
 HPDI_h = quantile(h,[0.05,0.95]);
 if p == 1 
  HPDI_alpha = quantile(alpha_rw,[0.05,0.95]); %jednorozmerny vektor ... neni potreba transpozice
 else
  HPDI_alpha = (quantile(alpha_rw',[0.05,0.95]))'; %jednorozmerny vektor ... neni potreba transpozice
 end


t = toc;

%% Prezentace vysledku, grafy
fprintf('Prior and posterior means for beta (std. deviations in parentheses)\n');
fprintf('                Prior       Posterior      NSE         Geweke CD             %3.1f%% HPDI \n',hperc*100);
fprintf('   ------------------------------------\n');

figure

for ii=1:k
 fprintf('Beta%1.0f     %12.4f   %12.4f   %10.4f   %10.4f   [%10.4f, %10.4f]\n',[ii beta_0(ii) beta_mean(ii) res_converg.NSE(ii) res_converg.CD(ii) HPDI_beta(ii,:)]);
 fprintf('         (%12.4f) (%12.4f)\n',[sqrt(var_beta_0(ii)) sqrt(var_beta_1(ii))]);

 subplot(3,2,ii)
 hist(beta(ii,:),50)
 xlabel(['\beta_' num2str(ii)])
end

 fprintf('h         %12.4f   %12.4f   %10.4f   %10.4f   [%10.4f, %10.4f]\n',[h_0 h_mean res_converg.NSE(k+1) res_converg.CD(k+1) HPDI_h]);
 fprintf('         (%12.4f) (%12.4f)\n', sqrt(var_h_0), sqrt(var_h_1));
 figure
  hist(h,50)
  xlabel('h')

figure
for ii=1:p
 fprintf('Alpha%1.0f         neinf.    %12.4f   %10.4f   %10.4f   [%10.4f, %10.4f]\n',[ii alpha_mean(ii) res_converg.NSE(k+1+ii) res_converg.CD(k+1+ii)  HPDI_alpha(ii,:)]);
 fprintf('                (N.A.)  (%12.4f)\n',sqrt(var_alpha_1(ii)));

    subplot(3,2,ii)
    hist(alpha_rw(ii,:),50)
    xlabel(['\alpha_' num2str(ii)])
end

fprintf('   ------------------------------------\n');
fprintf('Number of replications: %12.4f total \n',S-1);
fprintf('                        %12.4f discarded \n',S0-1);
fprintf('   ------------------------------------\n');

fprintf('   ------------------------------------\n');
fprintf('Accept.Ratio:  %8.4f\n', accept_ratio_rw)
fprintf('   ------------------------------------\n');

fprintf('   ------------------------------------\n');
fprintf('Errors:        %4.0f\n', error_check)
fprintf('   ------------------------------------\n');

fprintf('Time elapsed (s): %12.2f \n',t);
fprintf('***************************************************************************\n\n\n');
