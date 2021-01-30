%% Pindyck,Rubinfeld (1997) - example 6.4

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


%% Odhady Gibbsovym vzorkovacem za predpokladu nezname heteroskedasticity

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

beta = zeros(k,S+1);
h = zeros(1,S+1);
lambda = zeros(n,S+1);
nu_lambda_rw = zeros(1,S+1);

count_rw = 0; %èítaè akceptovaných kandidátù

%RW M-H algoritmus pro nu_lambda (jedna promenna)
postvar = 3;
d=0.1;
vscale_rw=d*postvar;
nu_lambda_0 = 5; %apriorni stredni hodnota pro lambda p(nu_lambda~G(nu_lambda_0,2))

beta(:,1) = beta_0;
h(:,1) = h_0;
lambda(:,1) = ones(n,1); %kov. matice Omega - predpoklad lambda_0_i = 1
nu_lambda_rw(:,1)=nu_lambda_0;

error_check = 0; %pocitadlo chyb pri Gibbsove vzorkovaci (uprava kovariancni matice)

%% Gibbsuv vzorkovac
 fprintf('Posteriorni simulace\n')
 fprintf('0%%     50%%      100%%\n')
 fprintf('|')
 %pomocna promenna pro graficke znazorneni prubehu
 pom_step = 0.05; %0.05 = posun po 5% 
 pom_graph = pom_step;  
 
for s=2:S
    Omega = diag(lambda(:,s-1).^-1); %kov. matice Omega 
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
    h_1 = (1/nu_1*((y-X*beta(:,s))'*invOmega*(y-X*beta(:,s))+nu_0*h_0^(-1)))^(-1);
    h(:,s)=gamm_rnd_Koop(h_1,nu_1,1);

    epsilon = y-X*beta(:,s); %vypocet epsilonu (na zaklade beta)
    for ii=1:n
        nu_lambda_1 = nu_lambda_rw(:,s-1)+1;
        lambda_1 = nu_lambda_1/(h(:,s)*epsilon(ii)^2+nu_lambda_1-1);
        lambda(ii,s)=gamm_rnd_Koop_2(lambda_1,nu_lambda_1,1); %Koop (6.25)
        %pouzita spolehlivejsi funkce pro generovani nah. cisel z gamma
        %rozdeleni
    end

%R-W M-H pro nu_lambda
    nu_lambda_can_rw = nu_lambda_rw(:,s-1) + norm_rnd(vscale_rw); %kandidat
    if nu_lambda_can_rw <=0 %nulova akceptacni pravdepodobnost pro zaporne nebo nulove kandidaty!
        log_accept_rw =-inf;
    else
    log_accept_rw = min(nu_lambda_post(nu_lambda_can_rw,nu_lambda_0,lambda(:,s))-nu_lambda_post(nu_lambda_rw(:,s-1),nu_lambda_0,lambda(:,s)),0);
    end
    if log_accept_rw > log(rand)
        nu_lambda_rw(:,s)=nu_lambda_can_rw;
        count_rw=count_rw+1;
    else
        nu_lambda_rw(:,s)=nu_lambda_rw(:,s-1);
    end
    
    %graficke zobrazeni prubehu po 5 %
    if s/S>=pom_graph    
     fprintf('|');
     pom_graph = pom_graph+pom_step;
    end
 end  
 fprintf('\n')

%vyhozeni prvnich S0 vzorku
beta(:,1:S0) = [];
h(:,1:S0) = [];
lambda(:,1:S0) = [];
nu_lambda_rw(:,1:S0) = [];

%% Vypocet posteriornich charakteristik
b_mean = mean(beta,2); %posteriorni stredni hodnota 
h_mean = mean(h);
lambda_mean = mean(lambda,2);
nu_lambda_mean = mean(nu_lambda_rw);
b_var = mean(beta.^2,2)-b_mean.^2; %posteriorni rozptyl
h_var = mean(h.^2)-h_mean.^2;
lambda_var = mean(lambda.^2,2)-lambda_mean.^2; %posteriorni rozptyl
nu_lambda_var = mean(nu_lambda_rw.^2,2)-nu_lambda_mean.^2;

accept_ratio_rw = count_rw/S;

%% KONVERGENCNI DIAGNOSTIKY - Gewekova CD pomoci funkce momentg
theta = [beta' h' nu_lambda_rw'];
res_converg = Geweke(theta);

%HPDI pro beta a dalsi parametry
 hperc=0.95; %HPDI
 HPDI_beta = (quantile(beta',[0.05,0.95]))';
 HPDI_h = quantile(h,[0.05,0.95]); %jednorozmerny vektor ... neni potreba transpozice
 HPDI_nu_lambda = quantile(nu_lambda_rw,[0.05,0.95]); 
 HPDI_lambda = (quantile(lambda',[0.05,0.95]))';
 
t = toc;

%% Prezentace vysledku, grafy
fprintf('Prior and posterior means for beta (std. deviations in parentheses)\n');
fprintf('                Prior       Posterior      NSE         Geweke CD             %3.1f%% HPDI \n',hperc*100);
fprintf('   ------------------------------------\n');

figure
for ii=1:k
fprintf('Beta%1.0f     %12.4f   %12.4f   %10.4f   %10.4f   [%10.4f, %10.4f]\n',[ii beta_0(ii) b_mean(ii) res_converg.NSE(ii) res_converg.CD(ii) HPDI_beta(ii,:)]);
fprintf('         (%12.4f) (%12.4f)\n',[sqrt(var_beta_0(ii)) sqrt(b_var(ii))]);

    subplot(3,2,ii)
    hist(beta(ii,:),50)
    xlabel(['\beta_' num2str(ii)])
end

 fprintf('h         %12.4f   %12.4f   %10.4f   %10.4f   [%10.4f, %10.4f]\n',[h_0 h_mean res_converg.NSE(k+1) res_converg.CD(k+1) HPDI_h]);
 fprintf('         (%12.4f) (%12.4f)\n', sqrt(var_h_0), sqrt(h_var));
 figure
  hist(h,50)
  xlabel('h')

figure

fprintf('nu_lambda %12.4f   %12.4f   %10.4f   %10.4f   [%10.4f, %10.4f]\n',[nu_lambda_0 nu_lambda_mean res_converg.NSE(k+1+1) res_converg.CD(k+1+1) HPDI_nu_lambda]);
fprintf('         (%12.4f) (%12.4f)\n',[nu_lambda_0 sqrt(nu_lambda_var)]);

    hist(nu_lambda_rw,50)
    xlabel('\nu_\lambda')

figure
 boxplot(lambda')
 xlabel('n')
    

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
