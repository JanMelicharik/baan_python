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


%% Odhady Gibbsovym vzorkovacem, bez predpokladu heteroskedasticity

y = Q;
X = [ones(size(y)) YEAR P log(INC)];

[n,k] = size(X);

%% prior hyperparameters
beta_0 = [-30 0.1 -0.5 0.5]';


h_0 = 1/0.1^2;
nu_0 = 5;
var_beta_0 = [18^2 0.1^2 0.1^2 0.25^2]';
var_h_0 = 2*h_0^2/nu_0;
V_0 = diag(var_beta_0);

S0 = 15000+1;
S1 = 15000;
S = S0+S1;

beta = zeros(k,S);
h = zeros(1,S);
SD_nom = zeros(k,S); %citatel pro BF pomoci Savage-Dickey

%% Gibbsuv vzorkovac
h(1) = h_0;
beta(:,1) = beta_0;

nu_1 =n+nu_0; %v ramci cyklu se jiz nemeni

 fprintf('Posteriorni simulace\n')
 fprintf('0%%     50%%      100%%\n')
 fprintf('|')
 %pomocna promenna pro graficke znazorneni prubehu
 pom_step = 0.05; %0.05 = posun po 5% 
 pom_graph = pom_step;   
 
for s=2:S
   V_1 = inv(inv(V_0)+h(1,s-1)*X'*X);
   beta_1 = V_1*(inv(V_0)*beta_0+h(1,s-1)*X'*y);
    %b_Gibbs(:,i)=mvnrnd(b1,V1)';
   beta(:,s)=beta_1+norm_rnd(V_1); %spolehlivejsi generator!!!
   h_1 = (1/nu_1*((y-X*beta(:,s))'*(y-X*beta(:,s))+nu_0*h_0^-1))^(-1);
   h(:,s)=gamm_rnd_Koop(h_1,nu_1,1);
    
   %citatel pro Savage-Dickey ratio M1: beta_j=0
    for j = 1:k
      SD_nom(j,s) = norm_pdf(0,beta_1(j,1),V_1(j,j)); 
    end

    %graficke zobrazeni prubehu po 5 %
    if s/S>=pom_graph    
     fprintf('|');
     pom_graph = pom_graph+pom_step;
    end
 end  
 fprintf('\n')

%% vyhozeni prvnich S0 vzorku
beta(:,1:S0) = [];
h(:,1:S0) = [];
SD_nom(:,1:S0)=[];

b_mean = mean(beta,2); %posteriorni stredni hodnota 
h_mean = mean(h);
b_var = mean(beta.^2,2)-b_mean.^2; %posteriorni rozptyl
h_var = mean(h.^2)-h_mean.^2;

%% Savage-Dickey ratio (Bayes factor)
%jmenovatel pro Savage-Dickey ratio M1: beta_j=0
SD_denom = zeros(k,1);
for j = 1:k
    SD_denom(j,1) = norm_pdf(0,beta_0(j,1),V_0(j,j));
end
BF = mean(SD_nom,2)./SD_denom; 

%% KONVERGENCNI DIAGNOSTIKY - Gewekova CD pomoci funkce Geweke.m
theta = [beta' h'];
res_converg = Geweke(theta);

t = toc;

%% Prezentace vysledku, grafy
fprintf('Prior and posterior means for beta (std. deviations in parentheses)\n');
fprintf('                Prior       Posterior      NSE         Geweke CD     Post. Odds\n');
fprintf('                                                                     (beta_j=0)\n');
fprintf('   ------------------------------------\n');

figure
for ii=1:k
 fprintf('Beta%1.0f     %12.4f   %12.4f   %10.4f  %10.4f   %10.4f\n',[ii beta_0(ii) b_mean(ii) res_converg.NSE(ii) res_converg.CD(ii) BF(ii)]);
 fprintf('         (%12.4f) (%12.4f)\n',[sqrt(var_beta_0(ii)) sqrt(b_var(ii))]); 

 subplot(3,2,ii)
 hist(beta(ii,:),50)
 xlabel(['\beta_' num2str(ii)])

end

fprintf('h           %12.4f   %12.4f  %10.4f  %10.4f  \n',[h_0 h_mean res_converg.NSE(k+1) res_converg.CD(k+1)]);
fprintf('           (%12.4f) (%12.4f)\n',sqrt(var_h_0),sqrt(h_var));
    figure
    hist(h,50)
    xlabel('h')



fprintf('   ------------------------------------\n');
fprintf('Number of replications: %12.4f total \n',S-1);
fprintf('                        %12.4f discarded \n',S0-1);
fprintf('   ------------------------------------\n');
fprintf('Time elapsed (s): %12.2f \n',t);
fprintf('***************************************************************************\n\n\n');
