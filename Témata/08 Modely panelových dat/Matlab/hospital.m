clear
close all
clc

%% Nastaveni cest a dat
addpath('.\Support');


load data_hospital.mat

%% Data
hospital_ID = data_hospital(:,1);   %... hospital ID
year        = data_hospital(:,2);   %... year
costs       = data_hospital(:,3);   %... log of hospital operating costs (in thousands of dollars)
beds        = data_hospital(:,4);   %... log of number of beds in hospital
inpatient   = data_hospital(:,5);   %... log of number of inpatients visit
outpatient  = data_hospital(:,6);   %... log of number of outpatient visit
case_mix    = data_hospital(:,7);   %... log of case mix index (higher values = difficult cases)
K           = data_hospital(:,8);   %... log of capital stock (in thousand dollars)
nonprofit   = data_hospital(:,9);   %... 1 for non-profit hospitals (=0 otherwise)
forprofit    = data_hospital(:,10); %... 1 for for-profit hospitals (=0 otherwise)


%% *********************************************

N = 382;  %pocet nemocnic
T = 5; %pocet let
k = 8;  %pocet parametru

%% Tvorba datovych matic

y = costs; %logaritmus nakladu
y_i = reshape(y,T,N); %y_i pro i=1...N (po sloupcich)

X = [ones(N*T,1) beds inpatient outpatient case_mix K nonprofit forprofit];
X_i = zeros(T,k,N); %X_i pro i=1...N (po treti dimenzi)
for s = 1:N
   X_i(:,:,s) = X((T*s-T)+1:T*s,:); 
end
X_i_til = X_i(:,2:end,:); %X_i s vlnovkou (bez urovnove konstanty)

X_ast = zeros(T*N,(N+k-1)); 
for s = 1:N
   X_ast((T*s-T)+1:T*s,s)=ones(T,1);
   X_ast((T*s-T)+1:T*s,end-(k-1)+1:end)=X_i_til(:,:,s);
end

tic

%% Pooled model
%"standardni" Gibbsuv vzorkovac - jen vetsi matice a vektory
%prior hyperparameters
    beta_0 = [1 1 1 1 1 1 1 1]';

    h_0 = 1/0.5^2;
    var_beta_0 = [1^2 1^2 1^2 1^2 1^2 1^2 1^2 1^2]';
    nu_0 = 1;
    var_h_0 = 2*h_0^2/nu_0;
    V_0 = diag(var_beta_0);

    S0 = 1000+1;
    S1 = 10000;
    S = S0+S1;

beta = zeros(k,S);
h = zeros(1,S);

%% Gibbsuv vzorkovac
h(1) = h_0;

 fprintf('Posteriorni simulace\n')
 fprintf('0%%     50%%      100%%\n')
 fprintf('|')
 %pomocna promenna pro graficke znazorneni prubehu
 pom_step = 0.05; %0.05 = posun po 5% 
 pom_graph = pom_step;  
 
nu_1 =T*N+nu_0;

for s=2:S

   V_1 = inv(inv(V_0)+h(1,s-1)*X'*X);
   beta_1 = V_1*(inv(V_0)*beta_0+h(1,s-1)*X'*y);
   beta(:,s)=beta_1+norm_rnd(V_1);  %podminena hustota pro beta
   h_1 = (1/nu_1*((y-X*beta(:,s))'*(y-X*beta(:,s))+nu_0*1/h_0))^-1;
   h(:,s)=gamm_rnd_Koop(h_1,nu_1,1); %podminena hustota pro h
 
   if s/S>=pom_graph    
     fprintf('|');
     pom_graph = pom_graph+pom_step;
    end
end

fprintf('\n');
fprintf('Hotovo!\n');


%vyhozeni prvnich S0 vzorku
    beta(:,1:S0) = [];
    h(:,1:S0) = [];

    beta_mean = mean(beta,2); %posteriorni stredni hodnota 
    h_mean = mean(h);
    var_beta = mean(beta.^2,2)-beta_mean.^2; %posteriorni rozptyl
    var_h = mean(h.^2)-h_mean.^2;

%KONVERGENCNI DIAGNOSTIKY - pro uspornost neresime
tt = toc;

%Prezentace vysledku, grafy
fprintf('Pooled model\n');
fprintf('Prior and posterior means for beta (std. deviations in parentheses)\n');
fprintf('                Prior       Posterior      \n');
fprintf('                                           \n');
fprintf('   ------------------------------------\n');

for s=1:k
    fprintf('Beta%1.0f     %12.4f   %12.4f   \n',[s beta_0(s) beta_mean(s)]);
    fprintf('           (%12.4f) (%12.4f)\n',[sqrt(var_beta_0(s)) sqrt(var_beta(s))]);
end

    fprintf('h         %12.4f   %12.4f   \n',[h_0 h_mean]);
    fprintf('           (%12.4f) (%12.4f)\n',[sqrt(var_h_0) sqrt(var_h)]);


fprintf('   ------------------------------------\n');
fprintf('Number of replications: %12.4f total \n',S-1);
fprintf('                        %12.4f discarded \n',S0-1);
fprintf('   ------------------------------------\n');
fprintf('Time elapsed (s): %12.2f \n',tt);
fprintf('***************************************************************************\n\n\n');











%% Individual effects model - nehierarchicka apriorni hustota
%prior hyperparameters
    alpha_0 = ones(N,1)';       %variabilni urovnove konstanty napric aerolinkami
    beta_0 = [alpha_0 1 1 1 1 1 1 1]';     %stejne ostatni parametry (sklonu)

    h_0 = 1/1^2;
    
    var_alpha_0 = ones(N,1)';
    var_beta_0 = [var_alpha_0 1^2 1^2 1^2 1^2 1^2 1^2 1^2]';
    nu_0 = 1;
    var_h_0 = 2*h_0^2/nu_0;
    V_0 = diag(var_beta_0);

    S0 = 500+1;
    S1 = 500;
    S = S0+S1;

beta = zeros(N+k-1,S);
h = zeros(1,S);

%% Gibbsuv vzorkovac
h(1) = h_0;
nu_1 =T*N+nu_0;

 fprintf('Posteriorni simulace\n')
 fprintf('0%%     50%%      100%%\n')
 fprintf('|')
 %pomocna promenna pro graficke znazorneni prubehu
 pom_step = 0.05; %0.05 = posun po 5% 
 pom_graph = pom_step;  
 
for s=2:S
   V_1 = inv(inv(V_0)+h(1,s-1)*X_ast'*X_ast);
   beta_1 = V_1*(inv(V_0)*beta_0+h(1,s-1)*X_ast'*y);
   beta(:,s)=beta_1+norm_rnd(V_1);  %podminena hustota pro beta
   pom = 0;
   for j=1:N
      pom = pom+(y_i(:,j)-beta(j,s)*ones(T,1)-X_i_til(:,:,j)*beta(N+1:N+(k-1),s))'*(y_i(:,j)-beta(j,s)*ones(T,1)-X_i_til(:,:,j)*beta(N+1:N+(k-1),s));
   end
   h_1 = (1/nu_1*(pom+nu_0*1/h_0))^-1;
   h(:,s)=gamm_rnd_Koop(h_1,nu_1,1); %podminena hustota pro h
   
    if s/S>=pom_graph    
     fprintf('|');
     pom_graph = pom_graph+pom_step;
    end
end

fprintf('\n');
fprintf('Hotovo!\n');

%% vyhozeni prvnich S0 vzorku
    beta(:,1:S0+1) = [];
    h(:,1:S0+1) = [];

    beta_mean = mean(beta,2); %posteriorni stredni hodnota 
    h_mean = mean(h);
    var_beta = mean(beta.^2,2)-beta_mean.^2; %posteriorni rozptyl
    var_h = mean(h.^2)-h_mean.^2;

%KONVERGENCNI DIAGNOSTIKY - pro uspornost neresim 
tt = toc;

%Prezentace vysledku, grafy
fprintf('Individual Effects Model - nehierarchicka apriorni hustota\n');
fprintf('Prior and posterior means for beta (std. deviations in parentheses)\n');
fprintf('                Prior       Posterior      \n');
fprintf('                                           \n');
fprintf('   ------------------------------------\n');

for s=1:N
    fprintf('Alfa%1.0f     %12.4f   %12.4f   \n',[s beta_0(s) beta_mean(s)]);
    fprintf('           (%12.4f) (%12.4f)\n',[sqrt(var_beta_0(s)) sqrt(var_beta(s))]);
end
for s=N+1:N+(k-1)
    fprintf('Beta%1.0f     %12.4f   %12.4f   \n',[s-N+1 beta_0(s) beta_mean(s)]);
    fprintf('           (%12.4f) (%12.4f)\n',[sqrt(var_beta_0(s)) sqrt(var_beta(s))]);
end

    fprintf('h         %12.4f   %12.4f   \n',[h_0 h_mean]);
    fprintf('           (%12.4f) (%12.4f)\n',[sqrt(var_h_0) sqrt(var_h)]);


fprintf('   ------------------------------------\n');
fprintf('Number of replications: %12.4f total \n',S-1);
fprintf('                        %12.4f discarded \n',S0-1);
fprintf('   ------------------------------------\n');
fprintf('Time elapsed (s): %12.2f \n',tt);
fprintf('***************************************************************************\n\n\n');







%% Individual effects model - hierarchicka apriorni hustota
%prior hyperparameters
    alpha_0 = ones(N,1)';       %variabilni urovnove konstanty napric aerolinkami
    beta_0 = [1 1 1 1 1 1 1]';     %stejne ostatni parametry (sklonu)

    h_0 = 1/1^2;
    
    var_alpha_0 = ones(N,1)';
    var_beta_0 = [1^2 1^2 1^2 1^2 1^2 1^2 1^2]';
    nu_0 = 1;
    var_h_0 = 2*h_0^2/nu_0;
    V_0 = diag(var_beta_0);

    mu_alpha_0 = 5;
    h_alpha_0 = 1/(2^2);
    V_alpha_0 = 5;
    nu_alpha_0 = 1;
    
    S0 = 500+1;
    S1 = 500;
    S = S0+S1;

beta = zeros(k-1,S);
h = zeros(1,S);
alpha = zeros(N,S);

V_alpha = zeros(1,S);
mu_alpha = zeros(1,S);

%% Gibbsuv vzorkovac
h(1) = h_0;
V_alpha(1) = V_alpha_0;
mu_alpha(1) = mu_alpha_0;
alpha(:,1) = alpha_0;

nu_1 =T*N+nu_0;
nu_alpha_1 =nu_alpha_0+N;
   
 fprintf('Posteriorni simulace\n')
 fprintf('0%%     50%%      100%%\n')
 fprintf('|')
 %pomocna promenna pro graficke znazorneni prubehu
 pom_step = 0.05; %0.05 = posun po 5% 
 pom_graph = pom_step;  
 
for s=2:S

pom = 0;
   for j=1:N
      pom = pom+X_i_til(:,:,j)'*X_i_til(:,:,j);
   end
V_1 = inv(inv(V_0)+h(:,s-1)*pom);

pom = 0;
   for j=1:N
      pom = pom+X_i_til(:,:,j)'*(y_i(:,j)-alpha(j,s-1)*ones(T,1));
   end
beta_1 = V_1*(inv(V_0)*beta_0+h(:,s-1)*pom);
beta(:,s)=beta_1+norm_rnd(V_1); %podminena hustota pro beta

   pom = 0;
   for j=1:N
      pom = pom+(y_i(:,j)-alpha(j,s-1)*ones(T,1)-X_i_til(:,:,j)*beta(:,s))'*(y_i(:,j)-alpha(j,s-1)*ones(T,1)-X_i_til(:,:,j)*beta(:,s));
   end
   h_1 = (1/nu_1*(pom+nu_0*1/h_0))^(-1);
h(:,s)=gamm_rnd_Koop(h_1,nu_1,1); %podminena hustota pro h

%podminene hustoty pro kazde alpha_j
for j=1:N
   V_1_j = (V_alpha(:,s-1)*h(:,s)^-1)/(T*V_alpha(:,s-1)+h(:,s)^-1);
   alpha_1 = (V_alpha(:,s-1)*(y_i(:,j)-X_i_til(:,:,j)*beta(:,s))'*ones(T,1)+h(:,s)^-1*mu_alpha(:,s-1))/(T*V_alpha(:,s-1)+h(:,s)^-1);
   alpha(j,s) = alpha_1+norm_rnd(V_1_j);  %Koop 7.16 
end

%podminena hustota pro mu_alfa
sigma_2_alpha_1 = V_alpha(:,s-1)*1/h_alpha_0/(V_alpha(:,s-1)+N*1/h_alpha_0);
mu_alpha_1 = (V_alpha(:,s-1)*mu_alpha_0+1/h_alpha_0*sum(alpha(:,s)))/(V_alpha(:,s-1)+N*1/h_alpha_0);
mu_alpha(:,s) = mu_alpha_1+norm_rnd(sigma_2_alpha_1); %Koop 7.17

%podminena hustota pro inv(V_alfa)
   pom = 0;
   for j=1:N
      pom = pom+(alpha(j,s)-mu_alpha(:,s))^2;
   end
V_alpha_1 = (pom+V_alpha_0*nu_alpha_0)/nu_alpha_1;
V_alpha(:,s)=inv(gamm_rnd_Koop(inv(V_alpha_1),nu_alpha_1,1)); %Koop 7.18


    if s/S>=pom_graph    
     fprintf('|');
     pom_graph = pom_graph+pom_step;
    end
end

fprintf('\n');
fprintf('Hotovo!\n');



%vyhozeni prvnich S0 vzorku
    alpha(:,1:S0) = [];
    beta(:,1:S0) = [];
    h(1:S0) = [];

    beta_mean = mean(beta,2); %posteriorni stredni hodnota 
    alpha_mean = mean(alpha,2);
    h_mean = mean(h);
    var_beta = mean(beta.^2,2)-beta_mean.^2; %posteriorni rozptyl
    var_alpha = mean(alpha.^2,2)-alpha_mean.^2;
    var_h = mean(h.^2)-h_mean.^2;

%KONVERGENCNI DIAGNOSTIKY - pro uspornost neresim 
tt = toc;

%Prezentace vysledku, grafy
fprintf('Individual Effects Model - hierarchicka apriorni hustota\n');
fprintf('Prior and posterior means for beta (std. deviations in parentheses)\n');
fprintf('                Prior       Posterior      \n');
fprintf('                                           \n');
fprintf('   ------------------------------------\n');

for s=1:N
    fprintf('Alfa%1.0f     %12.4f   %12.4f   \n',[s alpha_0(s) alpha_mean(s)]);
    fprintf('           (%12.4f) (%12.4f)\n',[sqrt(var_alpha_0(s)) sqrt(var_alpha(s))]);
end
for s=1:(k-1)
    fprintf('Beta%1.0f     %12.4f   %12.4f   \n',[s+1 beta_0(s) beta_mean(s)]);
    fprintf('           (%12.4f) (%12.4f)\n',[sqrt(var_beta_0(s)) sqrt(var_beta(s))]);
end

    fprintf('h         %12.4f   %12.4f   \n',[h_0 h_mean]);
    fprintf('           (%12.4f) (%12.4f)\n',[sqrt(var_h_0) sqrt(var_h)]);


fprintf('   ------------------------------------\n');
fprintf('Number of replications: %12.4f total \n',S-1);
fprintf('                        %12.4f discarded \n',S0-1);
fprintf('   ------------------------------------\n');
fprintf('Time elapsed (s): %12.2f \n',tt);
fprintf('***************************************************************************\n\n\n');
