clear
close all
clc

data = [
    1955    114043  8310    182113
    1956    120410  8529    193749
    1957    129187  8738    205192
    1958    134705  8952    215130
    1959    139960  9171    225021
    1960    150511  9569    237026
    1961    157897  9527    248897
    1962    165286  9662    260661
    1963    178491  10334   275466
    1964    199457  10981   295378
    1965    212323  11746   315715
    1966    226977  11521   337642
    1967    241194  11540   363599
    1968    260881  12066   391847
    1969    277498  12297   422382
    1970    296530  12955   455049
    1971    306712  13338   484677
    1972    329030  13738   520553
    1973    354057  15924   561531
    1974    374977  14154   609825];

obs     = data(:,1); %pozorovani
GDP     = data(:,2); %GDP Mexika v milionech pesos roku 1960
Labor   = data(:,3); %v tis. obyvatel
Capital = data(:,4); %v milionech pesos roku 1960


%% Pridani cesty k support funkcim
addpath('.\Support')
    
%% Priprava dat (do podoby indexu)
y = GDP/(mean(GDP));    %index GDP, kde prumer = 1
%prevedeni promennych do podoby indexu
X = [ones(size(y)) Labor/mean(Labor) Capital/mean(Capital)];
N = length(y);

%% Apriorni hustoty a apriorni hyperparametry
%p(gamma) ~ N(gamma_0,V_0)
gamma_0 = [1
           0.5
           0.5
           1];
k = length(gamma_0); %pocet parametru modelu
V_0 =diag([0.5^2
          0.25^2
          0.25^2
          0.5^2]);

%p(h)~G(h_0,nu_0)
h_0 = 1/0.5^2;
nu_0 = 5;
      
%% Metropolis within Gibbs - nastaveni
 S = 50000+1; %+1 ... pocatecni hodnota
 S0 = 30000+1;
 S1 = S-S0; %pocet ponechanych vzorku

%ukladani vzorku
 gamma = zeros(k,S);
 gamma(:,1) = gamma_0;
 h = zeros(1,S);
 
%Nastaveni RW M-H algoritmu
%kandidatska hustota ~ N(gamma(s-1),Sigma)
 d = 0.5; %skalovaci konstanta
 Sigma = d*eye(k); %prvotni nastaveni kov. matice
                   %M-H algoritmu
 %Sigma = d*[
 %  0.0680   -0.0343   -0.0284   -0.0024
 %  -0.0343    0.0449   -0.0021    0.0037
 %  -0.0284   -0.0021    0.0341   -0.0015
 %  -0.0024    0.0037   -0.0015    0.2144
 %  ];
 count = 0; %pocitadlo akceptovanych vzorku
 
%% Metropolis within Gibbs - beh
for s=2:S
  %a) podminena hustota p(h|gamma,y)~G(h_1,nu_1)
  nu_1 = nu_0+N; %(5.23) z Koopa
  fx = gamma(1,s-1)*(gamma(2,s-1)*X(:,2).^gamma(4,s-1) ...
      +gamma(3,s-1)*X(:,3).^gamma(4,s-1)).^(1/gamma(4,s-1));
  h_1 = (1/nu_1*((y-fx)'*(y-fx)+nu_0*1/h_0))^-1; %(5.23)
  h(s) = gamm_rnd_Koop(h_1,nu_1,1); %vygenerovany vzorek
    
  %b) podminena hustota p(gamma|h,y) -> logaritmus
  %jadrove podminene hustoty = externi funkce log_post_CES.m
  %(5.24)
  %generovani kandidat z kandidatske hustoty
   gamma_ast = gamma(:,s-1)+norm_rnd(Sigma);
  %spocitame logaritmus akceptacni pravdepodobnosti
   log_accept = min(...
       log_post_CES(y,X,gamma_ast,h(s),gamma_0,V_0)...
       -log_post_CES(y,X,gamma(:,s-1),h(s),gamma_0,V_0),...
       0);
  %rozhodnuti o akceptaci
  if log_accept > log(rand)
      gamma(:,s) = gamma_ast;
      count = count+1;
  else
      gamma(:,s) = gamma(:,s-1);
  end    
end

%% Pomocne radky pro nastaveni matice Sigma
%nastavime S = 5000+1;
%disp(count/(S-1)); %zobrazeni prumerne akceptace
%disp(cov(gamma')); %odhad posteriorni kovariancni
%disp(mean(gamma,2)); %odhad posteriorni stredni hodnoty



%% Prezentace vysledku a posteriorni analyza
%a) vyhozeni S0 prvnich vzorku
 gamma(:,1:S0) = [];
 h(:,1:S0) = [];
 
%b) vypocet aposteriornich momentu (a apriornich)
  E_gamma = mean(gamma,2); %prumer po radcich
  D_gamma = var(gamma,1,2); %rozptyl (populacni) po radcich
  E_h = mean(h);
  D_h = var(h,1);
 
 %apriorni momenty
  E_gamma_0 = gamma_0;
  D_gamma_0 = diag(V_0);
  E_h_0 = h_0;
  D_h_0 = 2*E_h_0^2/nu_0;
%%
 %Vykresleni vysledku na obrazovku
 fprintf('\n Parametr  E_prior std_prior  E_post std_post\n')
 fprintf('===============================================\n')
 for ii=1:length(gamma_0)
   fprintf('gamma_%1.0f     %6.4f %6.4f %6.4f %6.4f\n',...
       ii,E_gamma_0(ii),sqrt(D_gamma_0(ii)),...
       E_gamma(ii),sqrt(D_gamma(ii)))      
 end
   fprintf('h           %6.4f %6.4f %6.4f %6.4f\n',...
       E_h_0,sqrt(D_h_0),E_h,sqrt(D_h))

   
%% Vypocet marginalni verohodnosti modelu metodou Gelfanda a Deye
%I) zjednodusena variantu
 GD_simple = 0;
for s = 1:S1
 GD_simple = GD_simple+1/S1*1/lik_CES(y,X,gamma(:,s),h(s));  
end

fprintf('\n\n')
fprintf('Marginalni verohodnost CES produkcni funkce\n')
fprintf('       (zjednodusena varianta)\n')
fprintf('marg_lik = %6.4f    log_marg_lik = %6.4f\n',...
    1/GD_simple,log(1/GD_simple))
%%
%II) plna varianta metody GD
 theta = [gamma;h]; %slouceni matic vzorku parametru do jedne matice
 theta_hat = mean(theta,2); %sloupcovy vektor strednich hodnot
 Sigma_hat = cov(theta');
 kk = length(theta_hat);
 pp = 0.01; %pro (1-p) procentni kvantil chi2 rozdeleni
 chi_pp = chis_inv(1-pp,kk);
 
 %simulace integracni konstanty pro omezenou apriorni hustotu
 count_g = 0;
 for s=1:S
   pom = gamma_0+norm_rnd(V_0); %p(gamma)~N(gamma_0,V_0)
   count_g = count_g+(min(pom)>0);
 end
 %integracni konstanta = 1/(podil splneni omezeni)
 int_c = 1/(count_g/S);
 
 GD = 0; %statistika Gelfanda a Deye
 for s=1:S1
    Theta = (theta(:,s)-theta_hat)'*inv(Sigma_hat)*...
        (theta(:,s)-theta_hat);
    %fce. hustoty p,eueneho vicerozmerneho rozdeleni
    f_theta = 1/(1-pp)*1/(2*pi)^(kk/2)*...
        det(Sigma_hat)^(-1/2)*...
        exp(-1/2*Theta)*(Theta<=chi_pp);
    prior = prior_CES(gamma(:,s),h(s),gamma_0,V_0,h_0,nu_0);
    like = lik_CES(y,X,gamma(:,s),h(s));
    GD = GD+f_theta/(int_c*prior*like)*1/S1;    
 end
 
fprintf('\n\n')
fprintf('Marginalni verohodnost CES produkcni funkce\n')
fprintf('       (plna varianta)\n')
fprintf('marg_lik = %6.4f    log_marg_lik = %6.4f\n',...
    1/GD,log(1/GD))


%% Srovnani skutecnych a modelovych momentu
E_y_ast = zeros(1,S1); %vektor modelovych strednich hodnot
D_y_ast = zeros(1,S1); %vektor modelovych rozptylu
std_y_ast = zeros(1,S1); %vektor modelovych sm. odchylek

%simulace (generovani) umelych dat
for s = 1:S1
   fx = gamma(1,s)*(gamma(2,s)*X(:,2).^gamma(4,s)...
       +gamma(3,s)*X(:,3).^gamma(4,s)).^(1/gamma(4,s));
   y_ast = fx+randn*sqrt(1/h(s));
   E_y_ast(1,s) = mean(y_ast);
   D_y_ast(1,s) = var(y_ast);
   std_y_ast(1,s) = std(y_ast);
end

%vypocet predikcnich (jednostrannych) p-hodnot
%a) pro stredni hodnotu
p_E = (sum(E_y_ast<mean(y))/S1)*(sum(E_y_ast<mean(y))/S1<=0.5)...
    +(1-sum(E_y_ast<mean(y))/S1)*(sum(E_y_ast<mean(y))/S1>0.5);
%b) pro rozptyl
p_D = (sum(D_y_ast<var(y))/S1)*(sum(D_y_ast<var(y))/S1<=0.5)...
    +(1-sum(D_y_ast<var(y))/S1)*(sum(D_y_ast<var(y))/S1>0.5);
%c) pro rozptyl
p_std = (sum(std_y_ast<std(y))/S1)*(sum(std_y_ast<std(y))/S1<=0.5)...
    +(1-sum(std_y_ast<std(y))/S1)*(sum(std_y_ast<std(y))/S1>0.5);

fprintf('\n Predikèní p-hodnoty \n');
fprintf('p_E = %4.3f p_D = %4.3f p_std = %4.3f\n',...
    p_E,p_D,p_std);

%Graficke zobrazeni predikcnich p-hodnot resp.
%simulovanych a datovych momentu
figure
 subplot(3,1,1)
 histogram(E_y_ast,50)
 hold on
 plot(mean(y),0,'*r','LineWidth',6)
 title('Simulovane stredni hodnoty')
 
 subplot(3,1,2)
 histogram(D_y_ast,50)
 hold on
 plot(var(y),0,'*r','LineWidth',6)
 title('Simulovane rozptyly')
 
 subplot(3,1,3)
 histogram(std_y_ast,50)
 hold on
 plot(std(y),0,'*r','LineWidth',6)
 title('Simulovane smerodatne odchylky')
