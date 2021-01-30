clear
close all
clc

%% Pridani cesty k podpurnym funkcim v adresari Support
% adresar pro vlastni funkce (pokud nejsou ve stejnem adresari jako skript)
 addpath('E:\MUNI\projekt_bayes\Témata\Support\Matlab')

%% Nacteni dat
%pozorovani za 231 pracovnich dni v roce 2006)

% time      ... cas cesty
% depart    ... cas odchodu (v minutach od 6:30)
% reds      ... pocet cerevenych svetel na semaforu
% trains    ... pocet vlaku, ktere musi nechat projek na Murrumbeena prejezdu


load data_commute.mat

%% Nastaveni apriornich hyperparametru a Gibbsova vzorkovace
% Apriornimi hyperparametry
% p(beta)~N(beta_0, V_0)
% p(h)~G(h_0,nu_0)
 beta_0 = [30
           1
           3
           5];
 V_0 = diag([7.5^2;0.25^2;1^2;2^2]);
 nu_0 = 40;
 s2_0 = 10^2;
 h_0 = 1/s2_0;
 
% Definice modelu
 y = time;
 X = [ones(size(y)) depart reds trains];

% Nastaveni Gibbsova vzorkovace
 S = 50000+1;   %celkovy pocet generovanych vzorku + pocatecni hodnota
 S_0 = 30000+1; %pocet vyhozenych vzorku
 S_1 = S-S_0;   %pocet ponechanych vzorku

 beta = zeros(length(beta_0),S);    %vzorky pro beta
 h = zeros(1,S);                    %vzorky pro h
  
% nastaveni pocatecnich hodnot
 beta(:,1) = beta_0;
 h(1,1) = h_0;
 
% Dalsi deklarace vektoru, napr. pro Savage-Dickey pomer hustot
% ...

%% Gibbsuv vzorkovac
for s=2:S
 %1. blok Gibbsova vzorkovace
 %podminena hustota p(beta|h,y)~N(beta_1,V_1)
  V_1 = inv(inv(V_0)+h(1,s-1)*(X'*X)); %(4.4) dle Koop (2003)
  beta_1 = V_1*(inv(V_0)*beta_0+h(1,s-1)*(X'*y)); %(4.5) dle Koop (2003)
  
  beta(:,s) = beta_1+norm_rnd(V_1); %(4.7) dle Koop (2003)
  
 %2. blok Gibbsova vzorkovace
 %podminena hustota p(h|beta,y)~G(h_1,nu_1)
  nu_1 = length(y)+nu_0;            %(4.9)
  h_1 = (1/nu_1*((y-X*beta(:,s))'*(y-X*beta(:,s))+nu_0*1/h_0))^-1; %(4.10)
 
  h(1,s) = gamm_rnd_Koop(h_1,nu_1,1); %(4.8)
end

%% Posteriorni analyza
% vyhozeni prvnich S_0 vzorku
 beta(:,1:S_0) = [];
 h(:,1:S_0) = [];

 
% graficke zobrazeni konvergence
 k = 100;   %delka kroku
 figure
 for ii=1:length(beta_0)
    subplot(3,2,ii)
    plot(beta(ii,1:k:end));
    ylabel(['\beta_',num2str(ii)])
 end
  subplot(3,2,ii+1)
  plot(h(1:k:end))
  ylabel('h')

% Gewekova konvergencni diagnostika
 CD_beta = Geweke(beta');
 CD_h = Geweke(h');
 
%% Prezentace vysledku
%apriorni str. hodnoty a sm. odchylky
%beta_0, h_0 - apriorni stredni hodnoty
std_beta_0 = sqrt(diag(V_0)); %vektor apriornich sm. odchylek
std_h_0 = sqrt(2*h_0^2/nu_0); %apriorni sm. odchylka pro h

%posteriorni str. hodnoty a sm. odchylky
mean_beta_1 = mean(beta,2); %sloupcovy vektor radkovych prumeru
mean_h_1 = mean(h); %vystup = skalarni velicina
std_beta_1 = sqrt(mean(beta.^2,2)-mean_beta_1.^2);  
std_h_1 = sqrt(mean(h.^2)-mean_h_1.^2);

%Vystup na obrazovku
fprintf('Parametr  prior m.  prior std.  post m. post std. CD\n')
fprintf('====================================================\n')
for ii=1:length(beta_0)
fprintf('Beta %1.0f   %6.4f  %6.4f %6.4f %6.4f %6.4f\n',...
    ii,beta_0(ii),std_beta_0(ii),mean_beta_1(ii),...
    std_beta_1(ii),CD_beta.CD(ii)) 
end
fprintf('h        %6.4f  %6.4f %6.4f %6.4f %6.4f\n',...
    h_0,std_h_0,mean_h_1,std_h_1,CD_h.CD) 
    






