%% Ilustrace Gibbsova vzorkovaèe
%Cvièení 4, pøíklad 1
clc
clear
close all

%% Definování parametrù a poèáteèního nastavení G. vzorkovaèe
%stredni hodnoty
mu = [0
      0];
%kovariancni matice
rho = 0.5;
Sigma = [1 rho
         rho 1];
     
S = 10000; %poèet generovaných vzorkù
S1 = 5000; %poèet ponechaných vzorkù
S0 = S-S1; %poèet vyhozených vzorkù

%matice pro generované vzorky (matice samých nul)
theta = zeros(2,S+1);
%poèáteèní hodnota pro bìh Gibbsova vzorkovaèe = theta(0)
theta(:,1) = [-1000
              1000];
          
%% Gibbsùv vzorkovaè
for s=2:S+1
  %generovani theta_1/theta_2 ~ N(mu_12,Sigma_12)
   mu_12 = mu(1)+rho*(theta(2,s-1)-mu(2));
   Sigma_12 = 1-rho^2;
   theta(1,s) = randn*sqrt(Sigma_12)+mu_12;
  %generovani theta_2/theta_1 ~ N(mu_21,Sigma_21)
   mu_21 = mu(2)+rho*(theta(1,s)-mu(1));
   Sigma_21 = 1-rho^2;
   theta(2,s) = randn*sqrt(Sigma_21)+mu_21; 
end

%% Zobrazime prvnich k kroku Gibbsova vzorkovace
k = 10;
figure
 subplot(2,2,1)
 plot(theta(1,1:k),theta(2,1:k),'*-')
 xlabel('\theta_1')
 ylabel('\theta_2')
 title(['Prvnich ',num2str(k),' kroku Gibbsova vzorkovace'])
 
%% Vyhodime prvnich S0+1 vzorku a vypocet statistik sdruzeneho rozdeleni
theta(:,1:S0+1) = [];
theta_mean = mean(theta,2); %vypocet stredni hodnoty pres radky
theta_cov = cov(theta'); %kovariancni matice

disp('Vektor strednich hodnota')
disp(theta_mean)

disp('Kovariancni matice')
disp(theta_cov)

%% Graficke zobrazeni posteriorni hustoty
 subplot(2,2,2)
 plot(theta(1,:),theta(2,:),'.')
 xlabel('\theta_1')
 ylabel('\theta_2')
 title(['Sdruzena hustota na zaklade ',num2str(S1),' vzorku'])
 
%graficke overeni konvergence (vykresleni kazdeho k-teho vzorku)
 k = 100;
 subplot(2,2,3)
 plot(theta(1,1:k:end))
 ylabel('\theta_1')
 title(['Konvergence na zaklade ',num2str(k),'-te replikace'])

 subplot(2,2,4)
 plot(theta(2,1:k:end))
 ylabel('\theta_2')
 title(['Konvergence na zaklade ',num2str(k),'-te replikace'])


     
  