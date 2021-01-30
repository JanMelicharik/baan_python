%% Ilustrace Gibbsova vzorkova�e
%Cvi�en� 4, p��klad 1
clc
clear
close all

%% Definov�n� parametr� a po��te�n�ho nastaven� G. vzorkova�e
%stredni hodnoty
mu = [0
      0];
%kovariancni matice
rho = 0.5;
Sigma = [1 rho
         rho 1];
     
S = 10000; %po�et generovan�ch vzork�
S1 = 5000; %po�et ponechan�ch vzork�
S0 = S-S1; %po�et vyhozen�ch vzork�

%matice pro generovan� vzorky (matice sam�ch nul)
theta = zeros(2,S+1);
%po��te�n� hodnota pro b�h Gibbsova vzorkova�e = theta(0)
theta(:,1) = [-1000
              1000];
          
%% Gibbs�v vzorkova�
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


     
  