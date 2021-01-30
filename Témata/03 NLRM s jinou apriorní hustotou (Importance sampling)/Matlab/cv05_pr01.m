%Generovani momentu a vzorku N(0,1) z U(a,b)
 clear
 close all
 clc

%% Nastaveni cest
 addpath('.\Support')

%% Posteriorni hustotu N(0,1) + importance function U(a,b)
%nastaveni importance sampler
 S = 100000;
 w = zeros(1,S); %radkovy vektor vah
 theta = zeros(1,S); %vektor vzorku z importance function
 
%importance function U(a,b)
 a = -100;
 b = 100;

%% Importance sampling
for s = 1:S
   %1. generovani vzorku z importance function, q(theta)
   theta(s) = rand*(b-a)+a;
   %2. vypocet vah (na zaklade jadrovych hustot)
   p_theta = exp(-1/2*theta(s)^2); %jadrova hustota
   q_theta = 1/(b-a); %plna funkce hustoty (staci i jadro = 1)
   w(s) = p_theta/q_theta;    
end

%% Vypocet posteriorni stredni hodnoty a rozptylu
 E_theta = w*theta'/sum(w);
 E_theta_2 = w*(theta'.^2)/sum(w); %vazena E(theta^2)
 D_theta = E_theta_2-E_theta^2; %vazeny rozptyl
 w_avg = mean(w); %prumerna vaha

 fprintf('Importance sampling s vyuzitim U(%5.0f,%5.0f)\n',a,b)
 fprintf('E_theta   D_theta    w_avg          S\n')
 fprintf('%6.4f    %6.4f     %6.4f  %10.0f\n',E_theta,D_theta,w_avg,S)

%% Vazeny bootstrap pro ziskani vyberu z posteriorni hodnoty
 boot_smpl = my_boot(theta,length(theta),w);

figure
 subplot(2,1,1)
 histogram(theta,50)
 title('Vyber z importance function')
 
 subplot(2,1,2)
 histogram(theta(boot_smpl),50)
 title('Vyber z bootstrapovaneho rozdeleni')
 
% Vypocet strednich hodnot a rozptylu z bootstrapovaneho rozdeleni
 E_boot_theta = mean(theta(boot_smpl));
 D_boot_theta = var(theta(boot_smpl));
 
 fprintf('Momenty z boostrapovaneho rozdeleni\n')
 fprintf('E_boot_theta   D_boot_theta\n')
 fprintf('%6.4f           %6.4f     \n',E_boot_theta,D_boot_theta)

 

