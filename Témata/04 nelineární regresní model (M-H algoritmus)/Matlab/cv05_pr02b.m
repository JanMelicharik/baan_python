clear
close all
clc

%% Nastaveni Random-Walk Chain M-H algoritmu
S = 50000+1; %celkovy pocet generovanych vzorku + pocatecni vzorek
S0 = 30000+1; %celkovy pocet vyhozenych vzorku
S1 = S-S0; %pocet ponechanych vzorku k dalsi analyze

theta = zeros(1,S); %radkovy vektor pro uchovavani vzorku
theta(1) = 0; %pocatecni vzorek - theta(0)
count = 0; %pocitadlo poctu akceptovanych vzorku

% kandidatska hustota N(theta(s-1),c^2)
c = 5;

%% Random Walk Chain M-H algoritmus
for s=2:S
   %a) generovani kandidata theta_ast
   theta_ast = randn*c+theta(s-1);
   
   %b) spocitame akceptacni pravdepodobnost
   alpha = min(exp(-1/2*abs(theta_ast) ...
               +1/2*abs(theta(s-1))),1);
   
   %c) rozhodnuti o akceptaci nebo zamitnuti kandidata
   if alpha > rand
       theta(s) = theta_ast;
       count = count+1;
   else
       theta(s) = theta(s-1);
   end
end

%% Prezentace vysledku
theta(1:S0) = []; %vyhozeni prvnich S0 vzorku

%graficke zobrazeni konvergence
figure
plot(theta(1:500:end))

%pocitani momentu
E_theta = mean(theta);
D_theta = mean(theta.^2)-E_theta^2;
avg_count = count/(S-1); %prumerna mira akceptace

%zobrazeni vysledku na obrazovku
fprintf('E_theta   D_theta  Prumerna akceptace \n')
fprintf('%6.4f    %6.4f   %6.4f\n', E_theta, D_theta,avg_count)

%Histogram rozdeleni
figure
histogram(theta,50)





