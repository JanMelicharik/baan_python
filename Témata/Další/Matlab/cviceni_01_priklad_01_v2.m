%% Muj druhy skript - modifikace hry ALice s Bobem
clear           %vymaze vsechny promenne z pameti
close all       %zavre vsechna okna (s obrazky)
clc             %vymaze command window (prikazove okno)

%% 1. Hra Alice a Boba - deklarace promennych a parametru
S = 1000000;    %pocet opakovani cele hry (pocet simulaci)
r = 0;          %pocitadlo stavu 19:10 (t.j. kolikrat nastane stav 19:10 v prospech Alice)
E_Bwin = 0;     %stredni hodnota vyhry Boba pri stavu 19:10 (inicializace)

%% 2. Simulace hry
for s=1:S       %v cyklu pro 1,2,3,..., S opakujeme
  p = rand;     %pocatecni rozdeleni stolu
                %generujeme p~U(0,1)
  %nahodne zahrani 29 her
  y = rand(29,1); %vektor rozmeru 29x1 s nezavislymi U(0,1)

  %overime zda-li nastal stav 19:10 => Alice musi vyhrat prave 19x
  %jestli y < p ....vrati vektor 0/1, kde 1 = vyhra Alice
  %sum(y<p) ... sum secte prvky vektoru, t.j. vypocte pocet vyher Alice 
  %             -> porovname s 19 (jestli je po 29 hrach stav 19:10, pak  
  %             Alice vyhrala prave 19-krat)
  if (sum(y<p)==19)
      r = r+1; %pocitadlo stavu 19:10 se posune o jedna
      %scitame pravdepodobnosti vyhry Boba
      %(z predchozich behu pri stavu 19:10)
      E_Bwin = E_Bwin + (1-p)^10;
  end
end

%% Vypocet a prezentace pozadovanych statistik
%Stredni (ocekavana) hodnota pravdepodobnosti vyhry Boba pri
%stavu 19:10
 E_Bwin = E_Bwin/r;
%Stredni (ocekavana) hodnota pravdepodobnosti vyhry Alice pri
%stavu 19:10
 E_Awin = 1-E_Bwin;

%Jednoduchy vypis vysledku na obrazovku
disp('Prumerna pravdepodobnost vyhry Boba pri stavu 19:10')
disp(E_Bwin)
disp('Prumerna pravdepodobnost vyhry Alice pri stavu 19:10')
disp(E_Awin)
disp('Ferovy podil sanci A:B pro rozdeleni vyhry za stavu 19:10')
disp(E_Awin/E_Bwin)
disp('Pocet stavu 19:10')
disp(r)
disp('Relativni zastoupeni poctu stavu 19:10 na celkovem poctu simulaci')
disp(r/S)

%Graficke zobrazeni vysledku (vetsi znacky a uprava mezi os)
figure
plot(1,E_Bwin,'*r','LineWidth',3,'MarkerSize',16) % Stredni hodnota pravdepodobnosti vyhry Boba (jako bod s x-ovou souradnici=1)
hold on
plot(2,E_Awin,'+g','LineWidth',3,'MarkerSize',16) % Stredni hodnota pravdepodobnosti vyhry Alice (jako bod s x-ovou souradnici=2)
xlim([0,3])         % omezeni rozsahu osi x na interval (0,3)
ylim([-0.1,1.1])    % omezeni rozsahu osi y na interval (-0.1, 1.1)
title('Muj druhy graf')

