%% Muj prvni skript
clear           %vymaze vsechny promenne z pameti
close all       %zavre vsechna okna (s obrazky)
clc             %vymaze command window (prikazove okno)

%% 1. Hra Alice a Boba - deklarace promennych a parametru
S = 1000000;    %pocet opakovani cele hry (pocet simulaci)
r = 0;          %pocitadlo stavu 5:3 (t.j. kolikrat nastane stav 5:3 v prospech Alice)
E_Bwin = 0;     %stredni hodnota vyhry Boba pri stavu 5:3 (inicializace)

%% 2. Simulace hry
for s=1:S       %v cyklu pro 1,2,3,..., S opakujeme
  p = rand;     %pocatecni rozdeleni stolu
                %generujeme p~U(0,1)
  %nahodne zahrani 8 her
  y = rand(8,1); %vektor rozmeru 8x1 s nezavislymi U(0,1)
  
  %overime zda-li nastal stav 5:3 => Alice musi vyhrat prave 5x
  %jestli y < p ....vrati vektor 0/1, kde 1 = vyhra Alice
  %sum(y<p) ... sum secte prvky vektoru, t.j. vypocte pocet vyher Alice 
  %             -> porovname s 5 (jestli je po 8 hrach stav 5:3, pak Alice 
  %             vyhrala prave 5-krat)
  if (sum(y<p)==5)
      r = r+1; %pocitadlo stavu 5:3 se posune o jedna
      %scitame pravdepodobnosti vyhry Boba
      %(z predchozich behu pri stavu 5:3)
      E_Bwin = E_Bwin + (1-p)^3;
  end
end

%% Vypocet a prezentace pozadovanych statistik
%Stredni (ocekavana) hodnota pravdepodobnosti vyhry Boba pri
%stavu 5:3
 E_Bwin = E_Bwin/r;
%Stredni (ocekavana) hodnota pravdepodobnosti vyhry Alice pri
%stavu 5:3
 E_Awin = 1-E_Bwin;

%Jednoduchy vypis vysledku na obrazovku
disp('Prumerna pravdepodobnost vyhry Boba pri stavu 5:3')
disp(E_Bwin)
disp('Prumerna pravdepodobnost vyhry Alice pri stavu 5:3')
disp(E_Awin)
disp('Ferovy podil sanci A:B pro rozdeleni vyhry za stavu 5:3')
disp(E_Awin/E_Bwin)
disp('Pocet stavu 5:3')
disp(r)
disp('Relativni zastoupeni poctu stavu 5:3 na celkovem poctu simulaci')
disp(r/S)

%Graficke zobrazeni vysledku
figure
plot(1,E_Bwin,'*r') % Stredni hodnota pravdepodobnosti vyhry Boba (jako bod s x-ovou souradnici=1)
hold on
plot(2,E_Awin,'+g') % Stredni hodnota pravdepodobnosti vyhry Alice (jako bod s x-ovou souradnici=2)
xlim([0,3])     % omezeni rozsahu osi x na interval (0, 3)
title('Muj prvni graf')

