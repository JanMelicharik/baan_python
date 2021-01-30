%% Cviceni 02, priklad 2
clear           %vymaze vsechny promenne z pameti
close all       %zavre vsechna okna (s obrazky)
clc             %vymaze command window (prikazove okno)

%% Definovani vektoru s pocty generovanych vzorku
MC = [10 100 10000]; %budeme generovat postupne 10, 100, a nakonec 10 000 vzorku

%% a) MC simulace z N(1,4)
% Pomocny vypis na obrazovku (umozni nam sledovat postup vyhodnocovani kodu), \n zpusobi, ze kurzor se posune na dalsi radek
fprintf('Simulace rozdeleni N(1,4)\n')
fprintf('=========================\n')
figure %otevre okno pro kresleni
for s = 1:length(MC) % cyklus pro jednotlive prvky vektoru MC
    %generovani vzorku z N(1,4)
    % randn() je funkce pro generovani nahodnych vyberu ze standardizovaneho normalniho rozdeleni
    smpl = 1+2*randn(MC(s),1); %rozptyl = 4
                               %sm. odchylka = 2
    %pocitani statistik
    E = mean(smpl);     %vyberovy prumer jako odhad stredni hodnoty
    D = mean(smpl.^2)-E^2; %vyberovy rozptyl
    
    %vystup na obrazovku (formatovany), vypise pocet simulovanych vzorku, 
    %spocteny vyberovy prumer, rozptyl, a smerodatnou odchylku
    fprintf('MC = %8.0f E = %6.4f D = %6.4f sd = %6.4f\n',...
        MC(s),E,D,sqrt(D))
    % %8.0f ... rezervuje 8 znaku a 0 desetinnych mist pro cislo

    %vykresleni rozdeleni pomoci histogramu
    subplot(3,1,s)
    histogram(smpl)
    title(['MC = ',num2str(MC(s))])
end

%% b) MC simulace z U(2,5)
fprintf('\n\n')
fprintf('Simulace rozdeleni U(2,5)\n')
fprintf('=========================\n')
figure %otevre okno pro kresleni
for s = 1:length(MC)
    %generovani vzorku z U(2,5)
    % rand() je funkce pro generovani nahodnych vyberu z uniformniho rozdeleni
    smpl = 2+(5-2)*rand(MC(s),1); 
    %pocitani statistik
    E = mean(smpl);     %vyberovy prumer jako odhad stredni hodnoty
    D = mean(smpl.^2)-E^2;  %vyberovy rozptyl
    fprintf('MC = %8.0f E = %6.4f D = %6.4f sd = %6.4f\n',...
        MC(s),E,D,sqrt(D))
    subplot(3,1,s)
    histogram(smpl)
    title(['MC = ',num2str(MC(s))])
end

%% c) MC simulace z G(2,10) ... dle Koopa
fprintf('\n\n')
fprintf('Simulace rozdeleni G(2,10)\n')
fprintf('=========================\n')
figure %otevre okno pro kresleni
for s = 1:length(MC)
    %generovani vzorku z G(2,10)
    % gamrnd_Koop() je nami definovana funkce pro generovani nahodnych
    % vyberu z Gamma rozdeleni dle definice v Koop(2003)
    smpl = gamrnd_Koop(2,10,MC(s),1); 
    %pocitani statistik
    E = mean(smpl);     %vyberovy prumer jako odhad stredni hodnoty
    D = mean(smpl.^2)-E^2; %vyberovy rozptyl
    fprintf('MC = %8.0f E = %6.4f D = %6.4f sd = %6.4f\n',...
        MC(s),E,D,sqrt(D))
    subplot(3,1,s)
    histogram(smpl)
    title(['MC = ',num2str(MC(s))])
end
