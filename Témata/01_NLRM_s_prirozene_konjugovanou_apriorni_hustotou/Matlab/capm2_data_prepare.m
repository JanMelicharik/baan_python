%% Nacteni datoveho souboru capm2_data.mat, priprava datoveho souboru capm2.mat a vykresleni dat do grafu
clear           %vymaze vsechny promenne z pameti
close all       %zavre vsechna okna (s obrazky)
clc             %vymaze command window (prikazove okno)

%% Nacteni dat z excelovskeho souboru
 [NUM,TXT,RAW] = xlsread('capm2_data.xlsx');

%cell array obsahujici textove retezce nazvu promennych
 var_label = TXT(1,2:end);
 time_label = TXT(2:end,1);

 for ii=1:length(var_label)
     data.(var_label{ii}) = NUM(:,ii);
     %vytvoreni strukturni promenne 'data'
     %obsahujici data a nazvy
     %pro jednotlive promenne
 end

%ulozeni dat do datoveho souboru Matlabu 
 save capm2.mat var_label time_label data

%% Vykresleni vybrane promenne s vlastni osou x
%priklad obrazku s vlastnimi osami
 figure
 plot(data.MSFT,'LineWidth',2.5)
 index = 1:10:length(data.MSFT);
    %vytvoreni nove rady znacek na osi x (xticku)
 set(gca,'XTick',index)
    %prirazeni popisku xtickum
 set(gca,'XTickLabel',time_label(index))
    %otoceni popisek o 270 stupnu
 set(gca,'XTickLabelRotation',270) 
   
%% Obrazek s vyuzitim datovych typu Matlabu
startDate = datenum('01-01-1995'); %pocatecni datum
endDate = datenum(2005,1,12); %posledni datum

%vytvoreni datovych bodu 'xData' odpovidajici poctu mesicu
xData = linspace(startDate,endDate,length(data.MSFT));

figure
plot(xData,data.MSFT)
set(gca,'XTick',xData) %nastaveni poctu XTick 
                       %na pocet bodu xData
datetick('x',12)       %konverze oznaceni x-tick na mesice
xlim([startDate-30 endDate+30]); %uprava limitu osy x
