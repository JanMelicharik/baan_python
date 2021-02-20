%% Tobit model - data o postu ziskanych medaily (souhrnny panel)
 clear
 close all
 clc

% Pridani cesty k podpurnym funkcim v adresari Support
% adresar pro vlastni funkce (pokud nejsou ve stejnem adresari jako skript)
 addpath('.\Support')

% Nacteni data 
% structure array data_olympics, 1610 pozorovani
%{   
    country     country code
    year        olympics year
    gdp         gross domestic product, 1995 dollars
    pop         population
    gold        number of gold medals won
    silver      number of silver medals won
    bronze      number of bronze medals won
    medaltot    total number of medals won
    host        = 1 if host country
    planned     = 1 if non-soviet planned
    soviet      = 1 if soviet
    
Data source: Andrew B. Bernard and Meghan R. Busse "Who wins the olympic games: Economic resources and medal totals,"
             The Review of Economics and Statistics, February 2004, 86(1), 413-417

Prevzato z Hill et al. (2007) doplneno o promennou:

    share        = podil ziskanych medail v letech 1988, 1992 a 1996
    (share = medaltot/738.*(year==88)+medaltot/815.*(year==92)+medaltot/842.*(year==96)
%}

load olympics.mat

%% Priprava dat
% ziskani sloupcoveho vektoru ze structure array: [...]'
 y_raw = [data_olympics.share]';
 X_raw = [ones(size(y_raw)) log([data_olympics.gdp])' log([data_olympics.pop])'];

% ocisteni o chybejici hodnoty
 ind_nonan = ~isnan(y_raw);
 y = y_raw(ind_nonan);
 X = X_raw(ind_nonan,:);

%% Nastaveni apriornich hyperparametru a Gibbsova vzorkovace
% Apriornimi hyperparametry
% p(beta)~N(beta_0, V_0)
% p(h)~G(h_0,nu_0)
 beta_0 = [0
           0
           0
           ];
 V_0 = diag([0.2^2;0.01^2;0.01^2]);
 nu_0 = 100;
 s2_0 = 0.1^2;
 h_0 = 1/s2_0;

 % Nastaveni Gibbsova vzorkovace
 S = 50000+1;   %celkovy pocet generovanych vzorku + pocatecni hodnota
 S_0 = 30000+1; %pocet vyhozenych vzorku
 S_1 = S-S_0;   %pocet ponechanych vzorku

 beta = zeros(length(beta_0),S);    %vzorky pro beta
 h = zeros(1,S);                    %vzorky pro h
 y_ast = zeros(length(y),S);        %vzorky pro y_ast (latentni data)
  
% nastaveni pocatecnich hodnot
 beta(:,1) = beta_0;
 h(1,1) = h_0;
 y_ast(:,1) = y;
 

%% Gibbsuv vzorkovac

% graficky ukazatel postupu simulace a zapnuti pocitadla casu
 tic
    
  fprintf('\t Gibbsuv vzorkovac - tobit model \n\n') %\t vytvori v retezci tabulator
  fprintf('0%%      50%%       100%%\n')
  fprintf('********************\n')
  fprintf('%c',1421)
  pom_step = 0.05;       %0.05 = posun po 5% (krok pocitadla
  pom_graph = pom_step;  %pocitadlo prubehu simulace


for s=2:S
 %2. blok Gibbsova vzorkovace
 %podminena hustota p(beta|h,y_ast)~N(beta_1,V_1)
  V_1 = inv(inv(V_0)+h(1,s-1)*(X'*X));                %(4.4) dle Koop (2003)
  beta_1 = V_1*(inv(V_0)*beta_0+h(1,s-1)*(X'*y_ast(:,s-1))); %(4.5) dle Koop (2003)
  
  beta(:,s) = beta_1+norm_rnd(V_1);                   %(4.7) dle Koop (2003)
  
 %2. blok Gibbsova vzorkovace
 %podminena hustota p(h|beta,y_ast)~G(h_1,nu_1)
  nu_1 = length(y)+nu_0;                              %(4.9)
  h_1 = (1/nu_1*(
     (y_ast(:,s-1)-X*beta(:,s))'*
     (y_ast(:,s-1)-X*beta(:,s))
     +nu_0*1/h_0))^-1; 
                                                      %(4.10)
 
  h(1,s) = gamm_rnd_Koop(h_1,nu_1,1);                 %(4.8)
  
 %3. blok Gibbsova vzorkovace
 %podminena hustoty pro p(y_ast|beta,h,y)
  y_ast(:,s) = y;
  %nalezeni y==0;
   ind = find(y==0);
  %generovani nah. cisel z prava omezeneho normalniho rozdeleni
   c_r = zeros(length(ind),1); %omezeni - nuly z prava
   y_ast(ind,s)=normrt_rnd(X(ind,:)*beta(:,s),ones(length(ind),1)*1/h(s),c_r);
 
  
  %Graficky prubeh simulace
   if s/S>=pom_graph    
    fprintf('%c',1421)
    pom_graph = pom_graph+pom_step;
   end 
   
end

 cas = toc; %zaznamenani casu od pocatku spusteni "tic"
 fprintf('\n')
 fprintf('********************\n\n')
 fprintf('Gibbsuv vzorkovac byl dokoncen!\n')
 fprintf('============================\n\n')

%% Posteriorni analyza
% vyhozeni prvnich S_0 vzorku
 beta(:,1:S_0) = [];
 h(:,1:S_0) = [];
 y_ast(:,1:S_0) = [];

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
 fprintf('Parametr  prior m.   prior std.  post m.     post std.   CD\n')
 fprintf('===========================================================\n')
 for ii=1:length(beta_0)
 fprintf('Beta %1.0f \t  %6.4f \t %6.4f \t %6.4f \t %6.4f \t %6.4f\n',...
    ii,beta_0(ii),std_beta_0(ii),mean_beta_1(ii),...
    std_beta_1(ii),CD_beta.CD(ii)) 
 end
 fprintf('h         %6.4f \t %6.4f \t %6.4f \t %6.4f \t %6.4f\n',...
    h_0,std_h_0,mean_h_1,std_h_1,CD_h.CD)
