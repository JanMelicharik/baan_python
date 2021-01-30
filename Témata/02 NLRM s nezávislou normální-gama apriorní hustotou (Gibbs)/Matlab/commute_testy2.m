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
 
%% Dalsi deklarace vektoru, napr. pro Savage-Dickey pomer hustot
% Jmenovatel pro SD pomer hustot - priklad d)
SDd_nom = zeros(1,S);
% Jmenovatel pro SD pomer hustot - priklad e)
SDe_nom = zeros(1,S);

count_post_restrict = 0; %pocitadlo nesplneni podminek apriornich
                         %restrikci na parametry

%% Gibbsuv vzorkovac
for s=2:S
 %1. blok Gibbsova vzorkovace
 %podminena hustota p(beta|h,y)~N(beta_1,V_1)
  V_1 = inv(inv(V_0)+h(1,s-1)*(X'*X)); %(4.4) dle Koop (2003)
  beta_1 = V_1*(inv(V_0)*beta_0+h(1,s-1)*(X'*y)); %(4.5) dle Koop (2003)
  
  check_restrict = 0; %nastaveni promenne overujici splneni
                      %podminky apriornich restrikci 
                      %na hodnoty parametru (0 = nesplneno)
  while check_restrict==0 %cyklus se opakuje dokud neni splnena
                          %apriorni restrikce na parametry
    beta(:,s) = beta_1+norm_rnd(V_1); %(4.7) dle Koop (2003)
    count_post_restrict = count_post_restrict+1;
    if prod(beta(:,s)>=0)==1 %soucin vyhodnoceni podminek kladnosti
        %na parametry musi byt roven 1 (podminky lze menit dle
        %potreby dle apriornich omezeni - nutno je v ramci
        %vypoctu S-D pomeru hustot zohlednit v apriornich hustotach
       check_restrict = 1; 
       count_post_restrict = count_post_restrict-1;
        %korekce pocitadla
    end
  end
 %2. blok Gibbsova vzorkovace
 %podminena hustota p(h|beta,y)~G(h_1,nu_1)
  nu_1 = length(y)+nu_0;            %(4.9)
  h_1 = (1/nu_1*((y-X*beta(:,s))'*(y-X*beta(:,s))+nu_0*1/h_0))^-1; %(4.10)
 
  h(1,s) = gamm_rnd_Koop(h_1,nu_1,1); %(4.8)
  
 %Cast pro vypocet citatelu SD pomer hustot
 %d) beta_depart = 0 (druhy prvek vektoru beta
 % p(b|h,y) = N(beta_1,V_1)
  V_1 = inv(inv(V_0)+h(1,s)*(X'*X));
  beta_1 = V_1*(inv(V_0)*beta_0+h(1,s)*(X'*y));
  SDd_nom(s) = normpdf(0,beta_1(2),sqrt(V_1(2,2)));
 %e) beta_reds = beta_trains
 % p(R*b|h,y) = N(R*beta_1,R*V_1*R')
  R = [0 0 1 -1];
  r = 0;
  SDe_nom(s) = mvnpdf(r,R*beta_1,R*V_1*R');    
end

%% Posteriorni analyza
% vyhozeni prvnich S_0 vzorku
 beta(:,1:S_0) = [];
 h(:,1:S_0) = [];
 SDd_nom(:,1:S_0) = [];
 SDe_nom(:,1:S_0) = [];
 
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
    
%% test hypotezy ze beta_red (beta_2) >=2
%pro ilustraci: zobrazeni rozdeleni beta_2 pomoci histogramu
 figure
 histogram(beta(3,:))
 title('Posterní hustota parametru \beta_2')
 xlabel('\beta_2')
%pravdepodobnost beta_2>=2
 pc = sum(beta(3,:)>=2)/S_1;
 fprintf('\n Pravdepodobnost beta_red >= 2 a odpovidajici Bayesuv faktor\n')
 fprintf('\n Pravdepodobnost = %4.3f     Bayesuv faktor = %6.3f\n',...
     pc,pc/(1-pc))

%% d) Test hypotezy beta_depart = 0 (druhy prvek)radek vektoru/matice parametru beta)
%Savage-Dickey pomer hustot
%Vyhodnotime p(beta_depart=0|M2)
%apriorni hustota p(beta_depart)~N(beta_0(2),V_0(2,2))
%jmenovatel SD pomeru hustot
%vyuzijeme funkci statistickeho toolboxu normpdf 
%pro pdf normalniho rozdeleni
 SDd_denom = normpdf(0,beta_0(2),sqrt(V_0(2,2)));
%simulace integracni konstanty omezen apriorni hustoty
 S_sim = 50000;
 count_prior_restrict = 0; %pocitadlo nesplneni apriornich omezene
 for ii=1:S_sim
    check_restrict = 0;
    while check_restrict==0
     pom = beta_0(2)+norm_rnd(V_0(2,2));
     count_prior_restrict = count_prior_restrict+1;
     if pom>=0
        check_restrict=1;
        count_prior_restrict = count_prior_restrict-1;
     end
    end
 end
%integracni konstanta apriorni hustoty zohlednujici apriorni omezeni
int_prior_restrict = 1/(1-count_prior_restrict/...
    (S_sim+count_prior_restrict));
%korekce jmenovatele S-D pomeru hustot
 SDd_denom = int_prior_restrict*SDd_denom;
%vypocet citatele na zaklade vystupu Gibbsova vzorkovace (doplneneho)
 E_SDd_nom = mean(SDd_nom);
%integracni konstanta posteriorni hustoty zohlednujici
%apriorni omezeni
 E_SDd_nom = 1/(1-count_post_restrict/(S+count_post_restrict-1))...
     *E_SDd_nom;
%Bayesuv faktor
 SD_d = E_SDd_nom/SDd_denom;
fprintf('\n Bayesuv faktor pro beta_depart = 0 \n')
fprintf('BF = %6.4f\n',SD_d)

%% e) Test beta_red = beta_trains (treti a ctvrty prvek vektoru beta)
% Linearni omezeni R*beta = r
%matice omezeni (na cely vektor parametru beta
 R = [0 0 1 -1];
 r = 0;
%jmenovatel pro SD pomer hustot
%apriorni hustota p(beta)~N(beta_0,V_0)
%apriorni hustota p(R*beta)~N(R*beta_0,R*V_0*R')
%vyuziti funkce hustoty pro vicerozmerne normalni rozdeleni
%ze statistickeho toolboxu - mvnpdf
 SDe_denom = mvnpdf(r,R*beta_0,R*V_0*R');  

%simulace integracni konstanty omezen apriorni hustoty
 S_sim = 50000;
 count_prior_restrict = 0; %pocitadlo nesplneni apriornich omezene
 for ii=1:S_sim
    check_restrict = 0;
    while check_restrict==0
     pom = beta_0+norm_rnd(V_0);
     count_prior_restrict = count_prior_restrict+1;
     if prod(pom>=0)==1
        check_restrict=1;
        count_prior_restrict = count_prior_restrict-1;
     end
    end
 end
%integracni konstanta apriorni hustoty zohlednujici apriorni omezeni
int_prior_restrict = 1/(1-count_prior_restrict/...
    (S_sim+count_prior_restrict));
%korekce jmenovatele S-D pomeru hustot
 SDe_denom = int_prior_restrict*SDe_denom; 
  
%vypocet citatele na zaklade vystupu Gibbsova vzorkovace (doplneneho)
 E_SDe_nom = mean(SDe_nom);
%korekce citatele (neomezeneho normalniho rozdelenio) skrze
%integracni konstantu zohlednujici apriorni rozdeleni
%neni potreba ji pocitat znovu (spocitana v ramci Gibbsova vzorkovace)
 E_SDe_nom = 1/(1-count_post_restrict/(S+count_post_restrict-1))...
     *E_SDe_nom;
%Bayesuv faktor
 SD_e = E_SDe_nom/SDe_denom;
fprintf('\n Bayesuv faktor pro beta_reds = beta_trains \n')
fprintf('BF = %6.4f\n',SD_e)

%% Vypocet intervalu nejvyssi posteriorni hustoty
HPDI_90 = quantile([beta' h'],[0.05 0.95]);
fprintf('\n 90 procentni HPDI \n')
par_names = {'beta_0' 'beta_depart' 'beta_reds' 'beta_trains' 'h'};
for ii=1:length(par_names)
 fprintf('%15s %6.4f %6.4f\n',par_names{ii},HPDI_90(:,ii))   
end

