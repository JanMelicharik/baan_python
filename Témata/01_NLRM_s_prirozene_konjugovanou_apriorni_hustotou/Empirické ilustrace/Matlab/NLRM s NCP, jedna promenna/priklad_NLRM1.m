clear all
close all
clc

%Empiricka ilustrace - kapitola 2 - Koop 2003

%simulace umelych dat pro empirickou ilustraci kapitoly 2 z Koop (2003) 
n=50;                       %pocet generovanych "pozorovani"
beta=2;                     %parametr beta
sigma=1;                    %smerodatna odchylka
error= sigma*randn(n,1);    %tvorba chyboveho clenu - N(0,sigma^2)
x=rand(n,1);               %tvorba vysvetlujici promenne - U(0,1)
y=beta*x+error;             %tvorba vysvetlovane promenne

%data odpovidajici vysledkum z ucebniho textu (chceme-li vlastni, originalni
%vysledky, nutno zakomentovat)
load data_NLRM1.mat; %obsahuje nagenerovane promenne x a y


%% Hyperparametry pro prirozene konjugovany prior
nu0=10;
b0=1.5;
s02=1;
h0 = 1/s02;     %apriorni presnost chyby 'h'
V0=.25;


%Data pro vykresleni apriorniho marginalniho rozdeleni pro beta (t-rozdeleni) 
alower=1;                                  %dolni hranice intervalu (osy x)
aupper=3;                                  %horni hranice intervalu (osy x)
incr=0.01;                                 %krok pro generovani "osy x"
bplot=alower:incr:aupper;                  %generovani bodu na ose x, pro ktere spocteme hustotu pravdepodobnosti
b_plotpri = my_tpdf(bplot,b0,s02*V0,nu0);  %hustoty pravdepodobnosti, vyuziti vlastni funkce pro funkci hustoty t-rozdeleni


%Posteriorni analyza
%Odhady metodou nejmensich ctvercu
bols = inv(x'*x)*x'*y;      %normalni rovnice
nu=n-1;
s2 = (y-x*bols)'*(y-x*bols)/nu;

%Posteriorni hyperparametry pro normalni-gama rozdeleni
xsquare=x'*x;
nu1=nu0+n;      %posteriorni pocet stupnu volnosti
V1 = 1/(inv(V0)+xsquare); 
b1 = V1*(inv(V0)*b0 + bols*xsquare);
nu1s12 = nu0*s02 + nu*s2 + ((bols-b0)^2)/(V0+(1/xsquare));
s12 = nu1s12/nu1;

bcov = V1*nu1s12/(nu1-2); %rozptyl (obecne kovariancni matice) pro posteriorni beta
bstd=sqrt(bcov); %smerodatna odchylka pro posteriorni beta

%Posteriorni stredni hodnota a rozptyl presnosti chyby h
hmean = 1/s12;
hvar=2/(nu1*s12^2);
hstd=sqrt(hvar);

%Predikce
xstar=.5; %predpovidana hodnota x
ystarmean = xstar*b1;
ystarV = (1+V1*xstar^2)*s12;
ystarvar = ystarV*nu1/(nu1-2);
ystarstd=sqrt(ystarvar);

%Logaritmus marginalni verohodnosti pro informativni prior
    intcon=gammaln(nu1/2) + .5*nu0*log(nu0*s02)- gammaln(.5*nu0) -.5*n*log(pi); %vztah pro integracni konstantu
    lmarglik=intcon + .5*log(V1/V0) - .5*nu1*log(nu1s12); %logaritmus marginalni verohodnosti

    lmarg1=lmarglik;

%Data pro vykresleni posteriorniho marginalniho rozdeleni pro beta (t-rozdeleni)
b1_plotpost = my_tpdf(bplot,b1,s12*V1,nu1);


%Vypis vysledku
fprintf('Posterior and Predikce pro Informativni Prior\n');
fprintf('                    Prior Beta  Posterior Beta    Prior h   Posterior h\n');
fprintf('Stredni hodnota        %6.3f      %6.3f          %6.3f     %6.3f\n',[b0 b1 1/s02  hmean]);
fprintf('Smerodatna odchylka    %6.3f      %6.3f          %6.3f     %6.3f\n',[sqrt((nu0*s02/(nu0-2)*V0)) bstd sqrt(2/(nu0*s02^2)) hstd])
fprintf('***************************************************************************\n');
fprintf('                     Predikce pro x=0.5\n');
fprintf('E[y*]      %6.3f   \n',ystarmean)
fprintf('std[y*]    %6.3f   \n',ystarstd)
fprintf('Log of Marginal likelihood:    %6.3f\n\n\n',lmarg1)




%% Hyperparametry pro neinformativni prior
%nu0=0;
%inv(V0)=0;

%Posteriorni analyza (zkopirovaní predchoziho postupu + uprava dle prioru
%Posteriorni hyperparametry pro normalni-gama rozdeleni
xsquare=x'*x;
nu1=n;
V1 = 1/xsquare;
b1 = bols;
nu1s12 = nu*s2;
s12 = nu1s12/nu1;

bcov = V1*nu1s12/(nu1-2); 
bstd=sqrt(bcov); %smerodatna odchylka pro posteriorni beta

%Posteriorni stredni hodnota a rozptyl presnosti chyby h
hmean = 1/s12;
hvar=2/(nu1*s12^2);
hstd=sqrt(hvar);

%Predikce
xstar=.5; %predpovidana hodnota x
ystarmean = xstar*b1;
ystarV = (1+V1*xstar^2)*s12;
ystarvar = ystarV*nu1/(nu1-2);
ystarstd=sqrt(ystarvar);
b1_plotpost_noninf = my_tpdf(bplot,b1,s12*V1,nu1);




%Vypis vysledku
fprintf('Posterior and Predikce pro Neinformativni Prior\n');
fprintf('                      Posterior Beta   Posterior h\n');
fprintf('Stredni hodnota          %6.3f          %6.3f\n',[b1 hmean]);
fprintf('Smerodatna odchylka      %6.3f          %6.3f\n',[bstd hstd])
fprintf('***************************************************************************\n');
fprintf('                     Predikce pro x=0.5\n');
fprintf('E[y*]      %6.3f   \n',ystarmean)
fprintf('std[y*]    %6.3f   \n\n\n',ystarstd)


figure
plot(bplot,b_plotpri./sum(b_plotpri),'--',bplot,b1_plotpost./sum(b1_plotpost),bplot,b1_plotpost_noninf./sum(b1_plotpost_noninf),':')
title('Figure 2.1: Marginální Prior a Posteriory pro \beta')
legend('Prior','Posterior','Likelihood')
xlabel('\beta')
ylabel('Hustota pravdepodobnosti')


%% Porovnani s modelem obsahujicim jen urovnovou konstantu
x=ones(n,1);

%Hyperparametry - stejne jako u neomezeneho modelu v uvodu

%posteriorni analyza
%Odhady metodou nejmensich ctvercu
bols = inv(x'*x)*x'*y;
nu=n-1;
s2 = (y-x*bols)'*(y-x*bols)/nu;
bolscov = s2*inv(x'*x);
bolssd=sqrt(bolscov);

%Posteriorni hyperparametry pro normalni-gama rozdeleni
xsquare=x'*x;
nu1=nu0+n;
V1 = 1/(inv(V0)+xsquare);
b1 = V1*(inv(V0)*b0 + bols*xsquare);
nu1s12 = nu0*s02 + nu*s2 + ((bols-b0)^2)/(V0+(1/xsquare));
s12 = nu1s12/nu1;

bcov = V1*nu1s12/(nu1-2); 
bstd=sqrt(bcov); %smerodatna odchylka pro posteriorni beta

%Posteriorni stredni hodnota a rozptyl presnosti chyby h
hmean = 1/s12;
hvar=2/(nu1*s12^2);
hstd=sqrt(hvar);

%Predikce
xstar=.5; %predpovidana hodnota x
ystarmean = xstar*b1;
ystarV = (1+V1*xstar^2)*s12;
ystarvar = ystarV*nu1/(nu1-2);
ystarstd=sqrt(ystarvar);

%Logaritmus marginalni verohodnosti pro informativni prior (pouze urovnova
%konstanta)
    intcon=gammaln(nu1/2) + .5*nu0*log(nu0*s02)- gammaln(.5*nu0) -.5*n*log(pi); %vztah pro integracni konstantu
    lmarglik=intcon + .5*log(V1/V0) - .5*nu1*log(nu1s12); %logaritmus marginalni verohodnosti

postodds=exp(lmarg1 - lmarglik); %podil sanci modelu jen s parametrem sklonu a jen s urovnovou konstantou


 
%Vypis vysledku
fprintf('Posterior and Predikce pro Informativni Prior a model pouze s urovnovou konstantou\n');
fprintf('                    Prior Beta  Posterior Beta    Prior h   Posterior h\n');
fprintf('Stredni hodnota        %6.3f      %6.3f          %6.3f     %6.3f\n',[b0 b1 1/s02  hmean]);
fprintf('Smerodatna odchylka    %6.3f      %6.3f          %6.3f     %6.3f\n',[sqrt((nu0*s02/(nu0-2)*V0)) bstd sqrt(2/(nu0*s02^2)) hstd])
fprintf('***************************************************************************\n');
fprintf('                     Predikce pro x=0.5\n');
fprintf('E[y*]      %6.3f   \n',ystarmean)
fprintf('std[y*]    %6.3f   \n',ystarstd)
fprintf('***************************************************************************\n');
fprintf('Log of Marginal likelihood:    %6.3f\n',lmarglik)
fprintf('***************************************************************************\n');
fprintf('Podil sanci:    %6.3f\n',postodds)







