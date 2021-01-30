clear 
close all
clc

%% Nastaveni cest a dat
addpath('.\Support'); %slozka je v nadrazenem adresari, proto ..\

%% Data
%OBS      GPA      TUCE     PSI      GRADE 
data = [
  1      2.66      20      0        0  
  2      2.89      22      0        0  
  3      3.28      24      0        0  
  4      2.92      12      0        0  
  5      4.00      21      0        1  
  6      2.86      17      0        0  
  7      2.76      17      0        0  
  8      2.87      21      0        0  
  9      3.03      25      0        0  
 10      3.92      29      0        1  
 11      2.63      20      0        0  
 12      3.32      23      0        0  
 13      3.57      23      0        0  
 14      3.26      25      0        1  
 15      3.53      26      0        0  
 16      2.74      19      0        0  
 17      2.75      25      0        0  
 18      2.83      19      0        0  
 19      3.12      23      1        0  
 20      3.16      25      1        1  
 21      2.06      22      1        0  
 22      3.62      28      1        1  
 23      2.89      14      1        0  
 24      3.51      26      1        0  
 25      3.54      24      1        1  
 26      2.83      27      1        1  
 27      3.39      17      1        1  
 28      2.67      24      1        0  
 29      3.65      21      1        1  
 30      4.00      23      1        1  
 31      3.10      21      1        0  
 32      2.39      19      1        1
 ];

%Greene (2002) - Table F21.1: Program Effectiveness, 32 Cross Section Observations
%Source: Spector and Mazzeo (1980).
%   * Obs = observation,
%   * TUCE = Test score on economics test,
%   * PSI = participation in program,
%   * GRADE = Grade increase (1) or decrease (0) indicator 
%   * GPA = grade point average

%The data listed in Appendix Table F21.1 were taken from a study by Spector and Mazzeo
%(1980), which examined whether a new method of teaching economics, the Personalized
%System of Instruction (PSI), significantly influenced performance in later economics courses.
%The “dependent variable” used in our application is GRADE, which indicates the whether
%a student’s grade in an intermediate macroeconomics course was higher than that in the
%principles course. The other variables are GPA, their grade point average; TUCE, the score
%on a pretest that indicates entering knowledge of the material; and PSI, the binary variable
%indicator of whether the student was exposed to the new teaching method. (Spector and
%Mazzeo’s specific equation was somewhat different from the one estimated here.)

%model:
%GRADE = b1 + b2*GPA + b3*TUCE + b4*PSI + eps

GPA   = data(:,2);
TUCE  = data(:,3);
PSI   = data(:,4);
GRADE = data(:,5);


%%

k = 4;      %pocet parametru
N = 32;     %poct pozorovani (32)

y = GRADE;
X = [ones(size(y)) GPA TUCE PSI];

%prior hyperparameters
beta_0 = [0 0 0 0]';
nu_0 = 2;

%pro identifikovatelnost parametru volime primo h=1
h = 1;

var_beta_0 = [10^2;10^2;10^2;10^2]';
V_0 = diag(var_beta_0);

y_ast_0 = y;

S0 = 5000+1; %Vyhozene replikace
S1 = 10000;
S = S0+S1; %Replikace

beta = zeros(k,S);
y_ast = zeros(N,S);


tic
%% Gibbs sampler (h=1)
I1 = find(y==1);
I2 = find(y==0);
   y_ast(I1,1)=normlt_rnd(X(I1,:)*beta(:,1),ones(length(I1),1),zeros(length(I1),1)); %nah. cislo z omezeneho norm. rozdeleni - omezeneho z leva nulou
   y_ast(I2,1)=normrt_rnd(X(I2,:)*beta(:,1),ones(length(I2),1),zeros(length(I2),1)); %nah. cislo z omezeneho norm. rozdeleni - omezeneho z prava nulou


 fprintf('Posterior simulation\n')
 fprintf('0%%     50%%      100%%\n')
 fprintf('|')
 %pomocna promenna pro graficke znazorneni prubehu
 pom_step = 0.05; %0.05 = posun po 5% 
 pom_graph = pom_step;  
   
 
for s=2:S
V_1 = inv(inv(V_0)+X'*X);
b_1 = V_1*(inv(V_0)*beta_0+X'*y_ast(:,s-1));

beta(:,s)=b_1+norm_rnd(V_1);


%{
%generovani latentnich dat (Koop 9.6)
for j=1:N
   if y(j)==1
   y_ast(j,i)=normlt_rnd(X(j,:)*beta(:,i),1,0); %nah. cislo z omezeneho norm. rozdeleni - omezeneho z leva nulou
   end
   if y(j)==0
   y_ast(j,i)=normrt_rnd(X(j,:)*beta(:,i),1,0); %nah. cislo z omezeneho norm. rozdeleni - omezeneho z prava nulou
   end
end
%}

   y_ast(I1,s)=normlt_rnd(X(I1,:)*beta(:,s),ones(length(I1),1),zeros(length(I1),1)); %nah. cislo z omezeneho norm. rozdeleni - omezeneho z leva nulou
   y_ast(I2,s)=normrt_rnd(X(I2,:)*beta(:,s),ones(length(I2),1),zeros(length(I2),1)); %nah. cislo z omezeneho norm. rozdeleni - omezeneho z prava nulou

    if s/S>=pom_graph    
     fprintf('|');
     pom_graph = pom_graph+pom_step;
    end
end
fprintf('\n');
fprintf('Done!\n');

%discard first S0 replications
beta(:,1:S0)=[];
y_ast(:,1:S0)=[];
 

beta_mean = mean(beta,2); %posterior mean 
y_ast_mean = mean(y_ast);
var_beta = mean(beta.^2,2)-beta_mean.^2; %posterior variance
var_y_ast = mean(y_ast.^2)-y_ast_mean.^2;


t = toc;

%Prezentace vysledku, grafy
fprintf('Prior and posterior means for beta (std. deviations in parentheses)\n');
fprintf('                Prior       Posterior      \n');
fprintf('                                           \n');
fprintf('   ------------------------------------\n');

for s=1:k
fprintf('Beta%1.0f     %12.4f   %12.4f   \n',[s beta_0(s) beta_mean(s)]);
fprintf('           (%12.4f) (%12.4f)\n',[sqrt(var_beta_0(s)) sqrt(var_beta(s))]);
end

fprintf('   ------------------------------------\n');
fprintf('Number of replications: %12.4f total \n',S-1);
fprintf('                        %12.4f discarded \n',S0-1);
fprintf('   ------------------------------------\n');
fprintf('Time elapsed (s): %12.2f \n',t);
fprintf('***************************************************************************\n\n\n');    




%% Pravdepodobnosti volby - Koop 9.7
%Vypocet pravdepodobnosti (predikce) zlepseni Pr(GRADE=1) na zvolenych
%charakteristikach
%GPA = variabilni, TUCE = prumer vzorku, PSI = variantne 0, 1 (nezarazen
%nebo zarazen do programu)

%Varianta I - pocitani se strednimi hodnotami posteriorni hustoty parametru
%- b_mean

GPA_pred = (2:0.01:4)'; %vypocet pro ruzne hodnoty predikovaneho GPA
nn = size(GPA_pred);
Pr_PSI0 = 1-norm_cdf(-[ones(nn) GPA_pred mean(TUCE)*ones(nn) 0*ones(nn)]*beta_mean);
Pr_PSI1 = 1-norm_cdf(-[ones(nn) GPA_pred mean(TUCE)*ones(nn) 1*ones(nn)]*beta_mean);

%mezni efekt ucasti v programu (rozdil pravdepodobnosti pri danem GPA pri
%PSI=0 a PSI=1)
margPSI = Pr_PSI1-Pr_PSI0;

%graficke zobrazeni - dovedene k dokonalosti :-)
figure
[AX,H1,H2] = plotyy(GPA_pred,[Pr_PSI0 Pr_PSI1],GPA_pred,margPSI,'plot');
set(get(AX(1),'Ylabel'),'String','Pr(GRADE=1) = efekt zlepseni','fontsize',12,'fontweight','b') 
set(get(AX(2),'Ylabel'),'String','Mezni efekt ucasti v programu','fontsize',12,'fontweight','b')

set(H1,'LineWidth',2)
set(H2,'LineWidth',1.5,'LineStyle','--')

xlabel('GPA (score)','fontsize',12,'fontweight','b')
title('Pravdepodobnosti zlepseni pri ucasti/neucasti v programu a odpovidajici mezni efekty','fontsize',12,'fontweight','b')
grid

set(AX(1),'FontSize',12)
set(AX(2),'FontSize',12)

text(2.9,0.57,'PSI = 1','fontsize',12,'fontweight','b')
text(3.6,0.27,'PSI = 0','fontsize',12,'fontweight','b')
text(3.4,0.58,'PSI - mezni efekt','fontsize',12,'fontweight','b')


%Varianta II - vyjadreni posteriorni hustoty pravdepodobnosti volby

GPA_pred = (2:0.01:4)'; %vypocet pro ruzne hodnoty predikovaneho GPA
nn = size(GPA_pred);

Pr_PSI0 = zeros(length(GPA_pred),S-S0);
Pr_PSI1 = zeros(length(GPA_pred),S-S0);


fprintf('Computing marginal effects of PSI \n');
fprintf('0%%       50%%      100%% \n');
fprintf('|');
 %pomocna promenna pro graficke znazorneni prubehu
 pom_step = 0.05; %0.05 = posun po 5% 
 pom_graph = pom_step;  
   
for s=1:S1
Pr_PSI0(:,s) = 1-norm_cdf(-[ones(nn) GPA_pred mean(TUCE)*ones(nn) 0*ones(nn)]*beta(:,s));
Pr_PSI1(:,s) = 1-norm_cdf(-[ones(nn) GPA_pred mean(TUCE)*ones(nn) 1*ones(nn)]*beta(:,s));

   if s/S1>=pom_graph    
     fprintf('|');
     pom_graph = pom_graph+pom_step;
    end
end
fprintf('\n');
fprintf('Done!\n');

%stredni hodnoty a HPDI pro pravdepodobnosti zlepseni (pri danem GPA pri
%PSI=0 a PSI=1)
Pr_PSI0_mean = mean(Pr_PSI0,2);
Pr_PSI1_mean = mean(Pr_PSI1,2);

[k,l]=size(Pr_PSI0);
%HPDI pro PSI0 a PSI1
 hperc=0.95; %HPDI

HPDI_PSI0 = (quantile(Pr_PSI0',[1-hperc,hperc]))';%pro PSI0
HPDI_PSI1 = (quantile(Pr_PSI1',[1-hperc,hperc]))';%pro PSI0


%% graficke zobrazeni - opet dovedene k dokonalosti
figure
axes('FontSize',12)
hold on
plot(GPA_pred,Pr_PSI0_mean,'b','LineWidth',2)


hold on
plot(GPA_pred,Pr_PSI1_mean,'g','LineWidth',2)

hold on
plot(GPA_pred,HPDI_PSI0(1),'b:',GPA_pred,HPDI_PSI0(2),'b:',GPA_pred,HPDI_PSI1(1),'g:',GPA_pred,HPDI_PSI1(2),'g:')


ylabel('Pr(GRADE=1) = efekt zlepseni','fontsize',12,'fontweight','b') 
xlabel('GPA (score)','fontsize',12,'fontweight','b')
title('Pravdepodobnosti zlepseni pri ucasti/neucasti v programu','fontsize',12,'fontweight','b')
grid

text(2.9,0.57,'PSI = 1','fontsize',12,'fontweight','b')
text(3.6,0.27,'PSI = 0','fontsize',12,'fontweight','b')
    