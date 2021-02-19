clear
close all
clc

%% Nastaveni cest a dat
addpath('..\Support'); %slozka je v nadrazenem adresari, proto ..\

load yankees.txt

yyear = yankees(:,1);
yPCT = yankees(:,2); %winning percentage in year t
yOBP = yankees(:,5); %team on-base percentage in year t
ySLG = yankees(:,6); %team slugging average in year t
yERA = yankees(:,8); %team earned run average in year t

load redsox.txt

ryear = redsox(:,1);
rPCT =  redsox(:,2); %winning percentage in year t
rOBP =  redsox(:,5); %team on-base percentage in year t
rSLG =  redsox(:,6); %team slugging average in year t
rERA =  redsox(:,8); %team earned run average in year t

tic

%% SUR model - tvorba struktury SUR modelu

N = length(yyear);
M = 2; %pocet rovnic;


yX = [ones(size(yyear)) yOBP ySLG yERA];
rX = [ones(size(ryear)) rOBP rSLG rERA];

yy = yPCT;
ry = rPCT;

k1 = size(yX,2);
k2 = size(rX,2);
k = k1+k2;

y = zeros(M*N,1);
X = zeros(M*N,k);
for ii=1:N
X(M*ii-1:M*ii,:) = [yX(ii,:) zeros(1,k2)
                  zeros(1,k1) rX(ii,:)];
y(M*ii-1:M*ii,1) = [yy(ii,1)
                  ry(ii,1)];
end




%% Apriorni hyperparametry a nastaveni Gibbsova vzorkovace

S0 = 15000+1; %po��te�n� podminka
S1 = 15000;
S = S0+S1;

beta_0 = zeros(k,1);
V_0 = 4*eye(k);
nu_0 = 0;
invH_0 = zeros(M);
H_0 = eye(M);

beta = zeros(k,S);
H = zeros(M,M,S); %trojrozmerna matice presnosti chyb (treti rozmer = beh Gibbsova vzorkovace)

beta(:,1) = beta_0;
H(:,:,1) = H_0;

%% Gibbsuv vzorkovac
 fprintf('Posteriorni simulace\n');
 fprintf('0%%       50%%      100%% \n');
 fprintf('|');
 %pomocna promenna pro graficke znazorneni prubehu
 pom_step = 0.05; %0.05 = posun po 5% 
 pom_graph = pom_step;
 



 
for s=2:S
 sumV1 = zeros(k,k);
 sumb1 = zeros(k,1);
 for ii=1:N
   sumV1 = sumV1+(X(M*ii-1:M*ii,:)')*H(:,:,s-1)*X(M*ii-1:M*ii,:); %cast sumy v Koop 6.50
   sumb1 = sumb1+(X(M*ii-1:M*ii,:)')*H(:,:,s-1)*y(M*ii-1:M*ii,1); %cast sumy v Koop 6.51)
 end

V_1 = inv(V_0)+sumV1; %neinvertovane V1 ... Koop 6.50 pro ucely matalbu pak nasledne deleni zleva (presnejsi vypocet inverze
beta_1 = V_1\(V_0\beta_0+sumb1); %inverze pomoci deleni zleva "\" ... Koop 6.51

%podminena hustota pro b
 beta(:,s) = beta_1+norm_rnd(inv(V_1)); %Koop 6.49

%podminena hustota pro H
 sumH = zeros(M,M);
 for ii=1:N
   sumH = sumH+ 
   (y(M*ii-1:M*ii,1)-X(M*ii-1:M*ii,:)*beta(:,s))
   *
   (y(M*ii-1:M*ii,1)-X(M*ii-1:M*ii,:)*beta(:,s))'; %cast sumy v Koop 6.54
 end
 nu_1 = N+nu_0; %Koop 6.52
 H_1 = inv(invH_0+sumH); %Koop 6.54

%podminena hustota pro H
 H(:,:,s) = wish_rnd(H_1,nu_1); %Koop 6.49





%graficke zobrazeni prubehu po 5 %
 if s/S>=pom_graph    
  fprintf('|');
  pom_graph = pom_graph+pom_step;
 end
end

%% vyhozeni prvnich S0 vzorku
beta(:,1:S0) = [];
H(:,:,1:S0) = [];
r = H(2,1,:)./(sqrt(H(1,1,:).*H(2,2,:))); %korelace nahodne slozky napric rovnicemi

b_mean = mean(beta,2); %posteriorni stredni hodnota 
H_mean = mean(H,3);
r_mean = mean(r);
b_var = mean(beta.^2,2)-b_mean.^2; %posteriorni rozptyl
H_var = mean(H.^2,3)-H_mean.^2; %prumerovani pres treti rozmer
r_var = mean(r.^2)-r_mean.^2;


%% HPDI pro beta a r (korelace mezi nahodnymi slozkami jednotlivych modelu)

hperc=0.95; %HPDI s vyuzitim funkce "quantile" (dostupna ve sdtatistickem i ekonometrickem toolboxu)
HPDI_beta = (quantile(beta',[0.05,0.95]))';
HPDI_r = quantile(r,[0.05,0.95]);

fprintf('\n');
fprintf('Hotovo!\n\n');

fprintf('Apriorni a posteriorni hustoty pro parametry beta a korelace nahodnych slozek\n');
fprintf('                  Prior       Posterior      %3.1f%% HPDI \n',hperc*100);
fprintf('   ------------------------------------\n');


fprintf('Rovnice pro Yankees\n');
fprintf('   ------------------------------------\n');
for kk=1:k1
fprintf('Beta%1.0f          %10.4f   %10.4f   [%10.4f, %10.4f]\n',[kk b_mean(kk) sqrt(b_var(kk)) HPDI_beta(kk,:)]);
end
fprintf('   ------------------------------------\n');
fprintf('Rovnice pro Red Sox\n');
fprintf('   ------------------------------------\n');
for kk=k1+1:k1+k2
fprintf('Beta%1.0f          %10.4f   %10.4f   [%10.4f, %10.4f]\n',[kk b_mean(kk) sqrt(b_var(kk)) HPDI_beta(kk,:)]);
end
fprintf('   ------------------------------------\n');
fprintf('Korelace nahodnych slozek mezi rovnicemi\n');
fprintf('   ------------------------------------\n');
fprintf('corr(e1,e2)    %10.4f   %10.4f   [%10.4f, %10.4f]\n',[r_mean sqrt(r_var) HPDI_r]);
fprintf('   ------------------------------------\n');
