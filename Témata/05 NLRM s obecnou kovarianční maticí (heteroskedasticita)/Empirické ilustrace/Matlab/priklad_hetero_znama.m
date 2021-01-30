clear all
close all
clc

%% Empiricka ilustrace - kapitola 6 - 6.3.2 - Koop 2003
tic %meric casu - pro zjisteni rychlosti programu


%% Nacteni dat
load hprice.txt

price = hprice(:,1); %prodejni cena domu
lot_size = hprice(:,2); %rozloha ve stopach ctverecnich
n_bed = hprice(:,3); %poèet ložnic
n_bath = hprice(:,4); %poèet koupelen
n_storey = hprice(:,5); %poèet pater

y = price;
X = [ones(size(y)) lot_size n_bed n_bath n_storey];

[n,k] = size(X);

%% Apriorni hyperparametry
b0 = [0 10 5000 10000 10000]';

s02 = 5000^2;
h0 = 1/s02;
varb0 = [10000^2 25 2500^2 5000^2 5000^2]';
nu0 = 5;
V0 = diag(varb0);
varh0 = 2*h0^2/nu0;

S = 50000;
S0 = 20000;
S1 = S-S0;


%promenne "zpusobujici" heteroskedasticitu
z = [lot_size n_bed n_bath n_storey];
p = size(z,2);


b_Gibbs = zeros(k,S);
h_Gibbs = zeros(1,S);
a_rw = zeros(p,S);






count_rw = 0; %èítaè akcepotvaných kandidátù

%po prvnich 5000 replikacich se vzal vysledny vektor rozptylu (a kovarianci) jako zaklad
%kovariancni matice; po te se dale ladi "d" pro ziskani zadouci
%akceptacni pravdepodobnosti;
%nebo - najde se modus posteriorni podminene hustoty pro alfa a prislusny
%hessian (podmineno apriornimi strednimi hodnotami, nebo jeste lepe
%posteriornimi strednimi hodnotami z prvnich 5000 replikaci

%postvar = eye(p);
%postvar = diag([0.132; 15.6423; 4.6199]); 
postvar = [     0.000000020294351  -0.000018452288594   0.000000678682051  -0.000004937071691
  -0.000018452288594   0.139027551523147   0.109426915494931  -0.073456296455522
   0.000000678682051   0.109426915494931   0.232998062000183  -0.110725086293211
  -0.000004937071691  -0.073456296455522  -0.110725086293211   0.098634611557309
];
c=0.3;
vscale_rw=c*postvar;
%a0 = [0.5;0.5 ;0.5; 0.5];
a0 = [   0.000546458183385
   0.580309604041334
   1.840716965714496
  -0.317194417927449];
%a0 = [-0.1707; -6.2956; 3.1309];
%a0= [1.4918; -7.4384; 4.3193];
Omega = diag((ones(n,1)+z*a0).^2); %kov. matice Omega
invOmega = inv(Omega);




fprintf('Posteriorni simulace\n');
fprintf('0%%       50%%      100%% \n');
fprintf('|');

%% Gibbsuv vzorkovac
% Obecny postup, ktery by jeste slo vypocetne zefektivnit s ohledem na to,
% ze Omega je diagonalni matice
b_Omega =inv(X'*invOmega*X)*X'*invOmega*y;
h_Gibbs0 = h0;
V1 = inv(inv(V0)+h_Gibbs0*X'*invOmega*X);
b1 = V1*(inv(V0)*b0+h_Gibbs0*X'*invOmega*X*b_Omega);
b_Gibbs(:,1)=b1+norm_rnd(V1);
%b_Gibbs(:,1)=mvnrnd(b1,V1)';

nu1 =n+nu0;
s12 = 1/nu1*((y-X*b_Gibbs(:,1))'*invOmega*(y-X*b_Gibbs(:,1))+nu0*s02);
    %A=nu1/2; %jinak definovana funkce gamma v ramci gamrnd, nez v Koopovi
    %B=2*inv(s12)/nu1;
%h_Gibbs(:,1)=gamrnd(A,B);
h_Gibbs(:,1)=gamm_rnd(1,1,.5*nu1,.5*nu1*s12);

%R-W M-H pro alpha
    a_can_rw = a0 + norm_rnd(vscale_rw); %kandidat
    log_accept_rw = min(a_post(a_can_rw,b_Gibbs(:,1),h_Gibbs(:,1),y,X,z)-a_post(a0,b_Gibbs(:,1),h_Gibbs(:,1),y,X,z),0);
    if log_accept_rw > log(rand)
        a_rw(:,1)=a_can_rw;
        count_rw=count_rw+1;
    else
        a_rw(:,1)=a0;
    end



for i=2:S
 Omega = diag((ones(n,1)+z*a_rw(:,i-1)).^2); %kov. matice Omega
 invOmega = inv(Omega);
 b_Omega =inv(X'*invOmega*X)*X'*invOmega*y;
 V1 = inv(inv(V0)+h_Gibbs(:,i-1)*X'*invOmega*X);
 b1 = V1*(inv(V0)*b0+h_Gibbs(:,i-1)*X'*invOmega*X*b_Omega);
 %(nesystematicky zasah -- nekdy nastava numericky problem, ze Matlab nevezme V1 jako pozitivne
 %definitni matici => ignorujeme chybu a "udelame z V1 pozitivne definitni
 %prictenim hodne maleho cisla 
 try
 b_Gibbs(:,i)=b1+norm_rnd(V1);
 catch
 b_Gibbs(:,i)=b1+norm_rnd(V1+0.000001);
 end
 %b_Gibbs(:,1)=mvnrnd(b1,V1)';

nu1 =n+nu0;
s12 = 1/nu1*((y-X*b_Gibbs(:,i))'*invOmega*(y-X*b_Gibbs(:,i))+nu0*s02);
    %A=nu1/2; %jinak definovana funkce gamma v ramci gamrnd, nez v Koopovi
    %B=2*inv(s12)/nu1;
%h_Gibbs(:,1)=gamrnd(A,B);
h_Gibbs(:,i)=gamm_rnd(1,1,.5*nu1,.5*nu1*s12);

%R-W M-H pro alpha
    a_can_rw = a_rw(:,i-1) + norm_rnd(vscale_rw); %kandidat
    log_accept_rw = min(a_post(a_can_rw,b_Gibbs(:,i),h_Gibbs(:,i),y,X,z)-a_post(a_rw(:,i-1),b_Gibbs(:,i),h_Gibbs(:,i),y,X,z),0);
    if log_accept_rw > log(rand)
        a_rw(:,i)=a_can_rw;
        count_rw=count_rw+1;
    else
        a_rw(:,i)=a_rw(:,i-1);
    end

        if mod(i/S*100,5) == 0;    
        fprintf('|');
        end
end
fprintf('\n');
fprintf('Done!\n');


%% Posteriorni analyza
%vyhozeni prvnich S0 vzorku
b_S1 = b_Gibbs(:,S0+1:S);
h_S1 = h_Gibbs(:,S0+1:S);
a_S1 = a_rw(:,S0+1:S);


b_mean = mean(b_S1,2); %posteriorni stredni hodnota 
h_mean = mean(h_S1);
a_mean = mean(a_S1,2);
b_var = mean(b_S1.^2,2)-b_mean.^2; %posteriorni rozptyl
h_var = mean(h_S1.^2)-h_mean.^2;
a_var = mean(a_S1.^2,2)-a_mean.^2;


accept_ratio_rw = count_rw/S;


%KONVERGENCNI DIAGNOSTIKY - Gewekova CD pomoci funkce momentg
alldraws = [b_S1' h_S1' a_S1'];
%The function momentg is taken from LeSage's toolbox
%it inputs all Gibbs draws and produces posterior
%mean, standard deviation, nse and rne
%it calculates what the book calls S(0) in various ways
%see momentg.m for more details
result = momentg(alldraws);
means=[result.pmean]';
stdevs=[result.pstd]';
nse=[result.nse]';
nse1=[result.nse1]';
nse2=[result.nse2]';
nse3=[result.nse3]';
%calculate Geweke convergence diagnostic based on first .1
%and last .4 of draws
idraw1= round(.1*S1); %prvnich 10%
result = momentg(alldraws(1:idraw1,:));
meansa=[result.pmean]';
nsea=[result.nse1]';

idraw2= round(.6*S1)+1; %poslednich 40%
result = momentg(alldraws(idraw2:S1,:));
meansb=[result.pmean]';
nseb=[result.nse1]';

CD = (meansa - meansb)./(nsea+nseb);



%HPDI pro beta, h a alpha
[k_b,l_b] = size(b_S1);
[k_a,l_a] = size(a_S1);

hperc=0.95; %HPDI
HPDIb_l = zeros(k_b,1);
HPDIb_u = zeros(k_b,1);
for j=1:k_b
    botperc=round((1-hperc)/2*l_b);
    topperc=round((1+hperc)/2*l_b);
    temp = sort(b_S1(j,:));
    HPDIb_u(j)=temp(topperc);
    HPDIb_l(j)=temp(botperc);
end

HPDIh_l = zeros(1,1);
HPDIh_u = zeros(1,1);
    botperc=round((1-hperc)/2*l_b);
    topperc=round((1+hperc)/2*l_b);
    temp = sort(h_S1);
    HPDIh_u=temp(topperc);
    HPDIh_l=temp(botperc);

HPDIa_l = zeros(k_a,1);
HPDIa_u = zeros(k_a,1);
for j=1:k_a
    botperc=round((1-hperc)/2*l_a);
    topperc=round((1+hperc)/2*l_a);
    temp = sort(a_S1(j,:));
    HPDIa_u(j)=temp(topperc);
    HPDIa_l(j)=temp(botperc);
end








t = toc;

%% Prezentace vysledku, grafy
fprintf('Prior and posterior means for beta (std. deviations in parentheses)\n');
fprintf('                Prior       Posterior      NSE         Geweke CD             %3.1f%% HPDI \n',hperc*100);
fprintf('   ------------------------------------\n');

figure

for i=1:k
fprintf('Beta%1.0f     %12.4f   %12.4f   %10.4f   %10.4f   [%10.4f, %10.4f]\n',[i b0(i) b_mean(i) nse(i) CD(i) HPDIb_l(i) HPDIb_u(i)]);
fprintf('           (%12.4f) (%12.4f)\n',[sqrt(varb0(i)) sqrt(b_var(i))]);

    [f,xi] = ksdensity(b_S1(i,:));
    subplot(3,2,i)
    plot(xi,f)
    xlabel(['\beta_' num2str(i)])
end

figure
fprintf('h          %12.4f   %12.4f   %10.4f   %10.4f   [%10.4f, %10.4f]\n',[h0 h_mean nse(k+1) CD(k+1) HPDIh_l HPDIh_u]);
fprintf('           (%12.4f)        (%12.4f)\n',[sqrt(varh0) sqrt(h_var)]);

    [f,xi] = ksdensity(h_S1);
    plot(xi,f)
    xlabel(['h'])



figure
for i=1:p
fprintf('Alfa%1.0f     %12.4f   %12.4f   %10.4f   %10.4f   [%10.4f, %10.4f]\n',[i a0(i) a_mean(i) nse(k+1+i) CD(k+1+i) HPDIa_l(i) HPDIa_u(i)]);
fprintf('           (N.A.)        (%12.4f)\n',sqrt(a_var(i)));

    [f,xi] = ksdensity(a_S1(i,:));
    subplot(3,2,i)
    plot(xi,f)
    xlabel(['\alpha_' num2str(i)])
end

fprintf('   ------------------------------------\n');
fprintf('Number of replications: %12.4f total \n',S);
fprintf('                        %12.4f discarded \n',S0);
fprintf('   ------------------------------------\n');

fprintf('   ------------------------------------\n');
fprintf('Accept.Ratio:  %8.4f\n\n', accept_ratio_rw)
fprintf('   ------------------------------------\n');

fprintf('Time elapsed (s): %12.2f \n',t);
fprintf('***************************************************************************\n\n\n');
