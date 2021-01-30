clear all
close all
clc

%% Empiricka ilustrace - kapitola 4 - 4.2.7 - Koop 2003

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
Xstar = [1 5000 2 2 1];

[n,k] = size(X);

%% Apriorni hyperparametry
b0 = [0 10 5000 10000 10000]';

s02 = 5000^2;
h0 = 1/s02;
varb0 = [10000^2 25 2500^2 5000^2 5000^2]';
nu0 = 5;
V0 = diag(varb0);

S = 10000;
S0 = 1000;
S1 = S-S0;



b_Gibbs = zeros(k,S);
h_Gibbs = zeros(1,S);
SD_nom = zeros(k,S); %citatel pro BF pomoci Savage-Dickey
ystar = zeros(1,S);

%% Gibbsuv vzorkovac
%prvni beh - je mozno dat vse do jednoho cyklu, kdy budeme mit S+1
%replikaci, kdy prvni vektor parametru odpovida nasim pocatecnim podminkam
h_Gibbs0 = h0;
V1 = inv(inv(V0)+h_Gibbs0*X'*X);
b1 = V1*(inv(V0)*b0+h_Gibbs0*X'*y);

b_Gibbs(:,1)=b1+norm_rnd(V1);
%b_Gibbs(:,1)=mvnrnd(b1,V1)';

nu1 =n+nu0;
s12 = 1/nu1*((y-X*b_Gibbs(:,1))'*(y-X*b_Gibbs(:,1))+nu0*s02);
    %A=nu1/2; %jinak definovana funkce gamma v ramci gamrnd, nez v Koopovi
    %B=2*inv(s12)/nu1;
%h_Gibbs(:,1)=gamrnd(A,B);
h_Gibbs(:,1)=gamm_rnd(1,1,.5*nu1,.5*nu1*s12);


fprintf('Posterior simulation \n');
fprintf('0%%       50%%      100%% \n');
fprintf('|');
for i=2:S
   V1 = inv(inv(V0)+h_Gibbs(1,i-1)*X'*X);
   b1 = V1*(inv(V0)*b0+h_Gibbs(1,i-1)*X'*y);
    %b_Gibbs(:,i)=mvnrnd(b1,V1)';
   b_Gibbs(:,i)=b1+norm_rnd(V1); %spolehlivejsi generator!!!
    %nu1 =n+nu0; %v ramci cyklu se jiz nemeni
   s12 = 1/nu1*((y-X*b_Gibbs(:,i))'*(y-X*b_Gibbs(:,i))+nu0*s02);
    %A=nu1/2; %jinak definovana funkce gamma v ramci gamrnd, nez v Koopovi
    %B=2*inv(s12)/nu1;
    %h_Gibbs(:,i)=gamrnd(A,B);
   h_Gibbs(:,i)=gamm_rnd(1,1,.5*nu1,.5*nu1*s12);
    
    %citatel pro Savage-Dickey ratio M1: beta_j=0
    for j = 1:k
           SD_nom(j,i) = norm_pdf(0,b1(j,1),V1(j,j)); 
    end
    
    %predikce
    %podmineny vyber z predikcni hustoty, podmineno betou a h
        ystar(:,i) = Xstar*b_Gibbs(:,i) + norm_rnd(1/h_Gibbs(:,i));
    
        %graficky detail pro zobrazeni prubehu simulace
    if mod(i/S*100,5) == 0;    
        fprintf('|');
    end
end
fprintf('\n');
fprintf('Done!\n');

%% vyhozeni prvnich S0 vzorku
b_S1 = b_Gibbs(:,S0+1:S);
h_S1 = h_Gibbs(:,S0+1:S);
SD_nom(:,1:S0)=[];
ystar(:,1:S0) =[]; 

b_mean = mean(b_S1,2); %posteriorni stredni hodnota 
h_mean = mean(h_S1);
ystar_mean = mean(ystar);
b_var = mean(b_S1.^2,2)-b_mean.^2; %posteriorni rozptyl
h_var = mean(h_S1.^2)-h_mean.^2;
ystar_var = mean(ystar.^2)-ystar_mean.^2;

%Savage-Dickey ratio (Bayes factor)
%jmenovatel pro Savage-Dickey ratio M1: beta_j=0
SD_denom = zeros(k,1);
for j = 1:k
    SD_denom(j,1) = norm_pdf(0,b0(j,1),V0(j,j));
end
BF = mean(SD_nom,2)./SD_denom; 


%% KONVERGENCNI DIAGNOSTIKY - Gewekova CD pomoci funkce momentg
alldraws = [b_S1' h_S1'];
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



%% KONVERGENCNI DIAGNOSTIKY - NSE a Gewekova CD pomoci funkce psd (power
%spectral density)
NSE_psd = zeros(k,1);
CD_psd = zeros(k,1);
for i=1:k
    x = b_S1(i,:);
    sx = pwelch(x);
    %sx = psd(x);
    NSE_psd(i,1) = sqrt(sx(1)/length(x)^2);
    xa = b_S1(i,1:idraw1);
    xc = b_S1(i,idraw2:S1);
    sxa = pwelch(xa);
    sxc = pwelch(xc);
    %sxa = psd(xa);
    %sxc = psd(xc);
    CD_psd(i,1) = (mean(xa)-mean(xc))/(sqrt(sxa(1)/length(xa)^2)+sqrt(sxc(1)/length(xc)^2));
end

t = toc;

%% Prezentace vysledku, grafy
fprintf('Prior and posterior means for beta (std. deviations in parentheses)\n');
fprintf('                Prior       Posterior      NSE          NSE      Geweke CD     Geweke CD     Post. Odds\n');
fprintf('                                                       (psd)                     (psd)       (beta_j=0)\n');
fprintf('   ------------------------------------\n');

figure

for i=1:k
fprintf('Beta%1.0f     %12.4f   %12.4f   %10.4f  %10.4f   %10.4f   %10.4f   %10.4f\n',[i b0(i) b_mean(i) nse(i) NSE_psd(i) CD(i) CD_psd(i) BF(i)]);
fprintf('           (%12.4f) (%12.4f)\n',[sqrt(varb0(i)) sqrt(b_var(i))]);

    [f,xi] = ksdensity(b_S1(i,:));
    subplot(3,2,i)
    plot(xi,f)
    xlabel(['\beta_' num2str(i)])

end
fprintf('   ------------------------------------\n');
fprintf('Number of replications: %12.4f total \n',S);
fprintf('                        %12.4f discarded \n',S0);
fprintf('   ------------------------------------\n');
fprintf('Time elapsed (s): %12.2f \n',t);
fprintf('***************************************************************************\n\n\n');


figure
hist(ystar,25)
xlabel('House price')






