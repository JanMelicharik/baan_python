clear all
close all
clc

%% Empiricka ilustrace - kapitola 4 -4.3.6 - Koop 2003


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

%prior hyperparameters
b0 = [0 10 5000 10000 10000]';

s02 = 5000^2;
h0 = 1/s02;
varb0 = [10000^2 25 2500^2 5000^2 5000^2]';
nu0 = 5;
V0 = diag(varb0)*(nu0-2)/(nu0*s02);

S = 10000;
w = ones(S,1); %preddefinovani vektoru vah pro importance sampling


%% posteriorni analyza
%Odhady metodou nejmensich ctvercu
bols = inv(X'*X)*X'*y;
nu = n-k;
s2 = (y-X*bols)'*(y-X*bols)/nu;

%Posteriorni hyperparametry normalniho-gama rozdeleni
Xsquare=X'*X;
nu1=nu0+n;
V1 = inv(inv(V0)+Xsquare);
b1 = V1*(inv(V0)*b0 + Xsquare*bols);
nu1s12 = nu0*s02 + nu*s2 + (bols-b0)'*inv(V0+inv(Xsquare))*(bols-b0);
s12 = nu1s12/nu1;

beta = zeros(k,S);
ystar = zeros(size(Xstar,1),S);


Sigma1 = s12*V1;
Vchol = (chol(Sigma1))'; %choleskeho dekompozice matice Sigma
Vcholpred = (chol(s12*(ones(size(Xstar,1))+Xstar*V1*Xstar')))';
Sigma0 = s02*V0;


%% vypocet vah pro importance sampling
for i = 1:S
  
    %draw a t(0,1,v1) then transform to yield draw of beta
    beta(:,i)=b1 + Vchol*tdis_rnd(k,nu1);
   
        if beta(2,i)<5
            w(i,1)=0;
        end
        if beta(3,i)<2500
            w(i,1)=0;
        end
        if beta(4,i)<5000
            w(i,1)=0;
        end
        if beta(5,i)<5000
            w(i,1)=0;
        end
        
    %draw from predictive, conditional on beta and h
    ystar(:,i) = Xstar*beta(:,i) + Vcholpred*tdis_rnd(size(Xstar,1),nu1);

end


%% Vazene prumerovani pro odhad stredni hodnoty a rozptylu
bmean = (beta*w)./sum(w);
bsquare = ((beta.^2)*w)./sum(w);
bvar = bsquare - bmean.^2;
bsd=sqrt(bvar);
predmean = (ystar*w)./sum(w);
psquare = ((ystar.^2)*w)./sum(w);
predvar = psquare - predmean.^2;
predsd=sqrt(predvar);


%% calculate numerical standard errors
NSE=zeros(k,1);
for j = 1:k
   temp1 = (w.*beta(j,:)' - bmean(j,1)*w).^2;
   NSE(j,1) = mean(temp1)/(mean(w))^2;
   NSE(j,1)=sqrt(NSE(j,1)/S);
end


%% Savage-Dickey density ratio 
prior = zeros(k,1);
post = zeros(k,1);
for j = 1:k
    prior(j,1) = my_tpdf(b0(j,1),b0(j,1),Sigma0(j,j),nu0);
    post(j,1) = my_tpdf(b0(j,1),b1(j,1),Sigma1(j,j),nu1);
end
bf = post./prior; %Bayesuv faktor pro model beta_j=beta0_j a neomezeny model

%now obtain clower and cupper for each restriction
%and correct bayes factor appropriately
%korekce pro beta_2...5 vzhledem ktomu ze se jedna o omezene t-rozdeleni
zscore = (5-b0(2,1))/sqrt(Sigma0(2,2));
clower = 1/(1-tdis_cdf(zscore,nu0)); %integracni konstanta pro jmenovatel (apriorni hustota)
zscore = (5-b1(2,1))/sqrt(Sigma1(2,2));
cupper = 1/(1-tdis_cdf(zscore,nu1)); %integracni konstanta pro citatel (posteriorni hustota)
bf(2,1)= cupper*bf(2,1)/clower; 
zscore = (2500-b0(3,1))/sqrt(Sigma0(3,3));
clower = 1/(1-tdis_cdf(zscore,nu0));
zscore = (2500-b1(3,1))/sqrt(Sigma1(3,3));
cupper = 1/(1-tdis_cdf(zscore,nu1));
bf(3,1)= cupper*bf(3,1)/clower;
zscore = (5000-b0(4,1))/sqrt(Sigma0(4,4));
clower = 1/(1-tdis_cdf(zscore,nu0));
zscore = (5000-b1(4,1))/sqrt(Sigma1(4,4));
cupper = 1/(1-tdis_cdf(zscore,nu1));
bf(4,1)= cupper*bf(4,1)/clower;
zscore = (5000-b0(5,1))/sqrt(Sigma0(5,5));
clower = 1/(1-tdis_cdf(zscore,nu0));
zscore = (5000-b1(5,1))/sqrt(Sigma1(5,5));
cupper = 1/(1-tdis_cdf(zscore,nu1));
bf(5,1)= cupper*bf(5,1)/clower;



%% Vypis vysledku - beta
fprintf('Posterior results for beta\n');
fprintf('                  Mean        Std. Dev        NSE     BF beta_j=beta0_j\n');
for i=1:k
fprintf('Beta %2u      %12.4f %12.4f %12.4f %12.4f       \n',[i bmean(i) bsd(i) NSE(i) bf(i)]);
fprintf('   ------------------------------------\n');
end
%Predikce
fprintf('ystar         %12.4f %12.4f \n',[predmean predsd]);
