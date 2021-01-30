function results = nlrm_ncp(y,X,b0,V0,s02,nu0,Xstar)
%Funkce pro normalni linearni regresni model s priorzene konjugovanou
%apriorni hustotou
%   y ... vektor vysvetlovane promenne (n x 1)
%   X ... matice vysvetlujicich promennych (n x k)
%   b0 ... vektor apriornich strednich hodnot parametru b [k x 1]
%   V0 ... apriorni kovariancni matice parametru b [k x k]
%   s02 ... apriorni stredni hodnota pro rozptyl nahodnych slozek [1 x 1]
%   nu0 ... apriorni pocet stupnu volnosti [1 x 1] -- druhy parametr
%   hustoty gamma dle Koop (2003)
%   Xstar ... matice vysveltujicich pro predikci mimo vzorek [T x k]


nIn = nargin;

if nIn < 6
  error('Not enough input arguments.');
elseif nIn > 7
  error('Too many input arguments.');
end

[n,k] = size(X);

covb0 = nu0*s02/(nu0-2)*V0;
bstd0 = sqrt(diag(covb0));
hstd0 = sqrt(2/(nu0*s02^2));
h0=1/s02;

%posteriorni analyza
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



bmean1 =  b1;
covb1 = V1*nu1s12/(nu1-2); 
bstd1 = sqrt(diag(covb1)); %smerodatne odchylky pro posteriorni beta

%Posteriorni stredni hodnota a rozptyl presnosti chyby h
hmean1 = 1/s12;
hvar1=2/(nu1*s12^2);
hstd1=sqrt(hvar1);

%Predikce
    results.pred = 0;
   if nIn == 7
    ystarmean = Xstar*b1;
    ystarV = (eye(size(Xstar,1))+Xstar*V1*Xstar')*s12;
    covystar = ystarV*nu1/(nu1-2);
    ystarstd=sqrt(diag(covystar));
    results.Xstar=Xstar;
    results.ystarmean = ystarmean;
    results.covystar = covystar;
    results.ystarstd=ystarstd;
    results.pred = 1;
   end

%logaritmus marginalni verohodnosti pro model s informativni apriorni
%hustotou
    intcon=gammaln(nu1/2) + .5*nu0*log(nu0*s02)- gammaln(.5*nu0) -.5*n*log(pi); %vztah pro integracni konstantu
    lmarglik=intcon + .5*log(det(V1)/det(V0)) - .5*nu1*log(nu1s12); %logaritmus marginalni verohodnosti
   
   
%RESULTS
%apriorni hyperparametry
results.y=y;
results.X=X;
results.b0=b0;
results.V0=V0;
results.s02=s02;
results.nu0=nu0;
results.h0 = h0;

results.covb0 = covb0;
results.bstd0 = bstd0;
results.hstd0 = hstd0;

%posteriorni hyperparametry
results.b1=b1;
results.V1=V1;
results.s12=s12;
results.nu1=nu1;

results.bmean1 = bmean1;

results.covb1 = covb1;
results.bstd1 = bstd1;
results.hmean1 = hmean1;
results.hstd1 = hstd1;

results.lmarglik = lmarglik;





