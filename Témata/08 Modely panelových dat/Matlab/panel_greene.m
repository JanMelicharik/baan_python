clear
close all
clc

%% Nastaveni cest a dat
addpath('.\Support'); 

%Table F7.1 - Greene (2002) - Econometric analysis - 5th ed. 
%I   T        C             Q         PF            LF
data = [...
1    1     1140640       .952757    106650       .534487
1    2     1215690       .986757    110307       .532328
1    3     1309570      1.091980    110574       .547736
1    4     1511530      1.175780    121974       .540846
1    5     1676730      1.160170    196606       .591167
1    6     1823740      1.173760    265609       .575417
1    7     2022890      1.290510    263451       .594495
1    8     2314760      1.390670    316411       .597409
1    9     2639160      1.612730    384110       .638522
1   10     3247620      1.825440    569251       .676287
1   11     3787750      1.546040    871636       .605735
1   12     3867750      1.527900    997239       .614360
1   13     3996020      1.660200    938002       .633366
1   14     4282880      1.822310    859572       .650117
1   15     4748320      1.936460    823411       .625603
2    1      569292       .520635    103795       .490851
2    2      640614       .534627    111477       .473449
2    3      777655       .655192    118664       .503013
2    4      999294       .791575    114797       .512501
2    5     1203970       .842945    215322       .566782
2    6     1358100       .852892    281704       .558133
2    7     1501350       .922843    304818       .558799
2    8     1709270      1.000000    348609       .572070
2    9     2025400      1.198450    374579       .624763
2   10     2548370      1.340670    544109       .628706
2   11     3137740      1.326240    853356       .589150
2   12     3557700      1.248520   1003200       .532612
2   13     3717740      1.254320    941977       .526652
2   14     3962370      1.371770    856533       .540163
2   15     4209390      1.389740    821361       .528775
3    1      286298       .262424    118788       .524334
3    2      309290       .266433    123798       .537185
3    3      342056       .306043    122882       .582119
3    4      374595       .325586    131274       .579489
3    5      450037       .345706    222037       .606592
3    6      510412       .367517    278721       .607270
3    7      575347       .409937    306564       .582425
3    8      669331       .448023    356073       .573972
3    9      783799       .539595    378311       .654256
3   10      913883       .539382    555267       .631055
3   11     1041520       .467967    850322       .569240
3   12     1125800       .450544   1015610       .589682
3   13     1096070       .468793    954508       .587953
3   14     1198930       .494397    886999       .565388
3   15     1170470       .493317    844079       .577078
4    1      145167       .086393    114987       .432066
4    2      170192       .096740    120501       .439669
4    3      247506       .141500    121908       .488932
4    4      309391       .169715    127220       .484181
4    5      354338       .173805    209405       .529925
4    6      373941       .164272    263148       .532723
4    7      420915       .170906    316724       .549067
4    8      474017       .177840    363598       .557140
4    9      532590       .192248    389436       .611377
4   10      676771       .242469    547376       .645319
4   11      880438       .256505    850418       .611734
4   12     1052020       .249657   1011170       .580884
4   13     1193680       .273923    951934       .572047
4   14     1303390       .371131    881323       .594570
4   15     1436970       .421411    831374       .585525
5    1       91361       .051028    118222       .442875
5    2       95428       .052646    116223       .462473
5    3       98187       .056348    115853       .519118
5    4      115967       .066953    129372       .529331
5    5      138382       .070308    243266       .557797
5    6      156228       .073961    277930       .556181
5    7      183169       .084946    317273       .569327
5    8      210212       .095474    358794       .583465
5    9      274024       .119814    397667       .631818
5   10      356915       .150046    566672       .604723
5   11      432344       .144014    848393       .587921
5   12      524294       .169300   1005740       .616159
5   13      530924       .172761    958231       .605868
5   14      581447       .186670    872924       .594688
5   15      610257       .213279    844622       .635545
6    1       68978       .037682    117112       .448539
6    2       74904       .039784    119420       .475889
6    3       83829       .044331    116087       .500562
6    4       98148       .050245    122997       .500344
6    5      118449       .055046    194309       .528897
6    6      133161       .052462    307923       .495361
6    7      145062       .056977    323595       .510342
6    8      170711       .061490    363081       .518296
6    9      199775       .069027    386422       .546723
6   10      276797       .092749    564867       .554276
6   11      381478       .112640    874818       .517766
6   12      506969       .154154   1013170       .580049
6   13      633388       .186461    930477       .556024
6   14      804388       .246847    851676       .537791
6   15     1009500       .304013    819476       .525775];
%I   T        C             Q         PF            LF

%Cost Data for U.S. Airlines, 90 Oservations On 6 Firms For 15 Years, 1970-1984
%Source: These data are a subset of a larger data set provided to the author by Professor Moshe Kim.
%They were originally constructed by Christensen Associates of Madison, Wisconsin.

%    * I = Airline,
%    * T = Year,
%    * Q = Output, in revenue passenger miles, index number,
%    * C = Total cost, in $1000,
%    * PF = Fuel price,
%    * LF = Load factor, the average capacity utilization of the fleet. 

%{
a panel data study of a group of U.S. airlines. We will
fit a simple model for the total cost of production:

ln cost_it = b1 + b2 ln output_it + b3 ln fuel-price_it + b4 load-factor_it + eps_it .

Output is measured in “revenue passenger miles.” The load factor is a rate of capacity
utilization; it is the average rate at which seats on the airline’s planes are filled. More complete
models of costs include other factor prices (materials, capital) and, perhaps, a quadratic term
in log output to allow for variable economies of scale.
%}

%*********************************************

N = 6;  %pocet aerolinii
T = 15; %pocet let (1970-1984)
k = 4;  %pocet parametru

%% Tvorba datovych matic

y = log(data(:,3)); %logaritmus nakladu
y_i = reshape(y,T,N); %y_i pro i=1...N (po sloupcich)

X = [ones(N*T,1) log(data(:,4)) log(data(:,5)) data(:,6)];
X_i = zeros(T,k,N); %X_i pro i=1...N (po treti dimenzi)
for s = 1:N
   X_i(:,:,s) = X((T*s-T)+1:T*s,:); 
end
X_i_til = X_i(:,2:end,:); %X_i s vlnovkou (bez urovnove konstanty)

X_ast = zeros(T*N,(N+k-1)); 
for s = 1:N
   X_ast((T*s-T)+1:T*s,s)=ones(T,1);
   X_ast((T*s-T)+1:T*s,end-(k-1)+1:end)=X_i_til(:,:,s);
end

tic

%% Pooled model
%"standardni" Gibbsuv vzorkovac - jen vetsi matice a vektory
%prior hyperparameters
    beta_0 = [1 1 1 -1]';

    h_0 = 1/0.5^2;
    var_beta_0 = [1^2 1^2 1^2 1^2]';
    nu_0 = 1;
    var_h_0 = 2*h_0^2/nu_0;
    V_0 = diag(var_beta_0);

    S0 = 1000+1;
    S1 = 10000;
    S = S0+S1;

beta = zeros(k,S);
h = zeros(1,S);

%% Gibbsuv vzorkovac
h(1) = h_0;

 fprintf('Posteriorni simulace\n')
 fprintf('0%%     50%%      100%%\n')
 fprintf('|')
 %pomocna promenna pro graficke znazorneni prubehu
 pom_step = 0.05; %0.05 = posun po 5% 
 pom_graph = pom_step;  
 
nu_1 =T*N+nu_0;

for s=2:S

   V_1 = inv(inv(V_0)+h(1,s-1)*X'*X);
   beta_1 = V_1*(inv(V_0)*beta_0+h(1,s-1)*X'*y);
   beta(:,s)=beta_1+norm_rnd(V_1);  %podminena hustota pro beta
   h_1 = (1/nu_1*((y-X*beta(:,s))'*(y-X*beta(:,s))+nu_0*1/h_0))^-1;
   h(:,s)=gamm_rnd_Koop(h_1,nu_1,1); %podminena hustota pro h
 
    %graficke zobrazeni prubehu po 5 %
    if s/S>=pom_graph    
     fprintf('|');
     pom_graph = pom_graph+pom_step;
    end
end

fprintf('\n');
fprintf('Hotovo!\n');


%vyhozeni prvnich S0 vzorku
    beta(:,1:S0) = [];
    h(:,1:S0) = [];

    beta_mean = mean(beta,2); %posteriorni stredni hodnota 
    h_mean = mean(h);
    var_beta = mean(beta.^2,2)-beta_mean.^2; %posteriorni rozptyl
    var_h = mean(h.^2)-h_mean.^2;

%KONVERGENCNI DIAGNOSTIKY - pro uspornost neresime
tt = toc;

%Prezentace vysledku, grafy
fprintf('Pooled model\n');
fprintf('Prior and posterior means for beta (std. deviations in parentheses)\n');
fprintf('                Prior       Posterior      \n');
fprintf('                                           \n');
fprintf('   ------------------------------------\n');

for s=1:k
    fprintf('Beta%1.0f     %12.4f   %12.4f   \n',[s beta_0(s) beta_mean(s)]);
    fprintf('           (%12.4f) (%12.4f)\n',[sqrt(var_beta_0(s)) sqrt(var_beta(s))]);
end

    fprintf('h         %12.4f   %12.4f   \n',[h_0 h_mean]);
    fprintf('           (%12.4f) (%12.4f)\n',[sqrt(var_h_0) sqrt(var_h)]);


fprintf('   ------------------------------------\n');
fprintf('Number of replications: %12.4f total \n',S-1);
fprintf('                        %12.4f discarded \n',S0-1);
fprintf('   ------------------------------------\n');
fprintf('Time elapsed (s): %12.2f \n',tt);
fprintf('***************************************************************************\n\n\n');











%% Individual effects model - nehierarchicka apriorni hustota
%prior hyperparameters
    alpha_0 = ones(N,1)';       %variabilni urovnove konstanty napric aerolinkami
    beta_0 = [alpha_0 1 1 -1]';     %stejne ostatni parametry (sklonu)

    h_0 = 1/0.5^2;
    
    var_alpha_0 = ones(N,1)';
    var_beta_0 = [var_alpha_0 1^2 1^2 1^2]';
    nu_0 = 1;
    var_h_0 = 2*h_0^2/nu_0;
    V_0 = diag(var_beta_0);

    S0 = 1000+1;
    S1 = 10000;
    S = S0+S1;

beta = zeros(N+k-1,S);
h = zeros(1,S);

%% Gibbsuv vzorkovac
h(1) = h_0;
nu_1 =T*N+nu_0;

 fprintf('Posteriorni simulace\n')
 fprintf('0%%     50%%      100%%\n')
 fprintf('|')
 %pomocna promenna pro graficke znazorneni prubehu
 pom_step = 0.05; %0.05 = posun po 5% 
 pom_graph = pom_step;  
 
for s=2:S
   V_1 = inv(inv(V_0)+h(1,s-1)*X_ast'*X_ast);
   beta_1 = V_1*(inv(V_0)*beta_0+h(1,s-1)*X_ast'*y);
   beta(:,s)=beta_1+norm_rnd(V_1);  %podminena hustota pro beta
   pom = 0;
   for j=1:N
      pom = pom+(y_i(:,j)-beta(j,s)*ones(T,1)-X_i_til(:,:,j)*beta(N+1:N+(k-1),s))'*(y_i(:,j)-beta(j,s)*ones(T,1)-X_i_til(:,:,j)*beta(N+1:N+(k-1),s));
   end
   h_1 = (1/nu_1*(pom+nu_0*1/h_0))^-1;
   h(:,s)=gamm_rnd_Koop(h_1,nu_1,1); %podminena hustota pro h
   
    %graficke zobrazeni prubehu po 5 %
    if s/S>=pom_graph    
     fprintf('|');
     pom_graph = pom_graph+pom_step;
    end
end

fprintf('\n');
fprintf('Hotovo!\n');

%% vyhozeni prvnich S0 vzorku
    beta(:,1:S0) = [];
    h(:,1:S0) = [];

    beta_mean = mean(beta,2); %posteriorni stredni hodnota 
    h_mean = mean(h);
    var_beta = mean(beta.^2,2)-beta_mean.^2; %posteriorni rozptyl
    var_h = mean(h.^2)-h_mean.^2;

%KONVERGENCNI DIAGNOSTIKY - pro uspornost neresim 
tt = toc;

%Prezentace vysledku, grafy
fprintf('Individual Effects Model - nehierarchicka apriorni hustota\n');
fprintf('Prior and posterior means for beta (std. deviations in parentheses)\n');
fprintf('                Prior       Posterior      \n');
fprintf('                                           \n');
fprintf('   ------------------------------------\n');

for s=1:N
    fprintf('Alfa%1.0f     %12.4f   %12.4f   \n',[s beta_0(s) beta_mean(s)]);
    fprintf('           (%12.4f) (%12.4f)\n',[sqrt(var_beta_0(s)) sqrt(var_beta(s))]);
end
for s=N+1:N+(k-1)
    fprintf('Beta%1.0f     %12.4f   %12.4f   \n',[s-N+1 beta_0(s) beta_mean(s)]);
    fprintf('           (%12.4f) (%12.4f)\n',[sqrt(var_beta_0(s)) sqrt(var_beta(s))]);
end

    fprintf('h         %12.4f   %12.4f   \n',[h_0 h_mean]);
    fprintf('           (%12.4f) (%12.4f)\n',[sqrt(var_h_0) sqrt(var_h)]);


fprintf('   ------------------------------------\n');
fprintf('Number of replications: %12.4f total \n',S-1);
fprintf('                        %12.4f discarded \n',S0-1);
fprintf('   ------------------------------------\n');
fprintf('Time elapsed (s): %12.2f \n',tt);
fprintf('***************************************************************************\n\n\n');







%% Individual effects model - hierarchicka apriorni hustota
%prior hyperparameters
    alpha_0 = ones(N,1)';       %variabilni urovnove konstanty napric aerolinkami
    beta_0 = [1 1 -1]';     %stejne ostatni parametry (sklonu)

    h_0 = 1/0.5^2;
    
    var_alpha_0 = ones(N,1)';
    var_beta_0 = [1^2 1^2 1^2]';
    nu_0 = 1;
    var_h_0 = 2*h_0^2/nu_0;
    V_0 = diag(var_beta_0);

    mu_alpha_0 = 5;
    h_alpha_0 = 1/(2^2);
    V_alpha_0 = 5;
    nu_alpha_0 = 1;
    
    S0 = 1000+1;
    S1 = 10000;
    S = S0+S1;

beta = zeros(k-1,S);
h = zeros(1,S);
alpha = zeros(N,S);

V_alpha = zeros(1,S);
mu_alpha = zeros(1,S);

%% Gibbsuv vzorkovac
h(1) = h_0;
V_alpha(1) = V_alpha_0;
mu_alpha(1) = mu_alpha_0;
alpha(:,1) = alpha_0;

nu_1 =T*N+nu_0;
nu_alpha_1 =nu_alpha_0+N;
 
 fprintf('Posteriorni simulace\n')
 fprintf('0%%     50%%      100%%\n')
 fprintf('|')
 %pomocna promenna pro graficke znazorneni prubehu
 pom_step = 0.05; %0.05 = posun po 5% 
 pom_graph = pom_step;  
 
for s=2:S

pom = 0;
   for j=1:N
      pom = pom+X_i_til(:,:,j)'*X_i_til(:,:,j);
   end
V_1 = inv(inv(V_0)+h(:,s-1)*pom);

pom = 0;
   for j=1:N
      pom = pom+X_i_til(:,:,j)'*(y_i(:,j)-alpha(j,s-1)*ones(T,1));
   end
beta_1 = V_1*(inv(V_0)*beta_0+h(:,s-1)*pom);
beta(:,s)=beta_1+norm_rnd(V_1); %podminena hustota pro beta

   pom = 0;
   for j=1:N
      pom = pom+(y_i(:,j)-alpha(j,s-1)*ones(T,1)-X_i_til(:,:,j)*beta(:,s))'*(y_i(:,j)-alpha(j,s-1)*ones(T,1)-X_i_til(:,:,j)*beta(:,s));
   end
   h_1 = (1/nu_1*(pom+nu_0*1/h_0))^(-1);
h(:,s)=gamm_rnd_Koop(h_1,nu_1,1); %podminena hustota pro h

%podminene hustoty pro kazde alpha_j
for j=1:N
   V_1_j = (V_alpha(:,s-1)*h(:,s)^-1)/(T*V_alpha(:,s-1)+h(:,s)^-1);
   alpha_1 = (V_alpha(:,s-1)*(y_i(:,j)-X_i_til(:,:,j)*beta(:,s))'*ones(T,1)+h(:,s)^-1*mu_alpha(:,s-1))/(T*V_alpha(:,s-1)+h(:,s)^-1);
   alpha(j,s) = alpha_1+norm_rnd(V_1_j);  %Koop 7.16 
end

%podminena hustota pro mu_alfa
sigma_2_alpha_1 = V_alpha(:,s-1)*1/h_alpha_0/(V_alpha(:,s-1)+N*1/h_alpha_0);
mu_alpha_1 = (V_alpha(:,s-1)*mu_alpha_0+1/h_alpha_0*sum(alpha(:,s)))/(V_alpha(:,s-1)+N*1/h_alpha_0);
mu_alpha(:,s) = mu_alpha_1+norm_rnd(sigma_2_alpha_1); %Koop 7.17

%podminena hustota pro inv(V_alfa)
   pom = 0;
   for j=1:N
      pom = pom+(alpha(j,s)-mu_alpha(:,s))^2;
   end
V_alpha_1 = (pom+V_alpha_0*nu_alpha_0)/nu_alpha_1;
V_alpha(:,s)=inv(gamm_rnd_Koop(inv(V_alpha_1),nu_alpha_1,1)); %Koop 7.18


    %graficke zobrazeni prubehu po 5 %
    if s/S>=pom_graph    
     fprintf('|');
     pom_graph = pom_graph+pom_step;
    end
end

fprintf('\n');
fprintf('Hotovo!\n');



%vyhozeni prvnich S0 vzorku
    alpha(:,1:S0) = [];
    beta(:,1:S0) = [];
    h(1:S0) = [];

    beta_mean = mean(beta,2); %posteriorni stredni hodnota 
    alpha_mean = mean(alpha,2);
    h_mean = mean(h);
    var_beta = mean(beta.^2,2)-beta_mean.^2; %posteriorni rozptyl
    var_alpha = mean(alpha.^2,2)-alpha_mean.^2;
    var_h = mean(h.^2)-h_mean.^2;

%KONVERGENCNI DIAGNOSTIKY - pro uspornost neresim 
tt = toc;

%Prezentace vysledku, grafy
fprintf('Individual Effects Model - hierarchicka apriorni hustota\n');
fprintf('Prior and posterior means for beta (std. deviations in parentheses)\n');
fprintf('                Prior       Posterior      \n');
fprintf('                                           \n');
fprintf('   ------------------------------------\n');

for s=1:N
    fprintf('Alfa%1.0f     %12.4f   %12.4f   \n',[s alpha_0(s) alpha_mean(s)]);
    fprintf('           (%12.4f) (%12.4f)\n',[sqrt(var_alpha_0(s)) sqrt(var_alpha(s))]);
end
for s=1:(k-1)
    fprintf('Beta%1.0f     %12.4f   %12.4f   \n',[s+1 beta_0(s) beta_mean(s)]);
    fprintf('           (%12.4f) (%12.4f)\n',[sqrt(var_beta_0(s)) sqrt(var_beta(s))]);
end

    fprintf('h         %12.4f   %12.4f   \n',[h_0 h_mean]);
    fprintf('           (%12.4f) (%12.4f)\n',[sqrt(var_h_0) sqrt(var_h)]);


fprintf('   ------------------------------------\n');
fprintf('Number of replications: %12.4f total \n',S-1);
fprintf('                        %12.4f discarded \n',S0-1);
fprintf('   ------------------------------------\n');
fprintf('Time elapsed (s): %12.2f \n',tt);
fprintf('***************************************************************************\n\n\n');








%% Random coefficients model - hierarchicka apriorni hustota
%prior hyperparameters
 beta_0 = zeros(k,N);
 for j=1:N
    beta_0(:,j) = [1 1 1 -1]';     %vsechny parametry variabilni
 end
    h_0 = 1/0.5^2;
    
    var_beta_0 = [1^2 1^2 1^2 1^2]';
    nu_0 = 1;
    var_h_0 = 2*h_0^2/nu_0;
    V_0 = diag(var_beta_0);

    mu_beta_0 = [1 1 1 -1]';
    Sigma_beta_0 = diag(var_beta_0);
    V_beta_0 = V_0;
    nu_beta_0 = 1;
    
    
    
    S0 = 1000+1;
    S1 = 10000;
    S = S1+S0;

beta = zeros(k,S,N);
h = zeros(1,S);

V_beta = zeros(k,k,S);
mu_beta = zeros(k,S);

%% Gibbsùv vzorkovaè
h(1) = h_0;
V_beta(:,:,1) = V_beta_0;
mu_beta(:,1) = mu_beta_0;

nu_beta_1 = N+nu_beta_0;
nu_1 = T*N+nu_0;

 fprintf('Posteriorni simulace\n')
 fprintf('0%%     50%%      100%%\n')
 fprintf('|')
 %pomocna promenna pro graficke znazorneni prubehu
 pom_step = 0.05; %0.05 = posun po 5% 
 pom_graph = pom_step;  

for i=2:S
 for j=1:N
    V_1 = inv(h(:,i-1)*X_i(:,:,j)'*X_i(:,:,j)+inv(V_beta(:,:,i-1)));
    beta_1 = V_1*(h(:,i-1)*X_i(:,:,j)'*y_i(:,j)+inv(V_beta(:,:,i-1))*mu_beta(:,i-1));
    beta(:,i,j)=beta_1+norm_rnd(V_1); %podminena hustota pro beta - Koop 7.25
 end
 
%podminena hustota pro mu_beta - Koop (7.26) a naseldujici
Sigma_beta_1 = inv(N*inv(V_beta(:,:,i-1))+inv(Sigma_beta_0));
   pom = 0;
   for j=1:N
      pom = pom+beta(:,i,j);
   end
mu_beta_1 = Sigma_beta_1*(inv(V_beta(:,:,i-1))*pom+inv(Sigma_beta_0)*mu_beta_0);
mu_beta(:,i) = mu_beta_1+norm_rnd(Sigma_beta_1);

   pom = 0;
   for j=1:N
      pom = pom+(beta(:,i,j)-mu_beta(:,i))*(beta(:,i,j)-mu_beta(:,i))';
   end
V_beta_1 = pom+V_beta_0;
V_beta(:,:,i) = inv(wish_rnd(inv(nu_beta_1*V_beta_1),nu_beta_1));

%podminena hustota pro h
   pom = 0;
   for j=1:N
      pom = pom+(y_i(:,j)-X_i(:,:,j)*beta(:,i,j))'*(y_i(:,j)-X_i(:,:,j)*beta(:,i,j));
   end
   h_1 = (1/nu_1*(pom+nu_0*1/h_0))^-1;
   h(:,i)=gamm_rnd_Koop(h_1,nu_1,1); %podminena hustota pro h

     %graficke zobrazeni prubehu po 5 %
    if i/S>=pom_graph    
     fprintf('|');
     pom_graph = pom_graph+pom_step;
    end
end

fprintf('\n');
fprintf('Hotovo!\n');



%vyhozeni prvnich S0 vzorku
    beta(:,1:S0,:) = [];
    h(:,1:S0) = [];
    mu_beta(:,1:S0) = [];

    beta_mean = zeros(k,N);
    for j = 1:N
    beta_mean(:,j) = mean(beta(:,:,j),2); %posteriorni stredni hodnota 
    end
    h_mean = mean(h);
    mu_beta_mean = mean(mu_beta,2);

    var_beta = zeros(k,N);
    for j = 1:N
    var_beta(:,j) = mean(beta(:,:,j).^2,2)-beta_mean(:,j).^2; %posteriorni rozptyl 
    end
    var_h = mean(h.^2)-h_mean.^2;
    var_mu_beta = mean(mu_beta.^2,2)-mu_beta_mean.^2;

    
%KONVERGENCNI DIAGNOSTIKY - pro uspornost neresime
tt = toc;

%Prezentace vysledku, grafy
fprintf('Random Coefficients Model - hierarchicka apriorni hustota\n');
fprintf('Posterior means for beta (std. deviations in parentheses)\n');
fprintf('   Firma     1            2            3            4            5            6\n');
fprintf('   -----------------------------------------------------------------------------------\n');

for i=1:k
    fprintf('beta%1.0f   %10.4f   %10.4f   %10.4f   %10.4f   %10.4f   %10.4f \n',[i beta_mean(i,:)]);
    fprintf('           (%10.4f) (%10.4f) (%10.4f) (%10.4f) (%10.4f) (%10.4f)\n',[sqrt(var_beta(i,:))]);
end

    fprintf('h       %10.4f   %10.4f   \n',[h_0 h_mean]);
    fprintf('           (%10.4f) (%10.4f)\n',[sqrt(var_h_0) sqrt(var_h)]);


fprintf('   ------------------------------------\n');
fprintf('Number of replications: %12.4f total \n',S-1);
fprintf('                        %12.4f discarded \n',S0-1);
fprintf('   ------------------------------------\n');
fprintf('Time elapsed (s): %12.2f \n',tt);
fprintf('***************************************************************************\n\n\n');




