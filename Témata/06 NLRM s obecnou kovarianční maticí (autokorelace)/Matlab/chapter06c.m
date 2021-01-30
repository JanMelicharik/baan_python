clear
close all
clc

%% Nastaveni cest a dat
addpath('..\Support'); %slozka je v nadrazenem adresari, proto ..\

load yankees.txt

year = yankees(:,1);
PCT = yankees(:,2); %winning percentage in year t
OBP = yankees(:,5); %team on-base percentage in year t
SLG = yankees(:,6); %team slugging average in year t
ERA = yankees(:,8); %team earned run average in year t

tic %pocitadlo casu behu

%% Odhad Gibbsovym vzorkovacem za predpokladu autokorelace (AR(1))
%Datove matice
 y_full = PCT;
 X_full = [ones(size(y_full)) OBP SLG ERA];

[n,k] = size(X_full); %pocet pozorovani "n" a pocet parametru k

 S1 = 25000;
 S0 = 25000+1;
 S = S0+S1;

%neinformativni priory pro beta a h (nejsou tedy uvadeny)
 nu_0 = 0;
 inv_V_0 = zeros(k); %V0^-1

 p = 1; %rad AR procesu
 c = 1;
 rho_0 = 0;
 Vrho_0 = c*eye(p);

 beta = zeros(k,S);
 h = zeros(1,S);
 rho = zeros(p,S);

%% vytvoreni zpozdenych promennych + adekvatni zkraceni puvodni rady
y = y_full(2:end,:);
y_1 = y_full(1:end-1,:);
X = X_full(2:end,:);
X_1 = X_full(1:end-1,:);


%% Gibbsuv vzorkovac
%pocatecni hodnoty
 h(1,1) = 1;
 count_post = 0; %pocet zamitnutych rho; lze pak vyuzit v ramci 
%Savage-Dickey pro vypocet integracni konstanty (pro p=1 lze i analyticky)

 fprintf('Posteriorni simulace\n')
 fprintf('0%%     50%%      100%%\n')
 fprintf('|')
 %pomocna promenna pro graficke znazorneni prubehu
 pom_step = 0.05; %0.05 = posun po 5% 
 pom_graph = pom_step;  
 
for s=2:S
    ystar = y-rho(:,s-1)*y_1;
    Xstar = X-rho(:,s-1)*X_1;
    
    V_1 = inv(h(:,s-1)*Xstar'*Xstar);
    beta_1 = V_1*(h(:,s-1)*Xstar'*ystar);
    beta(:,s)=beta_1+norm_rnd(V_1);

    nu_1 =n;
    h_1 = (1/nu_1*((ystar-Xstar*beta(:,s))'*(ystar-Xstar*beta(:,s))))^(-1);
    h(:,s)=gamm_rnd_Koop(h_1,nu_1,1);

    %vypocet a sestrojeni matice E (nahodnych slozek a jejich zpozdeni)
    eps_full = y_full-X_full*beta(:,s);
    eps = eps_full(1+p:end);
    E = eps_full(1:end-1);

    Vrho_1 = inv(inv(Vrho_0)+h(:,s)*E'*E);
    rho_1 = Vrho_1*(inv(Vrho_0)*rho_0+h(:,s)*E'*eps);

    testRho = 0; %vyber rho jen ze stacionarni oblasti
    while testRho==0
     rho(:,s) = rho_1+norm_rnd(Vrho_1);
     count_post = count_post+1;
     if abs(rho(:,s))<1
        testRho=1;
        count_post = count_post-1;
     end
    end
    
    %graficke zobrazeni prubehu po 5 %
    if s/S>=pom_graph    
     fprintf('|');
     pom_graph = pom_graph+pom_step;
    end
 end  
 fprintf('\n')

%vyhozeni prvnich S0 vzorku
beta(:,1:S0) = [];
h(:,1:S0) = [];
rho(:,1:S0) = [];

%Vypocet posteriornich charakteristik
 b_mean = mean(beta,2); %posteriorni stredni hodnota 
 h_mean = mean(h);
 rho_mean = mean(rho);
 b_var = mean(beta.^2,2)-b_mean.^2; %posteriorni rozptyl
 h_var = mean(h.^2)-h_mean.^2;
 rho_var = mean(rho.^2)-rho_mean.^2;


%KONVERGENCNI DIAGNOSTIKY - Gewekova CD pomoci funkce Geweke.m
theta = [beta' h' rho'];
res_converg = Geweke(theta);

%HPDI pro beta
 hperc=0.95; %HPDI
 HPDI_beta = (quantile(beta',[0.05,0.95]))';
 HPDI_h = quantile(h,[0.05,0.95]);
 HPDI_rho = quantile(rho,[0.05,0.95]); %jednorozmerny vektor ... neni potreba transpozice

 t = toc;

%Prezentace vysledku, grafy
fprintf('Prior and posterior means for beta (std. deviations in parentheses)\n');
fprintf('                Prior       Posterior      NSE         Geweke CD             %3.1f%% HPDI \n',hperc*100);
fprintf('   ------------------------------------\n');

figure
for ii=1:k
fprintf('Beta%1.0f           neinf.   %12.4f   %10.4f   %10.4f   [%10.4f, %10.4f]\n',[ii b_mean(ii) res_converg.NSE(ii) res_converg.CD(ii) HPDI_beta(ii,:)]);
fprintf('                   (N.A.) (%12.4f)\n', sqrt(b_var(ii)));

    subplot(3,2,ii)
    hist(beta(ii,:),50)
    xlabel(['\beta_' num2str(ii)])
end

fprintf('h               neinf.   %12.4f   %10.4f   %10.4f   [%10.4f, %10.4f]\n',[h_mean res_converg.NSE(k+1) res_converg.CD(k+1) HPDI_h]);
fprintf('                   (N.A.) (%12.4f)\n', sqrt(h_var));
    figure
    hist(h,50)
    xlabel('h')

figure
for ii=1:p
fprintf('rho_%1.0f     %12.4f   %12.4f   %10.4f   %10.4f   [%10.4f, %10.4f]\n',[ii rho_0 rho_mean(ii) res_converg.NSE(k+1+ii) res_converg.CD(k+1+ii) HPDI_rho(ii,:)]);
fprintf('           (%12.4f)   (%12.4f)\n',[sqrt(Vrho_0(ii,ii)) sqrt(rho_var(ii))]);

    subplot(3,2,ii)
    hist(rho(ii,:),50)
    xlabel(['\rho_' num2str(ii)])
end



fprintf('   ------------------------------------\n');
fprintf('Number of replications: %12.4f total \n',S-1);
fprintf('                        %12.4f discarded \n',S0-1);
fprintf('   ------------------------------------\n');

fprintf('Time elapsed (s): %12.2f \n',t);
fprintf('***************************************************************************\n\n\n');
