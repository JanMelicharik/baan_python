function results = my_NLRM(y,X,beta_0,V_0,h_0,nu_0)
%Funkce pro bayesianskou analyzu NLRM (prirozene konjugovana apriorni hustota)
% y ... vektory vysvetlujici promenne N x 1
% X ... matice pl�nu N x k
% beta_0, V_0, h_0, nu_0 ... apriorn� hyperparametry
%                            z p(beta,h)~NG(beta_0,V_0,h_0,nu_0)
% beta_1, V_1, h_1, nu_1 ... aposteriorn� hyperparametry
%                            z p(beta,h|y)~NG(beta_0,V_0,h_0,nu_0)
% b0_cov, b0_std, h0_std ... apriorn� kovarian�n� matice
%			     a vektory apriorn�ch sm�rodatn�ch odchylek parametr�
% b1_cov, b1_std, h1_std ... posteriorn� kovarian�n� matice
%			     a vektory posteriorn�ch sm�rodatn�ch odchylek parametr�
% log_ML 		 ... logaritmus margin�ln� v�rohodnosti modelu

%% Po��t�n� charakteristik apriorn�ch hustot
 %apriorn� kovarian�n� matice pro beta
b0_cov = nu_0*h_0^(-1)/(nu_0-2)*V_0;	% analogie (3.16) z Koop (2003) 
b0_std = sqrt(diag(b0_cov));
 %apriorni sm. odchylka pro h
h0_std = sqrt(2*h_0/nu_0);		% analogie odmocniny z (3.19) z Koop (2003) 

%% Po��t�n� aposteriorn�ch hyperparametr�
 N = length(y); %pocet pozorovani
  %odhady OLS (ordinary least squares)
 b_OLS = (X'*X)^-1*X'*y;	% (3.5) z Koop (2003) 
 nu_OLS = N-size(X,2);		% (3.4) z Koop (2003) 
 s2_OLS = 1/nu_OLS*(y-X*b_OLS)'*(y-X*b_OLS);	% (3.6) z Koop (2003) 
 
 %Aposteriorn� hyperparametry
 V_1 =inv(inv(V_0)+X'*X); % (3.10)		% (3.9) z Koop (2003) 
 beta_1 = V_1*(inv(V_0)*beta_0+X'*X*b_OLS); 	% (3.11) z Koop (2003) 
 nu_1 = nu_0+N; 				% (3.12) z Koop (2003) 
 h_1 = nu_1*(nu_0*1/h_0 + nu_OLS*s2_OLS +...
     (b_OLS-beta_0)'*inv(V_0+inv(X'*X))*(b_OLS-beta_0))^-1;
			% (3.13) z Koop (2003), kdy h_1 = s2_1^-1
		 

 % Aposteriorn� kovarian�n� matice a sm�rodatn� odchylky
 b1_cov = nu_1*h_1^(-1)/(nu_1-2)*V_1;	% (3.16) z Koop (2003) 
 b1_std = sqrt(diag(b1_cov));

 %Aposteriorn� sm�rodatn� odchylka pro p�esnost chyby h
 h1_std = sqrt(2*h_1/nu_1);		% odmocnina z (3.19) z Koop (2003) 
 
 %logaritmus marginaln� v�rohodnosti
  %log konstanty a margin�ln� v�rohodnosti
  log_c = gammaln(nu_1/2)+nu_0/2*log(nu_0/h_0)...
      -gammaln(nu_0/2)-N/2*log(pi);	% logaritmus (3.35) z Koop (2003) 
  log_ML = log_c+1/2*(log(det(V_1))-log(det(V_0))) ...
      -nu_1/2*log(nu_1/h_1);		% logaritmus (3.34) z Koop (2003) 
 
  %% V�stup funkce do prom�nn� results
   %p�vodn� data
   results.y = y;
   results.X = X;
   
   %apriorn� hyperparametry a kovariance
   results.beta_0 = beta_0;
   results.h_0 = h_0;
   results.V_0 = V_0;
   results.nu_0 = nu_0;
   results.b_OLS = b_OLS;
   results.nu_OLS = nu_OLS;
   results.s2_OLS = s2_OLS;
   
   results.b0_cov = b0_cov;
   results.b0_std = b0_std;
   results.h0_std = h0_std;
   
   %aposteriorn� hyperparametry a kovariance
   results.beta_1 = beta_1;
   results.h_1 = h_1;
   results.V_1 = V_1;
   results.nu_1 = nu_1;
   
   results.b1_cov = b1_cov;
   results.b1_std = b1_std;
   results.h1_std = h1_std;
   
   %logaritmus margin�ln� v�rohodnosti
   results.log_c = log_c
   results.log_ML = log_ML;
   
end

