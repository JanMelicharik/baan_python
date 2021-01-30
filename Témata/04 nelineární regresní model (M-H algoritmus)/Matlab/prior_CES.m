function res = prior_CES(gamma,h,gamma_0,V_0,h_0,nu_0)
%Neomezena apriorni hustota (bez restrikci na kladnost gamma)
k = length(gamma_0);
prior_gamma = 1/(2*pi)^(k/2)*det(V_0)^(-1/2)...
    *exp(-1/2*(gamma-gamma_0)'*inv(V_0)*(gamma-gamma_0));
c_h = (2*h_0/nu_0)^(nu_0/2)*my_gamma(nu_0/2);
prior_h = 1/c_h*h^((nu_0-2)/2)*exp(-h*nu_0/(2*h_0));
res = prior_gamma*prior_h;
end

