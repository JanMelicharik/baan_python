function fprior = CES_prior(beta,h,b0,V0,h0,nu0)

k = length(b0);

pb = 1/(2*pi)^(k/2)*det(V0)^(-1/2)*exp(-1/2*(beta-b0)'*inv(V0)*(beta-b0));

cg = (2*h0/nu0)^(nu0/2)*gamma(nu0/2);
ph = cg^(-1)*h^((nu0/2)/2)*exp(-(h*nu0)/(2*h0));

fprior = pb*ph;