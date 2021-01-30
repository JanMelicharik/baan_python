function fpost = nu_lambda_post(nu_lambda,nu_lambda0,lambda)
N=length(lambda);


eta = 1/nu_lambda0+.5*sum(-log(lambda)+lambda);


fpost=N*nu_lambda/2*log(nu_lambda/2)-N*gammaln(nu_lambda/2)-eta*nu_lambda;
