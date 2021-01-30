function b = log_mvtpdf(y,mu,Sigma,nu)
%for multivariate t with arguments y,mu,Sigma,nu 
%evaluate the t density

k=length(y);

c = k*.5*log(pi)+gammaln(.5*nu)-.5*nu*log(nu)-gammaln(.5*(nu+k));

b = -c-.5*log(det(Sigma))-.5*(nu+k)*log(nu+(y-mu)'*inv(Sigma)*(y-mu));
