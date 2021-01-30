function b = my_tpdf(y,mu,s2,nu)
%for univariate t with arguments y,mu,s2,nu 
%evaluate the t density at for each element of y

c = .5*log(pi)+gammaln(.5*nu)-.5*nu*log(nu)-gammaln(.5*(nu+1));

b = -c-.5*log(s2)-.5*(nu+1)*log(nu+(y-mu).^2./s2);

b=exp(b);






