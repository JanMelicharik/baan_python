function plt_nlrm_ncp(results)
% funkce pro vykresleni apriornich a posteriornich marginalnich hustot

k=length(results.b0);

%data pro vykresleni apriorniho marginalniho rozdeleni pro beta (t-rozdeleni) 

sig0 = diag(results.s02*results.V0);
sig1 = diag(results.s12*results.V1);
figure

for num=1:k
alower=min([results.b0(num),results.b1(num)])-2*max([1,results.bstd0(num),results.bstd1(num)]);
aupper=max([results.b0(num),results.b1(num)])+2*max([1,results.bstd0(num),results.bstd1(num)]);
incr=0.01;
bplot=alower:incr:aupper;
b_plotpri = my_tpdf(bplot,results.b0(num),sig0(num),results.nu0);
%data pro vykresleni posteriorniho marginalniho rozdeleni pro beta (t-rozdeleni)
b1_plotpost = my_tpdf(bplot,results.b1(num),sig1(num),results.nu1);

subplot(round(sqrt(k)),ceil(sqrt(k)),num)     
plot(bplot,b_plotpri./sum(b_plotpri),'--',bplot',b1_plotpost./sum(b1_plotpost))
title(['Obrazek: Marginalni Prior a Posterior pro \beta_',num2str(num)])
legend('Prior','Posterior')
xlabel(['\beta_',num2str(num)])
ylabel('Hustota pravdepodobnosti')

end





function b = my_tpdf(y,mu,s2,nu)
%for univariate t with arguments y,mu,s2,nu 
%evaluate the t density at for each element of y


c = .5*log(pi)+gammaln(.5*nu)-.5*nu*log(nu)-gammaln(.5*(nu+1));

b = -c-.5*log(s2)-.5*(nu+1)*log(nu+(y-mu).^2./s2);

b=exp(b);