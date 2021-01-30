function HPDI=HPDI_nlrm_ncp(results,pdi)
% vypocet intervalu najvyssej posteriornej hustoty na zaklade posteriornych
% vysledkov
% results ... posteriorne vysledky
% pdi ... nastavenie parametru HPDI (pre 95% HPDI, pdi=0.95)
k = length(results.b1);
HPDI=zeros(k,2);


for i = 1:k
    HPDI(i,1) = results.b1(i)-tinv(1-(1-pdi)/2,results.nu1)*sqrt(results.covb1(i,i));
    HPDI(i,2) = results.b1(i)+tinv(1-(1-pdi)/2,results.nu1)*sqrt(results.covb1(i,i));
end

