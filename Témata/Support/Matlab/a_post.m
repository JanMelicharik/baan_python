function fpost = a_post(a,b,h,y,X,z)

Omega = diag((ones(size(y))+z*a).^2);

%log(det(Omega)) lze zapsat nasledovne, cimz se vyvarujeme numerickych chyb
%pri pocitani determinantu diagonalni matici obrovskych rozmeru
pom = sum(log(diag(Omega)));

f = y-X*b;


%fpost=-0.5*log(det(Omega))-0.5*h*f'*inv(Omega)*f;
fpost=-0.5*pom-0.5*h*f'*inv(Omega)*f;
