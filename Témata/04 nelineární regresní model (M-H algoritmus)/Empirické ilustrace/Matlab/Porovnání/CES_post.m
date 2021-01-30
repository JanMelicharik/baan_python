function fpost = CES_post(gamma,y,X,h,gamma0,V0)

f = gamma(1)*X(:,1)+(gamma(2)*X(:,2).^gamma(4)+gamma(3)*X(:,3).^gamma(4)).^(1/gamma(4));

%logaritmus posteriorni hustoty (jadra) krat -1 (vstupuje do minimalizace)
fpost=-((-h/2*(y-f)'*(y-f))+(-1/2*(gamma-gamma0)'*inv(V0)*(gamma-gamma0)));