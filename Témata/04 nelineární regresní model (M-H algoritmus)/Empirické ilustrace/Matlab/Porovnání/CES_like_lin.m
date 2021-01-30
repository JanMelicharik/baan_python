function flike = CES_like_lin(gamma,h,y,X)

f = gamma(1)*X(:,1)+gamma(2)*X(:,2)+gamma(3)*X(:,3);

N=length(y);
flike = h^(N/2)/(2*pi)^(N/2)*exp(-h/2*(y-f)'*(y-f));