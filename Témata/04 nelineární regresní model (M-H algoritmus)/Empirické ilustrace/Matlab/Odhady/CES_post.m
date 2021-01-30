function fpost = CES_post(gamma,y,X)

f = gamma(1)+(gamma(2)*X(:,2).^gamma(4)+gamma(3)*X(:,3).^gamma(4)).^(1/gamma(4));

fpost=-log(((y-f)'*(y-f)))*(-length(y)/2);
