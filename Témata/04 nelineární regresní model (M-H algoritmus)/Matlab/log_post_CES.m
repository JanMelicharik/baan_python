function res = log_post_CES(y,X,gamma,h,gamma_0,V_0)
%res ... logaritmus jadra posteriorni hustoty
%CES funkce
fx = gamma(1)*(gamma(2)*X(:,2).^gamma(4) +gamma(3)*X(:,3).^gamma(4)).^(1/gamma(4));
%logaritmus (5.24) z Koop
res = -1/2*h*(y-fx)'*(y-fx) -1/2*(gamma-gamma_0)'*inv(V_0)*(gamma-gamma_0);
%osetreni nulove hustoty pro parametry gamma
%=> restrikce na kladnost vesch parametru gamma
if min(gamma)<0
    res = -inf;
end

end
