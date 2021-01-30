function res = lik_CES(y,X,gamma,h)
%Verohodnostni funkce pro CES produkcni funkci
fx = gamma(1)*(gamma(2)*X(:,2).^gamma(4)+gamma(3)*X(:,3).^gamma(4)).^(1/gamma(4));
N = length(y);
log_res = -N/2*log(2*pi)+N/2*log(h)-h/2*(y-fx)'*(y-fx);
res = exp(log_res);
end

