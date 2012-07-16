function [A,l] = MethodContrast(r,g,ktype,epsilon,a0)

% Fit a GRBF to the log-PACF or log-PCCF using the method of contrast

sinterp = [epsilon:0.1:a0];
ginterp = interp1(r,g,sinterp);
c_hat = log(ginterp);

beta = fminsearch(@(par) fbeta(par,c_hat,ktype,sinterp),[c_hat(1),1]);
A = beta(1);
l = beta(2);

%Following Moller 1998


function y = fbeta(b,c_hat,ktype,s)

ds = mean(diff(s));
if strcmp(ktype,'RBF') 
   r_beta =  b(1)*exp(-s.^2./b(2));
elseif strcmp(ktype,'EXP') 
   r_beta = b(1)*exp(-abs(s)./b(2));
end
y = ds^2*sum(sum((c_hat - r_beta).^2));