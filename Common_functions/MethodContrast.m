function [A,l] = MethodContrast(r,g,ktype,epsilon,a0)

% Fit a GRBF to the log-PACF or log-PCCF using the method of contrast.
% (Moller 1998)
%
% Inputs: r (radius vector)
%         g (PACF or PCCF)
%         ktype ('RBF' or 'EXP')
%         epsilon (starting radius value for minimization)
%         a0 (final radius value for minimization)
%
%   Outputs: A (amplitude)
%            l (length scale)

sinterp = [epsilon:0.1:a0];
ginterp = interp1(r,g,sinterp);
c_hat = log(ginterp);

% Find minimum, fminsearch is fast enough for 2-D
beta = fminsearch(@(par) fbeta(par,c_hat,ktype,sinterp),[c_hat(1),1]);
A = beta(1);
l = beta(2);



function y = fbeta(b,c_hat,ktype,s)
% Function to minimize

ds = mean(diff(s));
if strcmp(ktype,'RBF') 
   r_beta =  b(1)*exp(-s.^2./b(2));
elseif strcmp(ktype,'EXP') 
   r_beta = b(1)*exp(-abs(s)./b(2));
end
y = ds^2*sum(sum((c_hat - r_beta).^2));