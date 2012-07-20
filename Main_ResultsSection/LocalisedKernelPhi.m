function [phi] = LocalisedKernelPhi(s1,s2,mu1,mu2,sigma21,sigma22)

% Evaluate the CGRBF centred on (mu1,mu2) with stds (sigma21,sigma22) on
% meshgrid matrices (s1,s2).

s1 = s1(1,:);
s2 = s2(:,1);
J = size(s1,2);
nx = length(mu1);
phi = zeros(J,J,nx);
for i = 1:nx
    beta1 = sqrt(pi/sigma21(i));
    beta2 = sqrt(pi/sigma22(i));
    l1 = 2*pi/beta1;
    l2 = 2*pi/beta2;
    [~,ilow1] = min(abs(s1 - (mu1(i) - l1)));   %Find limit of basis
    [~,ihigh1] = min(abs(s1 - (mu1(i) + l1)));  %Find limit of basis
    [~,ilow2] = min(abs(s2 - (mu2(i) - l2)));   %Find limit of basis
    [~,ihigh2] = min(abs(s2 - (mu2(i) + l2)));  %Find limit of basis
    s_on1 = s1(ilow1:ihigh1);
    s_on2 = s2(ilow2:ihigh2);
    [Delta1,Delta2] = meshgrid(beta1*abs(s_on1 - mu1(i)),beta2*abs(s_on2 - mu2(i)));
    phi(ilow2:ihigh2,ilow1:ihigh1,i) = ((2*pi - Delta1).*(1 + cos(Delta1)/2) + 3/2*sin(Delta1))./(3*pi).*((2*pi - Delta2).*(1 + cos(Delta2)/2) + 3/2*sin(Delta2))./(3*pi);
end
