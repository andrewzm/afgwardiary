function [Aest,lkernelest,g] = PairCrossCorrFunc(r,spikes1,spikes2,s)

% Filename: PairCrossCorrFunc.m
% Author: Andrew Zammit Mangion
% Date: February 2011
% Description: Find the PCCF of two spatial point processes
%
% Inputs: r (vector of radial distances to consider)
%         spikes1 (first point process)
%         spikes2 (second point process)
%         length of square describing domain
%
% Outputs: Aest (amplitude of fitted GRBF)
%          lkernelest (kernel width of fitted GRBF)
%          g (Nonparametric estimate of the PCCF)


% Initialize
Ntot1 = length(spikes1);
Ntot2 = length(spikes2);
g = zeros(1,length(r));
srange = s(end);
fintest1 = @(s) Ntot1/srange^2*ones(size(s,1),1);
fintest2 = @(s) Ntot2/srange^2*ones(size(s,1),1);
b = 1;
ktype = 'RBF';

% Start estimation over radii points in vector r
for i = 1:length(r) 
    mysum = 0;
    for j = 1:length(spikes1)
        radprop = findpercentage(repmat(spikes1(j,:),Ntot2,1),spikes2,[s(1) s(1)],[s(end) s(end)],r(i),b);
        w = 1./radprop;
        mysum = mysum + sum(kb1(repmat(r(i),Ntot2,1) - hypot(repmat(spikes1(j,1),Ntot2,1) - spikes2(:,1),repmat(spikes1(j,2),Ntot2,1) - spikes2(:,2)),min(b,r(i)))...
            ./(repmat(fintest1(spikes1(j,:)),Ntot2,1).*fintest2(spikes2)).*w);
    end
    g(i) = 1/(2*pi*r(i))*mysum/(s(end)^2);
end

%Estimation
epsilon = r(1); %lower limit
a0 = 5; %upper limit
[Aest,lkernelest] = MethodContrast(r,g,ktype,epsilon,a0);

function w = findpercentage(x,y,s1,s2,a,b)
% Correction for circle not being completely within square boundary.
w = zeros(size(x,1),1);
Npoints = 50;
for i = 1:size(x,1)
   r = hypot(x(i,1) - y(i,1),x(i,2) - y(i,2)); %find radius
   if (r == 0) || (abs(r - a) > b) %If point is itself or "dr" outside kernel BW set it v. high so 1/w = 0; 
       w(i) = 1e6;
   else 
       d1 = min(abs(s1(1) - x(i,1)),abs(s2(1) - x(i,1)));
       d2 = min(abs(s1(2) - x(i,2)),abs(s2(2) - x(i,2)));
       if r < d1 && r < d2 %circle inside box
           w(i) = 1;
       else %Find it numerically
           THETA=linspace(0,2*pi,Npoints);
           RHO=ones(1,Npoints)*r;
           %Find circle coordinates
           [X,Y] = pol2cart(THETA,RHO);
           X=X+x(i,1);
           Y=Y+x(i,2);
           %Count
           w(i) = sum(((X > s1(1)) & (X < s2(1)))&( (Y > s1(2)) & (Y < s2(2))))/Npoints;
           if (w(i) == 0) 
               w(i) = 1/Npoints; 
           end
       end
       %DEBUG
%        circle(x(1,:),r,1000)
%        rectangle('Position',[0,0,s2(1),s2(2)])
   end
end

function y = kb1(x,b)
% Epanechnikov kernel

y = 3/(4*b)*(1 - x.^2/b^2).*(abs(x(:,1))<b);
