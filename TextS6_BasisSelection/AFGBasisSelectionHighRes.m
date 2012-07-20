function AFGBasisSelectionHighRes()

%--------------------------------------------------------------------
% Filename: AFGBasisSelectionHighRes.m
% Authors: Andrew Zammit Mangion
% Date: March 14 2011
% Description: A principle way for basis function placement using a priori
% frequency-based knowledge
%
% Requires: Mapping toolbox
%           ProcessedWeekData.mat (a curated version of the AWD, no points
% outside Afghanistan and organised into weeks)
%
% Generates: AfgBasis.mat, used as a basis in main VB inference program
%--------------------------------------------------------------------

close all
clear all

%-----------------
% Setup constants
%-----------------

%Setup space
J = 101;
s = linspace(0,36,J);
l=3.5;
smax = s(end)-l; %These are spatial limits on basis representation
smin = l;
[s1,s2] = meshgrid(s,s);

%-------------------------------------------------------------------
% Field setup
%-------------------------------------------------------------------

%Basis decomposition
Basis.nx = 256;
% The following assures beating the Shannon's criterion by approx 1.25
[C1] = meshgrid(linspace(smin,smax,sqrt(Basis.nx)),linspace(smin,smax,sqrt(Basis.nx)));

Basis.mu1 = reshape(C1,1,[]);
Basis.mu2 = reshape(C1',1,[]);
Basis.tau(1:Basis.nx) = 1.8;
Basis.sigma2(1:Basis.nx) = pi./(Basis.tau.^2);
Basis.phi = LocalisedKernelPhi(s1,s2,Basis.mu1,Basis.mu2,Basis.sigma2,Basis.sigma2);

% -----------------------------------------------
% Get points
% -----------------------------------------------


load('ProcessedWeekData')
AllAFGevents = [];
for i = 1:length(spikeAllWeekData)
    AllAFGevents = [AllAFGevents;spikeAllWeekData(i).Coords];
end

FilterBasis(AllAFGevents,Basis,l,shift,scale);

%---------------------------
%Auxiliary functions
%---------------------------
function Basisnew = FilterBasis(spikes,Basis,l,shift,scale)
N = length(spikes);
i = 1;
rad = 1.3;

% Load shapefiles
Countrybounds = shaperead('admin1_poly_32.shp','UseGeoCoords',true);

% For each basis see whether it is significantly outside of Afghanistan. If
% it is remove the basis.
while(i <= Basis.nx)
    if ~inpolygon(Basis.mu1(i)/scale(1)+shift(1),Basis.mu2(i)/scale(2)+shift(2),Countrybounds.Lon',Countrybounds.Lat')
        if p_poly_dist(Basis.mu1(i)/scale(1)+shift(1),Basis.mu2(i)/scale(2)+shift(2),Countrybounds.Lon',Countrybounds.Lat') > 0.4
            Basis.nx = Basis.nx-1;
            Basis.mu1(i) = [];
            Basis.mu2(i) = [];
            Basis.sigma2(i) = [];
            Basis.phi(:,:,i) = [];
            i = i-1;
        end
    end
    i = i+1;
end


% For each basis see how many events are within its scope. If there are
% only a few corresponding to a background of exp(-3.5) or less then remove basis.
i=1;
while(i <= Basis.nx)
    distances = hypot(Basis.mu1(i) - spikes(:,1),Basis.mu2(i) - spikes(:,2));
    %     pntsinbasis(i) = sum(distances < 1.2);
    pntsinbasis(i) = sum(distances < rad);
    %     if pntsinbasis(i) < 24
    if pntsinbasis(i) < exp(-3.5)*313*pi*rad^2 %Bg.noise*numofweeks*area
        Basis.nx = Basis.nx-1;
        Basis.mu1(i) = [];
        Basis.mu2(i) = [];
        Basis.sigma2(i) = [];
        Basis.phi(:,:,i) = [];
        i = i-1;
    end
    i = i+1;
end
Basisnew = Basis;
save('AFGBasis','Basis')

function [phi] = LocalisedKernelPhi(s1,s2,mu1,mu2,sigma21,sigma22)

% Evaluate the CGRBF centred on (mu1,mu2) with stds (sigma21,sigma22) on
% meshgrid arrays (s1,s2).

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
    [~,ilow1] = min(abs(s1 - (mu1(i) - l1)));  %Find limit of kernel
    [~,ihigh1] = min(abs(s1 - (mu1(i) + l1)));  %Find limit of kernel
    [~,ilow2] = min(abs(s2 - (mu2(i) - l2)));  %Find limit of kernel
    [~,ihigh2] = min(abs(s2 - (mu2(i) + l2)));  %Find limit of kernel
    s_on1 = s1(ilow1:ihigh1);
    s_on2 = s2(ilow2:ihigh2);
    [Delta1,Delta2] = meshgrid(beta1*abs(s_on1 - mu1(i)),beta2*abs(s_on2 - mu2(i)));
    phi(ilow2:ihigh2,ilow1:ihigh1,i) = ((2*pi - Delta1).*(1 + cos(Delta1)/2) + 3/2*sin(Delta1))./(3*pi).*((2*pi - Delta2).*(1 + cos(Delta2)/2) + 3/2*sin(Delta2))./(3*pi);
end

