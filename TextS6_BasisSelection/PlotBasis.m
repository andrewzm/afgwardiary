function y = PlotBasis()

%--------------------------------------------------------------------
% Filename: AFGBasisSelectionHighRes.m
% Authors: Andrew Zammit Mangion
% Date: March 14 2011
% Description: Plots the basis functions on AFG map
%
% Requires: Mapping toolbox
%           ../Shapefiles
%           ../Common_functions (for drawing map)
%
% Generates: 
%--------------------------------------------------------------------


load('AFGBasis')

Countrybounds = shaperead('../Shapefiles/admin1_poly_32.shp','UseGeoCoords',true);
shift = [57.5,28.5];
scale = [2,3.3];

 %Upsample
s = linspace(0,36,201);
ds = mean(diff(s));
[s1,s2] = meshgrid(s);
Basis.phi = LocalisedKernelPhi(s1,s2,Basis.mu1,Basis.mu2,Basis.sigma2,Basis.sigma2);
strue1 = s./scale(1) + shift(1);
strue2 = s./scale(2) + shift(2);
[s1true,s2true] = meshgrid(strue1,strue2);
 
figure('Position',[100,100,600,400])
h = axes('Position',[0 0 1 1],'Visible','off');

axes('Position',[.1 .1 .9 .9])
DrawAFGmapnoprov
setm(ax,'Grid','on','Frame','on','meridianlabel','on','parallellabel','on','FontSize',14)
hold on
for i = 1:Basis.nx
    [C,h] = contourm(s2true,s1true,Basis.phi(:,:,i),[1 exp(-1^2/2)],'r');
    plotm(Basis.mu2(i)./scale(2) + shift(2),Basis.mu1(i)./scale(1) + shift(1),'kx')
%     set(h,'LevelStep',0.5)   
end
xlabel('Longitude','Fontsize',14)
ylabel('Latitude','Fontsize',14)
set(gca,'XTick',[])
set(gca,'YTick',[])

set(gcf,'PaperPositionMode','auto')
print -dpng -r600  BasisPlacement.png

load('../AfghanModelHighResWholePrec_Covariates');
figure('Position',[100,100,600,400])
h = axes('Position',[0 0 1 1],'Visible','off');
axes('Position',[.1 .1 .9 .9])
DrawAFGmapnoprov
setm(ax,'Grid','on','Frame','on','meridianlabel','on','parallellabel','on','FontSize',14)
hold on
for k = 1:length(Weeknum)-1
    plotm(spikeAllWeekData(k).Coords(:,2)./scale(2) + shift(2),spikeAllWeekData(k).Coords(:,1)./scale(1) + shift(1),'r.')
end
xlabel('Longitude','Fontsize',14)
ylabel('Latitude','Fontsize',14)
set(gca,'XTick',[])
set(gca,'YTick',[])

xlabel('Lon')
ylabel('Lat')

set(gcf,'PaperPositionMode','auto')
print -dpng -r600  ../../Figures/AllPoints.png

function [phi] = LocalisedKernelPhi(s1,s2,mu1,mu2,sigma21,sigma22)

% Evaluate the CGRBF centred on (mu1,mu2) with stds (sigma21,sigma22) on
% vectors (s1,s2).

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
  [~,ilow1] = min(abs(s1 - (mu1(i) - l1)));  %Find centre of kernel
  [~,ihigh1] = min(abs(s1 - (mu1(i) + l1)));  %Find centre of kernel
  [~,ilow2] = min(abs(s2 - (mu2(i) - l2)));  %Find centre of kernel
  [~,ihigh2] = min(abs(s2 - (mu2(i) + l2)));  %Find centre of kernel
  s_on1 = s1(ilow1:ihigh1);
  s_on2 = s2(ilow2:ihigh2);
  [Delta1,Delta2] = meshgrid(beta1*abs(s_on1 - mu1(i)),beta2*abs(s_on2 - mu2(i)));
  phi(ilow2:ihigh2,ilow1:ihigh1,i) = ((2*pi - Delta1).*(1 + cos(Delta1)/2) + 3/2*sin(Delta1))./(3*pi).*((2*pi - Delta2).*(1 + cos(Delta2)/2) + 3/2*sin(Delta2))./(3*pi);
end

