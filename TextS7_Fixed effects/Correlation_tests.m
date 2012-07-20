% Filename: Correlation_tests.m
% Author: Andrew Zammit Mangion
% Date: April 2012
%
% Description: Descriptive plots of various covariates against intensity
% estimates in the AWD. For user convenience, some maps have been digitised
% a priori, if the original files are required please contact the first
% author.
%
% Requirements: MAPPING TOOLBOX
%               Correlation_examine.m (plots relevant quantities)
%               Find_pop_density.m (finds population density in provinces)
%               Spatial_Nonparametric.mat (AWD spatial intensity map, can be generated with Spatial_nonparametric.m)
%               IN.mat (digital mask of AFG interior)
%               Dist_to_city.mat (distances to major cities - digitised)
%               Dist_to_Pak.mat (distances to Pakistan - digitised)
%               admin1_poly_32.shp (AFG border)
%               admin2_poly_32.shp (AFG provincial borders)
%               20120408123907_1134450312.tif (GEOTIF for elevation)
%               ../Shapefiles/ (Folder containing population shapefiles)
%               ../Common_functions/ (Folder containing some useful functions)


close all
clear all

addpath('../Common_functions')
load('Spatial_Nonparametric')   % Load nonparametric spatial intensity map

min_intensity = 0.01;   % Treat extremely low intensities as 0 intensity
ds = mean(diff(s1));

% Load shapefiles
Countrybounds = shaperead('../Shapefiles/admin1_poly_32.shp','UseGeoCoords',true);
Provbounds = shaperead('../Shapefiles/admin2_poly_32.shp','UseGeoCoords',true);
load('IN') % Load Afghanistan interior mask

% Nonparametric spatial plot
figure('Position',[100,100,700,400])
c = flipud(hot);
GRID = reshape(GRID(:).*IN,size(GRID)); %Only consider places within borders in subsequent analysis

% Plot spatial intensity map
DrawAFGmap
surfm(S2,S1,log(GRID))
shading interp
colormap(c)
colorbar
caxis([0 12])
setm(ax,'Grid','on','Frame','on','meridianlabel','on','parallellabel','on','FontSize',14)

% Check effect of elevation
figure('Position',[100,100,700,400])
DrawAFGmap
[el_im, cmap, el_bounds] = geotiffread('20120408123907_1134450312');
s1_el = linspace(el_bounds(1,1),el_bounds(2,1),size(el_im,2));
s2_el = linspace(el_bounds(1,2),el_bounds(2,2),size(el_im,1));
[S1_el,S2_el] = meshgrid(s1_el,s2_el);
el_im = double(flipud(el_im));

el_cont =  @(s1,s2) interp2(S1_el,S2_el,el_im,s1,s2);
el_im_reduced = el_cont(S1,S2);
surfm(S2,S1,reshape(el_im_reduced(:).*IN,size(el_im_reduced)));
colormap(c)
colorbar
setm(ax,'Grid','on','Frame','on','meridianlabel','on','parallellabel','on','FontSize',14)
text(0.185,0.545,'elevation (m)','rotation',90,'FontSize',16)
set(gcf,'PaperPositionMode','auto')
print -dpng -r300 Elevation_map.png
 
figure('Position',[100,100,700,400])
Correlation_examine(900:100:5000,el_im_reduced,GRID,min_intensity) % Check only between 900 and 5000m
xlabel('elevation (m)'); ylabel('log intensity');
set(gcf,'PaperPositionMode','auto')
print -dpng -r300 Elevation_corr.png

% Check effect of change in elevation
figure('Position',[100,100,700,400])
DrawAFGmap
el_gradient = atand(gradient(el_im,S2_el,S1_el)/110000);

el_grad_cont =  @(s1,s2) interp2(S1_el,S2_el,el_gradient,s1,s2);
el_grad_reduced = el_grad_cont(S1,S2);
surfm(S2,S1,reshape(abs(el_grad_reduced(:)).*IN,size(el_grad_reduced)));
colormap(c)
colorbar
setm(ax,'Grid','on','Frame','on','meridianlabel','on','parallellabel','on','FontSize',14)
text(0.185,0.515,'gradient (degrees)','rotation',90,'FontSize',16)
set(gcf,'PaperPositionMode','auto')
print -dpng -r300 Gradient_map.png


figure('Position',[100,100,700,400])
Correlation_examine(0:0.5:20,el_grad_reduced,GRID,min_intensity) % Only examine between 0 and 20 degrees gradient
xlabel('gradient (degrees)'); ylabel('log intensity')
set(gcf,'PaperPositionMode','auto')
print -dpng -r300 Gradient_corr.png

% Check effect of population density
Find_pop_density;
figure('Position',[100,100,700,400])
DrawAFGmap
surfm(S2,S1,reshape(Pop_density(:).*IN,size(GRID)));
caxis([0 100]); colormap(c); colorbar
setm(ax,'Grid','on','Frame','on','meridianlabel','on','parallellabel','on','FontSize',14)
text(0.185,0.51,'population density (pax/km^2)','rotation',90,'FontSize',16)
set(gcf,'PaperPositionMode','auto')
print -dpng -r300 PopDensity_map.png


figure('Position',[100,100,700,400])
 Correlation_examine([1,2,5,10,20,30,40,50,100,200,500,1000,15000],Pop_density,GRID,min_intensity)
xlabel('population density (pax/km^2)'); ylabel('log intensity')
set(gca,'xscale','log')
set(gcf,'PaperPositionMode','auto')
print -dpng -r300 PopDensity_corr.png

% Check effect of distance to major city
% Find_dist_to_city;  % Uncomment this line to generate Dist_to_city.mat
load('Dist_to_city')
figure('Position',[100,100,700,400])
DrawAFGmap
surfm(S2,S1,reshape(Dist_to_city(:).*IN,size(GRID)));
colormap(c)
colorbar
caxis([0 200])
setm(ax,'Grid','on','Frame','on','meridianlabel','on','parallellabel','on','FontSize',14)
plotm(Filtered_cities(:,2),Filtered_cities(:,1),'kx')
text(0.185,0.51,'distance to major city (km)','rotation',90,'FontSize',16)
set(gcf,'PaperPositionMode','auto')
print -dpng -r300 Dist_city_map.png

figure('Position',[100,100,700,400])
Correlation_examine(0:10:200,Dist_to_city,GRID,min_intensity)
xlabel('distance to major city (km)'); ylabel('log intensity')
set(gcf,'PaperPositionMode','auto')
print -dpng -r300 Dist_city_corr.png

% Distance to Pakistan border
% Find_dist_to_Pak();   % Uncomment this line to generate Dist_to_Pak.mat
figure('Position',[100,100,700,400])
load('Dist_to_Pak')
DrawAFGmap
surfm(S2,S1,reshape(Dist_to_Pak(:).*IN,size(GRID)));
colormap(c)
colorbar
caxis([0 200])
setm(ax,'Grid','on','Frame','on','meridianlabel','on','parallellabel','on','FontSize',14)
text(0.185,0.51,'distance to Pakistan (km)','rotation',90,'FontSize',16)
set(gcf,'PaperPositionMode','auto')
print -dpng -r300 DistPak_map.png

figure('Position',[100,100,700,400])
Correlation_examine(0:5:320,Dist_to_Pak,GRID,min_intensity)
xlabel('distance to Pakistan (km)'); ylabel('log intensity')
set(gcf,'PaperPositionMode','auto')
print -dpng -r300 DistPak_corr.png

