% PlotGrowthMap.m
% Description: Generates Video S1 in paper
% Requires: Mapping toolbox
%           AfghanModelHighResWholePrec_Covariates.mat (results)
%           IN.mat (digital mask of Afghanistan)
%           LocalisedKernelPhi.m (basis function definition)
%           ../Shapefiles
%           inpolygons.m (available from Matlab Exchange)
%           multiprod.m (available from Matlab Exchange)
%           posdist.m (available from Matlab Exchange)
%           mpgwrite (available from Matlab Exchange)
%           month.m (Financial toolbox but can be coded easily)
%
%   If Quantile_data.mat is not present, generate it by uncommenting the
%   line "Find_Quantiles" below. Note: this operation takes a long time.


clear all
close all

studydata = '../AfghanModelHighResWholePrec_Covariates'; 
load(studydata)
mu = Thetainfo(m).muest + Thetainfo(m).varmu/2;


%Upsample Basis
J2 = 201;
s = linspace(0,36,J2);
ds = mean(diff(s));
[s1,s2] = meshgrid(s);
strue1 = s./scale(1) + shift(1);
strue2 = s./scale(2) + shift(2);
[s1true,s2true] = meshgrid(strue1,strue2);
Basis.phi = LocalisedKernelPhi(s1,s2,Basis.mu1,Basis.mu2,Basis.sigma2,Basis.sigma2);
Countrybounds = shaperead('../Shapefiles/admin1_poly_32.shp','UseGeoCoords',true);
Provbounds = shaperead('../Shapefiles/admin2_poly_32.shp','UseGeoCoords',true);
% IN = inpolygons(s1true(:),s2true(:),Countrybounds.Lon',Countrybounds.Lat');
load('IN')

%--------------------------------
%FOURTH FIGURE IN TECHNICAL REPORT
%Growth map estimate
%--------------------------------
figure('Position',[100,100,700,600])
% Getmask2
c = flipud(hot);
% ezmillerCart(shift,scale)
hold on
test = sum(multiprod(Basis.phi,reshape(Thetainfo(m).thetaest,1,1,Basis.nx)),3);
test = test.*(test>0);
test = reshape(test(:).*IN,size(test));

%FACTOR of 6.6 LON Scale = 2, LAT SCALE = 3.3. 100km approx 0.9LAT*1.1LON approx 1
DrawAFGmap
surfm(strue2,strue1,exp(test))
% geoshow(ax,s1truevec',s2truevec',test);
shading interp
colormap(c)
% h = colorbar('EastOutside')
% set(h,'position',[0.68,0.21,0.03,0.25])
% set(h,'YTick',[]);
h = colorbar('SouthOutside')
set(h,'position',[0.5,0.1,0.3,0.03])
labels=get(h,'xticklabel'); % get the y axis labels
for i=1:size(labels,1)
    newlabel = sprintf('%0.1f',(str2num(labels(i,:))-1)*100);
    labels_modif(i,:)=[newlabel '%'];
end
set(h,'xticklabel',labels_modif);

hold off
% caxis([0 0.025])
xlabel('Lon','Fontsize',14)
ylabel('Lat','Fontsize',14)

print -dpng -r300 GrowthMap.png

%--------------------------------
%FIFTH FIGURE IN TECHNICAL REPORT
%Predictability map estimate
%--------------------------------
figure('Position',[100,100,700,600])
% Getmask2
c = flipud(hot);
% ezmillerCart(shift,scale)
hold on
% test = sum(multiprod(Basis.phi,reshape(1./Thetainfo(m).precisionest,1,1,Basis.nx)),3);
test = sum(multiprod(Basis.phi,reshape(1./diag(Thetainfo(m).Meanprecmat),1,1,Basis.nx)),3);
test = test.*(test>0);
test = reshape(test(:).*IN,size(test));

%FACTOR of 6.6 LON Scale = 2, LAT SCALE = 3.3. 100km approx 0.9LAT*1.1LON approx 1
DrawAFGmap
surfm(strue2,strue1,test)

% geoshow(ax,s1truevec',s2truevec',test);
shading interp
colormap(c)
caxis([0.055,0.07])
h = colorbar('EastOutside')
set(h,'position',[0.68,0.21,0.03,0.25])
set(h,'YTick',[]);
hold off
% caxis([0 0.025])
xlabel('Lon','Fontsize',14)
ylabel('Lat','Fontsize',14)

print -dpng -r600 VolatilityMap_unsrt.png

%--------------------------------
%INSETS OF FIGURE 4 IN TECHNICAL REPORT
%Growth map estimate Localised
%--------------------------------

% Find mean, median, lower and upper
% Find_Quantiles;
load('Quantile_data')
load('../AfghanDataAllDay','spikeAll')

%SANGIN: LAT-LON 32.08, 64.83
%100km side BOX Corners at 31.63,64.3 and 32.53,65.35
figure('Position',[100 200 600 250])
lcorner = [64.3,31.63];
ucorner = [65.35,32.53];
Analyzecity(strue1,strue2,lcorner,ucorner,spikeAll,int_lower,int_median,int_upper,ds)
set(gcf,'PaperPositionMode','auto')
print -dpng -r600 Sangin.png

%KUNDUZ: LAT-LON 36.76,68.83
%100km side BOX Corners at 36.31,68.27 and 37.21,69.39
figure('Position',[100 200 600 250])
lcorner = [68.27,36.31];
ucorner = [69.39,37.21];
Analyzecity(strue1,strue2,lcorner,ucorner,spikeAll,int_lower,int_median,int_upper,ds)
set(gcf,'PaperPositionMode','auto')
print -dpng -r600 Kunduz.png

%BAGHLAN: LAT-LON 36.13,69.24
%100km side BOX Corners at 36.31,68.27 and 37.21,69.39
figure('Position',[100 200 600 250])
lcorner = [68.12,35.68];
ucorner = [69.24,36.58];
Analyzecity(strue1,strue2,lcorner,ucorner,spikeAll,int_lower,int_median,int_upper,ds)
set(gcf,'PaperPositionMode','auto')
print -dpng -r600 Baghlan.png

%DARREH-YE-BUM: LAT-LON 35.13,63.46
%100km side BOX Corners at 34.68,62.91 and 35.58,64.01
figure('Position',[100 200 600 250])
lcorner = [62.91,34.68];
ucorner = [64.01,35.58];
Analyzecity(strue1,strue2,lcorner,ucorner,spikeAll,int_lower,int_median,int_upper,ds)
set(gcf,'PaperPositionMode','auto')
print -dpng -r600 Darreh.png

%KABUL: LAT-LON 34.57,69.13
%100km side BOX Corners at 34.12,68.58 and 35.02,69.68
figure('Position',[100 200 600 250])
lcorner = [68.58,34.12];
ucorner = [69.68,35.02];
Analyzecity(strue1,strue2,lcorner,ucorner,spikeAll,int_lower,int_median,int_upper,ds)
set(gcf,'PaperPositionMode','auto')
print -dpng -r600 Kabul.png
