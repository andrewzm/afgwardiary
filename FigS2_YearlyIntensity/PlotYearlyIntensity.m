% PlotYearlyIntensity.m
% Description: Plots Fig S2 in paper
% Requires: Mapping toolbox
%           AfghanModelHighResWholePrec_Covariates.mat (results)+
%           IN.mat (digital mask of Afghanistan)
%           LocalisedKernelPhi.m (basis function definition)
%           ../Shapefiles
%           inpolygons.m (available from Matlab Exchange)
%           multiprod.m (available from Matlab Exchange)
%           posdist.m (available from Matlab Exchange)
%           month.m (Financial toolbox but can be coded easily)


clear all
close all

studydata = '../AfghanModelHighResWholePrec_Covariates'; 
load(studydata)
mu = Thetainfo(m).muest + Thetainfo(m).varmu/2;


%Upsample Basis
s = linspace(0,36,201);
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

% Find population density
Afg_pop = shaperead('../Shapefiles/afg_population_CSO_2011_2012','UseGeoCoords',true);
Pop_density = zeros(size(s1));
for i = 1:length(Afg_pop)
X = inpolygons(s1true(:),s2true(:),Afg_pop(i).Lon',Afg_pop(i).Lat');
Pop_density(X == 1) = Afg_pop(i).T_P_BS/ Afg_pop(i).AREA_SQKM;
end

% Find distance to city
Afg_cities = shaperead('../Shapefiles/07_03_settlements','UseGeoCoords',true);
Filtered_cities = 0; j =1;
for i = 1:length(Afg_cities)
   if Afg_cities(i).TYPE == 2 ||  Afg_cities(i).TYPE == 3
       Filtered_cities(j,1) = Afg_cities(i).LON;
       Filtered_cities(j,2) = Afg_cities(i).LAT;
       j = j+1;
   end
end
Dist_to_city = zeros(size(s1true));
for i = 1:size(s1true,1)
    for j = 1:size(s1true,2)
        Dist_to_city(i,j) = min(posdist(s1true(1,j),s2true(i,1),Filtered_cities(:,1),Filtered_cities(:,2)));
    end
end

%--------------------------------------
%Intensity estimates
%--------------------------------------
figure('Position',[10,10,1200,650])
c = flipud(hot);
Covariates = Thetainfo(m).best(1)*Pop_density + Thetainfo(m).varb(1)*Pop_density.^2/2 ...
                 + Thetainfo(m).best(2)*Dist_to_city + Thetainfo(m).varb(2)*Dist_to_city.^2/2;

% Plot two intensity maps per year
for j = 1:12
    % i is the week number
    i = (j-1)*26 + 1;
    subplot(3,4,j)
    yearstr = num2str(2004+floor(i/52));
    [temp,monthstr] = month((i - 52*floor(i/52))*7);
    mystr = [monthstr,' ',yearstr];
%     text(5,32,1000,mystr,'FontSize',14);
    DrawAFGmap
    hold on
    myintensity = 6.6*exp(mu + Covariates).*exp(sum(multiprod(Basis.phi,reshape(Estinfo(i).xestRTS + diag(Estinfo(i).PRTS)./2,1,1,Basis.nx)),3));
    % The 6.6 is due to the scale operations in the warp (2*3.3 = 6.6)
    myintensity(IN == 0)= 0.001; % Set a minimum for plotting purposes
    
    % log plot
    surfm(s2true,s1true,log10(myintensity))
    shading interp
    myscale = [0.99 300];
    caxis(myscale)
    colormap(c);
    h=colorbar;
    ticks_wanted=unique([0.1,1,10,100,1000]);
    caxis(log10(myscale))
    set(h,'YTick',log10(ticks_wanted));
    set(h,'YTickLabel',ticks_wanted);
    colorbar off
    title(mystr,'FontSize',17)
    set(gca,'XTick',[])
    set(gca,'YTick',[])
end
h = colorbar;
ticks_wanted=unique([1,10,100,1000]);
caxis(log10(myscale))
set(h,'YTick',log10(ticks_wanted));
set(h,'YTickLabel',ticks_wanted);
set(h,'Position',[0.91,0.1,0.015,0.8])
set(h,'Fontsize',15)
text(0.27,0.6,'Mean report count per week per degree^2','rotation',90,'FontSize',15)
set(gcf,'PaperPositionMode','auto')
print -dpng -r600  Yearlyintensity.png
close all

