% GenerateIntensityVideo.m
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


clear all
close all

% Load results
studydata = '../AfghanModelHighResWholePrec_Covariates';
load(studydata)
mu = Thetainfo(m).muest;

%Upsample Basis for better resolution and setup space
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


%---------------------------------------
% MOVIE OF INTENSITY
%---------------------------------------
% Find fixed effects
Covariates = Thetainfo(m).best(1)*Pop_density + Thetainfo(m).varb(1)*Pop_density.^2/2 ...
    + Thetainfo(m).best(2)*Dist_to_city + Thetainfo(m).varb(2)*Dist_to_city.^2/2;
% Multiple frames for decreased frame rate
FRMult = 6;
nframes = size((length(Estinfo)-1)*FRMult);  %Defines length of animation
M = moviein(nframes); close
c = flipud(hot(1024));
c(end-120:end,:) = [];
count=1;

screen_size = get(0, 'ScreenSize');
f1 = figure;
set(f1, 'Position', [0 60 screen_size(3)-200 screen_size(4) ] );
weekinter = 365/7;
for i = 2:length(Estinfo)-1
    DrawAFGmap
    i
    hold on
    myintensity = 6.6*exp(mu + Covariates).*exp(sum(multiprod(Basis.phi,reshape(Estinfo(i).xestRTS + diag(Estinfo(i).PRTS)./2,1,1,Basis.nx)),3));
    % The 6.6 is a result of the affine warping (2*3.3 = 6.6)
    myintensity(IN == 0)= 0.001; % Lower bound for plotting purposes
    
    % log plot
    surfm(s2true,s1true,log10(myintensity))
    shading interp
    myscale = [1 50];
    caxis(myscale)
    colormap(c);
    h=colorbar;
    ticks_wanted=unique([1,10,20,30,50]);
    caxis(log10(myscale))
    set(h,'YTick',log10(ticks_wanted));
    set(h,'YTickLabel',ticks_wanted);
    set(h,'FontSize',20)
    set(h,'Position',[0.94,0.1,0.02,0.8])
    
    yearstr = num2str(2004+floor(i/weekinter));
    [temp,monthstr] = month((i - round(weekinter*floor(i/weekinter)))*7);
    mystr = ['\bf ',monthstr,' ',yearstr];
    textm(38.4,61,mystr,'FontSize',20);
    points = [spikeAllWeekData(i).Coords(:,1)/scale(1)+shift(1) spikeAllWeekData(i).Coords(:,2)/scale(2) + shift(2)];
    Inpoints = inpolygon(points(:,1),points(:,2),Countrybounds.Lon',Countrybounds.Lat'); % Only consider points inside AFG
    if sum(Inpoints == 0) > 0
        points(Inpoints == 0,:) = [];
    end
    plot3m(points(:,2),points(:,1),repmat(1e10,length(points),1),'g.','Linewidth',0.5)
    
    
    setm(ax,'Grid','on','Frame','on','meridianlabel','on','parallellabel','on','FontSize',14)
    set(gca,'FontSize',20)
    xlabel('Longitude','FontSize',25)
    ylabel('Latitude','FontSize',25)
    text(-0.165,0.5,'Expected report count per week per degree^2','rotation',90,'FontSize',16)
    
    set(gcf,'color',[1,1,1])
    
    hold off
    clear Editedframe Thisframe
    
    Thisframe = getframe(gcf);
    Editedframe = Thisframe;
   
    % Shift frame horizontally
    Editedframe.cdata(:,:,1) = circshift(Thisframe.cdata(:,:,1),[0,-50]);
    Editedframe.cdata(:,:,2) = circshift(Thisframe.cdata(:,:,2),[0,-50]);
    Editedframe.cdata(:,:,3) = circshift(Thisframe.cdata(:,:,3),[0,-50]);
    M(:,count:count+FRMult-1)=Editedframe;
    
    count = count+FRMult;
    pause(0.1)
    colorbar off
    i
end
save movie.mat M
% mov=close(mov); % closes the mov
mpgwrite(M,jet,'./VideoS1.mpg'); % Convert the movie to MPEG format

