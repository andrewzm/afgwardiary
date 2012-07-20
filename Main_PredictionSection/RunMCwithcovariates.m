% Filename: RunMCwithcovariates.m
% Author: Andrew Zammit Mangion
% Date: November 2011
% Description: Runs Monte Carlo runs with inferred model (see paper for
% details) to predict AOG activity in 2010.
%
%Requires: Mapping Toolbox
%          Statistical Toolbox
%          LocalisedKernelPhi.m
%          ../AfghanModelHighResWholePrec_Covariates (main data)
%          ../Shapefiles
%          inpolygons.m (available from MATLAB Exchange)
%          multiprod.m (available from MATLAB Exchange)
%          posdist.m (available from MATLAB Exchange)

% Load results from VB-Laplace
studydata = '../AfghanModelHighResWholePrec_Covariates'; 
load(studydata)
DisturbVar = inv(diag(diag(Thetainfo(m).Meanprecmat))); % Adopt marginal variance (precision is diagonally dominant)

N= 2000;                            % Number of samples
xsample2009 = zeros(Basis.nx,52);   % Samples for 52 weeks, state in 2009
fieldsample2009 = zeros(J,J,52);    % Sample field in 2009
intsample2009 = zeros(J,J,52);      % Sample intensity in 2009
xsample2010 = zeros(Basis.nx,52);   % Sample state in 2010
fieldsample2010 = zeros(J,J,52);    % Sample field in 2010
intsample2010 = zeros(J,J,52);      % Sample intensity in 2010

% Load Provincial Boundaries
Provbounds = shaperead('../Shapefiles/admin2_poly_32.shp','UseGeoCoords',true);

% Setup space
strue1 = s./scale(1) + shift(1);
strue2 = s./scale(2) + shift(2);
[s1true,s2true] = meshgrid(strue1,strue2);

% Provincial masks
Provmask = zeros(J^2,length(Provbounds));
for Pnum = 1:length(Provbounds)
    Provmask(:,Pnum) = inpolygons(s1true(:),s2true(:),Provbounds(Pnum).Lon',Provbounds(Pnum).Lat');
    Provname(Pnum).str = Provbounds(Pnum).PRV_NAME;
end
Provintensity2009 = zeros(length(Provbounds),N);    % Provincial intensity in 2009
Count2009 = zeros(length(Provbounds),N);            % Provincial count in 2009
Provintensity2010 = zeros(length(Provbounds),N);    % Provincial intensity in 2010
Count2010 = zeros(length(Provbounds),N);            % Provincial count in 2010
Growth = zeros(length(Provbounds),N);               % Sample provincial growth


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

% Fixed effects on intensity field
Covariates = Thetainfo(m).best(1)*Pop_density + Thetainfo(m).varb(1)*Pop_density.^2/2 ...
                 + Thetainfo(m).best(2)*Dist_to_city + Thetainfo(m).varb(2)*Dist_to_city.^2/2;

% Main Monte Carlo routine for intensity in each province            
for i = 1:N
    for j = 261:312
        %Sample x
        xsample2009(:,j-260) = mvnrnd(Estinfo(j).xestRTS,diag(diag(Estinfo(j).PRTS)));
        fieldsample2009(:,:,j-260) = sum(multiprod(Basis.phi,reshape(xsample2009(:,j-260),1,1,Basis.nx)),3);
        intsample2009(:,:,j-260) = exp(fieldsample2009(:,:,j-260) + Thetainfo(m).muest + Covariates);
    end
    TotIntensity2009 = sum(intsample2009,3);
    for Pnum = 1:length(Provbounds)
        Provintensity2009(Pnum,i) = sum(sum(reshape(TotIntensity2009(:).*Provmask(:,Pnum),length(s),[])))*ds^2;
        Count2009(Pnum,i) = poissrnd(Provintensity2009(Pnum,i));
    end
    xsample2010(:,1) = xsample2009(:,end);
    fieldsample2010(:,:,1) = fieldsample2009(:,:,end);
    intsample2010(:,:,1) = intsample2009(:,:,end);
    for j = 2:52
        xsample2010(:,j) = xsample2010(:,j-1) + Thetainfo(m).thetaest + mvnrnd(zeros(Basis.nx,1),DisturbVar)';
        fieldsample2010(:,:,j) = sum(multiprod(Basis.phi,reshape(xsample2010(:,j),1,1,Basis.nx)),3);
        intsample2010(:,:,j) = exp(fieldsample2010(:,:,j) + Thetainfo(m).muest +  Covariates);
    end
    TotIntensity2010 = sum(intsample2010,3);
    for Pnum = 1:length(Provbounds)
        Provintensity2010(Pnum,i) = sum(sum(reshape(TotIntensity2010(:).*Provmask(:,Pnum),length(s),[])))*ds^2;
        Count2010(Pnum,i) = poissrnd(Provintensity2010(Pnum,i));
        Growth(Pnum,i) = Count2010(Pnum,i)/Count2009(Pnum,i);
    end
    i
end

% Save results
save('GrowthtrialN2000WholePrecWithDistCovariates')
        