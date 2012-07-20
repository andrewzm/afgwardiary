% Generate the lower, median and upper quantiles from the VB-Laplace
% results

% Find population density
Afg_pop = shaperead('../Shapefiles/afg_population_CSO_2011_2012','UseGeoCoords',true);
Pop_density = zeros(size(s1));
for i = 1:length(Afg_pop)
X = inpolygons(s1true(:),s2true(:),Afg_pop(i).Lon',Afg_pop(i).Lat');
Pop_density(X == 1) = log10(Afg_pop(i).T_P_BS/ Afg_pop(i).AREA_SQKM);
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

% Initialize
Uest = zeros(J2,J2,N); Umedian = zeros(J2,J2,N); Uupper = zeros(J2,J2,N); Ulower = zeros(J2,J2,N);
intensityest = zeros(J2,J2,N); int_median = zeros(J2,J2,N); int_upper = zeros(J2,J2,N); int_lower  = zeros(J2,J2,N);

% Find required quantiles
for i = 1:N
    Uest(:,:,i) = sum(multiprod(Basis.phi,reshape(Estinfo(i).xestRTS+diag(Estinfo(i).PRTS)./2,1,1,Basis.nx)),3);
    Umedian(:,:,i) = sum(multiprod(Basis.phi,reshape(Estinfo(i).xestRTS,1,1,Basis.nx)),3);
    Uupper(:,:,i) = sum(multiprod(Basis.phi,reshape(Estinfo(i).xestRTS + 1.65*sqrt(diag(Estinfo(i).PRTS)),1,1,Basis.nx)),3);
    Ulower(:,:,i) = sum(multiprod(Basis.phi,reshape(Estinfo(i).xestRTS-1.65*sqrt(diag(Estinfo(i).PRTS)),1,1,Basis.nx)),3);
    
    Covariates = Thetainfo(m).best(1)*Pop_density + Thetainfo(m).varb(1)*Pop_density.^2/2 ...
                 + Thetainfo(m).best(2)*Dist_to_city + Thetainfo(m).varb(2)*Dist_to_city.^2/2;
    
    intensityest(:,:,i) = exp(mu + Covariates+ Uest(:,:,i));
    int_median(:,:,i) = exp(mu + Covariates+Umedian(:,:,i));
    int_upper(:,:,i) = exp(mu + Covariates+Uupper(:,:,i));
    int_lower(:,:,i) = exp(mu + Covariates+Ulower(:,:,i));
end

% Save data (only needs to be run once per set of results)
save('Quantile_data','int_lower','int_median','int_upper','intensityest')