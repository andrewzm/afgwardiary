% Find population density
Afg_pop = shaperead('../Shapefiles/afg_population_CSO_2011_2012','UseGeoCoords',true);
Pop_density = zeros(size(GRID));
for i = 1:length(Afg_pop)
    X = inpolygons(S1(:),S2(:),Afg_pop(i).Lon',Afg_pop(i).Lat');
    Pop_density(X == 1) = Afg_pop(i).T_P_BS/ Afg_pop(i).AREA_SQKM;
end