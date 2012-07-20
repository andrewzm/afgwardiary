% Find distance to Pakistan
Pak_border = shaperead('Pakistan/PAK_adm0','UseGeoCoords',true);
Dist_to_Pak = zeros(size(GRID));
Pak_lon = Pak_border.Lon(1:end);
Pak_lat = Pak_border.Lat(1:end);
for i = 1:size(GRID,1)
    for j = 1:size(GRID,2)
%         Dist_to_city(i,j) = min(my_dist(s1(j),s2(i),Filtered_cities(:,1),Filtered_cities(:,2)));
        Dist_to_Pak(i,j) = min(posdist(s1(j),s2(i),Pak_lon,Pak_lat));
    end
    i
end
