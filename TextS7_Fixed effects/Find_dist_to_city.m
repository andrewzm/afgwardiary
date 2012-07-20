% Find distance to city
Afg_cities = shaperead('Settlements/07_03_settlements','UseGeoCoords',true);
Filtered_cities = 0; j =1;
% Only consider "important" cities
for i = 1:length(Afg_cities)
   if Afg_cities(i).TYPE == 2 ||  Afg_cities(i).TYPE == 3
       Filtered_cities(j,1) = Afg_cities(i).LON;
       Filtered_cities(j,2) = Afg_cities(i).LAT;
       j = j+1;
   end
end
Dist_to_city = zeros(size(GRID));
for i = 1:size(GRID,1)
    for j = 1:size(GRID,2)
%         Dist_to_city(i,j) = min(my_dist(s1(j),s2(i),Filtered_cities(:,1),Filtered_cities(:,2)));
        Dist_to_city(i,j) = min(posdist(s1(j),s2(i),Filtered_cities(:,1),Filtered_cities(:,2)));
    end
end
