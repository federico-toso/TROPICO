function [distance] = distance_haversine (lat1, lon1, lat2, lon2, radius)

%{
INPUT:
* latX, lonX = 2-column vector of coordinates [rad]
* radius = of the planet [m]

OUPUT:
* distance = great circle distances between the input points, the matrix 
has 1 row per each latlon1, 1 column per each latlon2

%}

distance = zeros (size(lat1,1),size(lat2,1));
for i = 1 : size(lat1,1)
    for j = 1:size(lat2,1)
        delta_lat = lat2(j) - lat1(i);                            % difference in latitude
        delta_lon = lon2(j) - lon1(i);                            % difference in longitude

        a = sin(delta_lat/2)^2 + cos(lat1(i)) * cos(lat2(j)) * sin(delta_lon/2)^2;
        c = 2 * atan2(sqrt(a), sqrt(1-a));

        distance(i,j) = radius * c;                            	% distance in km
    end
end