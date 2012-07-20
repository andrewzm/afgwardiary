% Filename: Afg_Kernel_Estimation.m
% Author: Andrew Zammit Mangion
% Date: April 2012
% Description: Nonparametric intensity estimation for the AWD
% 
% Requires: FindSTVolume.m
%           ST_Kernek_estimation.m
%           ../AfghanDataAllDay.mat (AWD in Matlab format)

clear all
load('../AfghanDataAllDay')

%Sort into weeks
numofweeks = ceil(length(spikeAll)/7);
Weeknum = zeros(numofweeks,1);
shift = [57.5,28.5]; % Scale and warp such that Afghanistan is approx. on a 36 x 36 square
scale = [2,3.3];
i=1;
k=1;

while i < length(spikeAll)
   spikeAllWeekData(k).Coords = [];
   for j = 0:6  %Cluster in weeks and put on 30x30 square
    if length(spikeAll(i).Coords) > 0
        spikeAll(i).Coords(:,1) = (spikeAll(i).Coords(:,1) - shift(1))*scale(1);
        spikeAll(i).Coords(:,2) = (spikeAll(i).Coords(:,2) - shift(2))*scale(2);
        spikeAllWeekData(k).Coords = [spikeAllWeekData(k).Coords;spikeAll(i).Coords];
    end
    Daynum(i) = size(spikeAll(i).Coords,1);
    Weeknum(k) = Weeknum(k) + length(spikeAll(i).Coords);
    i = i+1;
   end
   hold on
   k = k+1;
end

Weeknum(end) = [];

% Kernel setup
sigmas1 = 0.5;
sigmas2 = 0.5;
sigmat =  2;
s1 = 0:0.1:36;
s2 = 0:0.1:36;

% Estimate ST intensity using nonparametric methods
GRID = ST_Kernel_estimation(spikeAllWeekData,sigmas1,sigmas2,sigmat,s1,s2,'F');

% Save data
save('AWD_Nonparametric_smooth','GRID')