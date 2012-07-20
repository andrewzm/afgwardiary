function AfghanCorrelationTests(filename)

% Filename: AfghanCorrelationTests.m
% Author: Andrew Zammit Mangion
% Date: June 2011
%
% Description: AFG PACF and PCCF analysis. This file does not need to be
% run - the output is available in CorrelationTests.mat
%
% Requires: ../AfghanDataAllDay (AWD in Matlab format)
%           ../Common_functions (for PACF/PCCF)
%
% Input: filename (output file name)

addpath('../Common_functions')

%Setup time
dt = 1 ;  
N = 302;
t = 0:dt:dt*(N-1);  

%Setup space
J = 101;
s = linspace(0,36,J);

load('../AfghanDataAllDay')
% Arrange into weeks
numofweeks = ceil(length(spikeAll)/7);
Weeknum = zeros(numofweeks,1);
shift = [57.5,28.5];  % Scale and shift so AFG roughly fits into a 36 x 36 square with bottom LH corner on origin
scale = [2,3.3];
i=1;
k=1;
% hold on
while i < length(spikeAll)
   spikeAllWeekData(k).Coords = [];
   for j = 0:6  %Cluster in weeks and put on 30x30 square
    if length(spikeAll(i).Coords) > 0
        spikeAll(i).Coords(:,1) = (spikeAll(i).Coords(:,1) - shift(1))*scale(1);
        spikeAll(i).Coords(:,2) = (spikeAll(i).Coords(:,2) - shift(2))*scale(2);
        spikeAllWeekData(k).Coords = [spikeAllWeekData(k).Coords;spikeAll(i).Coords];
    end
    spikeAllWeekData(k).Coords = CleanfromNans(spikeAllWeekData(k).Coords);
    spikeAllWeekData(k).Coords = CleanOutliers(spikeAllWeekData(k).Coords);
    Daynum(i) = size(spikeAll(i).Coords,1);
    Weeknum(k) = Weeknum(k) + length(spikeAll(i).Coords);
    i = i+1;
   end
   hold on
   k = k+1;
end
Weeknum(end) = [];
N = length(Weeknum);

%Pair correlation estimator (Moller Pg. 44, Brix Diggle Pg 828 Stoyan Pg 284)
r = 0.1:0.5:7;
Aest = zeros(N,1);
lkernelest = zeros(N,1);
gest = zeros(length(r),N);
for i = 1:N
    [Aest(i),lkernelest(i),gest(:,i)] = PairCorrFunc(r,spikeAllWeekData(i).Coords,s);
    i
end

%Pair cross-correlation estimator
Acrossest = zeros(N,1);
lkernelcrossest = zeros(N,1);
gcrossest = zeros(length(r),N);
for i = 2:N
    [Acrossest(i),lkernelcrossest(i),gcrossest(:,i)] = PairCrossCorrFunc(r,spikeAllWeekData(i).Coords,spikeAllWeekData(i-1).Coords,s);
    i
end

save(filename)

function spikes = CleanfromNans(spikes)
% Remove NaNs from events in dataset
k=1;
if ~isempty(spikes)
    while k <= size(spikes,1)
        if isnan(spikes(k,1)) || isnan(spikes(k,2))
            spikes(k,:) = [];
        else
            k = k+1;
        end
    end
end

function spikes = CleanOutliers(spikes)
% Create a rough polygon around AFG from which to exclude spikes (this
% function may be replaced with inpolygon for more accurate removal)
k=1;
if ~isempty(spikes)
    while k <= size(spikes,1)
        if spikes(k,1) < 0 || spikes(k,1) > 30 || spikes(k,2)<0 || spikes(k,2) > 30 ...
                || spikes(k,2) - 1.75*spikes(k,1) < -30
            spikes(k,:) = [];
        else
            k = k+1;
        end
    end
end

