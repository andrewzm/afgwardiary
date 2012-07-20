% Filename: Spatial_nonparametric.m
% Author: Andrew Zammit Mangion
% Date: April 2012
%
% Description: Creates a spatial intensity map of the AWD
%
% Requirements: 
%               FindVolume.m (finds volume of kernel)
%               ../AfghanDataAllDay.mat (AWD in Matlab format)


clear all

% Spatial map
load('../AfghanDataAllDay')

% Nonparametric kernel
sigmas1 = 0.08;
sigmas2 = 0.08;
GaussKernel = @(s1,s2,mu1,mu2) exp(-((s1-mu1).^2)./(2*sigmas1^2) - ((s2-mu2).^2)./(2*sigmas2^2));

N = length(spikeAll);
s1_fine = -2:0.01:2;  % Spatial quadrature points
s2_fine = s1_fine;


disp('Calculating default volume...')
Volume = FindVolume(s1_fine,s2_fine,GaussKernel);

% Setup space and time
s1 = 60:0.1:75;
s2 = 29:0.1:39;
[S1,S2] = meshgrid(s1,s2);
t = 1:N;

% Initialize
GRID = zeros(length(s2),length(s1));

for i = 1:N
    if ~isempty(spikeAll(i).Coords)
    spikeAll(i).Coords(isnan(spikeAll(i).Coords(:,1)),:) = []; % remove NaNs
    for j = 1:length(spikeAll(i).Coords(:,1))  % Place kernel at each point (ignore boundary effects for this map)
            GRID = GRID(:,:) + GaussKernel(S1,S2,spikeAll(i).Coords(j,1),spikeAll(i).Coords(j,2))/Volume;
    end
    end
    i
end

save('Spatial_nonparametric.mat')
