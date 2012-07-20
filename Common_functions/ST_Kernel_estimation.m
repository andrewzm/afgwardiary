function GRID = ST_Kernel_estimation(spikes,sigmas1,sigmas2,sigmat,s1,s2,spatcorrection)

% Filename: ST_Kernel_estimation.m
% Author: Andrew Zammit Mangion
% Date: April 2012
% Description: Nonparametric intensity estimation for spatiotemporal
% systems with ST boundary correction
% 
% Inputs:   spikes (set of x/y coordinates)
%           sigmas1 (std of kernel in s1)
%           sigmas2 (std of kernel in s2)
%           sigmas2 (std of kernel in t)
%           s1 (1D vector of boundary box in s1)
%           s2 (1D vector of boundary box in s2)
%           spatcorrection = 'T' or 'F' (enable or disable boundary correction)
%
% Requires: FindSTVolume.m
         
GaussKernel = @(s1,s2,t,mu1,mu2,mu3) exp(-((s1-mu1).^2)./(2*sigmas1^2) - ((s2-mu2).^2)./(2*sigmas2^2) -((t-mu3).^2)./(2*sigmat^2));

% Fine grid for numerical integration
s1_fine = linspace(-3*sigmas1,3*sigmas1,61);
s2_fine = linspace(-3*sigmas2,3*sigmas2,61);
t_fine = linspace(-3*sigmat,3*sigmat,1201);

N = length(spikes);
disp('Calculating default volume...')
Volume = FindSTVolume(s1_fine,s2_fine,t_fine,GaussKernel,[0 0 0]);

% Initialize
GRID = zeros(length(s1),length(s2),N);
[S1,S2] = meshgrid(s1,s2);

if spatcorrection == 'F'
    % Do not correct at spatial boundaries
    for i = 1:N
        disp(strcat('Smoothing out week ',num2str(i)));
        spikes(i).Coords(isnan(spikes(i).Coords(:,1)),:) = []; % Remove NaNs
        
        % Temporal edge correction
        if i <= max(t_fine)
            disp('Volume correction...')
            Volume_true = FindSTVolume(s1_fine,s2_fine,t_fine(find(t_fine == 1-i):end),GaussKernel,[0 0 0]);
        elseif i >= N - max(t_fine)
            disp('Volume correction...')
            Volume_true = FindSTVolume(s1_fine,s2_fine,t_fine(1:find(t_fine == N-i)),GaussKernel,[0 0 0]);
        else
            Volume_true = Volume;
        end
        
        for j = 1:length(spikes(i).Coords(:,1))
            for k = max(i-t_fine(end),1):min(i+t_fine(end),N)
                GRID(:,:,k) = GRID(:,:,k) + GaussKernel(S1,S2,k,spikes(i).Coords(j,1),spikes(i).Coords(j,2),i)/Volume_true;
            end
        end
    end
end


if spatcorrection == 'T'
    % Correct at both space and time boundaries
    for i = 1:N
        disp(strcat('Smoothing out week ',num2str(i)));
        spikes(i).Coords(isnan(spikes(i).Coords(:,1)),:) = [];
        
        for j = 1:length(spikes(i).Coords(:,1))
            if (i <= max(t_fine) || i >= N - max(t_fine) || ...
                    spikes(i).Coords(j,1) <= max(s1_fine) || spikes(i).Coords(j,1) >= s1(end) - max(s1_fine) || ...
                    spikes(i).Coords(j,2) <= max(s2_fine) || spikes(i).Coords(j,2) >= s2(end) - max(s2_fine))
                
                s1_fine_temp = linspace(max(s1(1), spikes(i).Coords(j,1)-3*sigmas1),...
                    min(s1(end), spikes(i).Coords(j,1) + 3*sigmas1),21);
                s2_fine_temp = linspace(max(s2(1), spikes(i).Coords(j,2)-3*sigmas2),...
                    min(s2(end), spikes(i).Coords(j,2) + 3*sigmas2),21);
                t_fine_temp = linspace(max(1, i-3*sigmat),...
                    min(N, i + 3*sigmat),101);
                Volume_true = FindSTVolume(s1_fine_temp,s2_fine_temp,t_fine_temp,GaussKernel,...
                               [spikes(i).Coords(j,1),spikes(i).Coords(j,2),i]);
%                 disp(strcat(num2str(j),':Volume correction...'))
            else
                
                Volume_true = Volume;
                
            end
            for k = max(i-5,1):min(i+5,N)
                GRID(:,:,k) = GRID(:,:,k) + GaussKernel(S1,S2,k,spikes(i).Coords(j,1),spikes(i).Coords(j,2),i)/Volume_true;
            end
        end
    end
end
