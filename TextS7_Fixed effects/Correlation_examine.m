function y = Correlation_examine(int,covariate,intensity,min_intensity)

% Check for association between covariate and intensity
%
% Inputs:   int (interval over which to check for association)
%           covariate (pop. density etc.)
%           intensity (nonparametric estimate of spatial intensity)
%           min_intensity (regions in space to consider have to have an intensity greater than min_intensity)

% Initialize
my_mean = zeros(size(int));
my_std = zeros(size(int));
my_lowerq = zeros(size(int));
my_upperq = zeros(size(int));

% For each interval in domain
for i = 1:length(int)
    if i > 1
      Selection = (abs(covariate)<int(i))&(abs(covariate)>=int(i-1))&(intensity > min_intensity);
    else
      Selection = (abs(covariate)<int(i))&(intensity > min_intensity);
    end
    % Find intensity mean and std corresponding to this interval of the fixed effect
    my_mean(i) =mean(log(intensity(Selection)));
    my_std(i) = std(log(intensity(Selection)));%/sqrt(sum(sum(Selection)));
    my_lowerq(i) = quantile(log(intensity(Selection)),0.25);
    my_upperq(i) = quantile(log(intensity(Selection)),0.75);
end

% Plot results
plot(int(1:end),my_mean(1:end),'k'); hold on; plot([int(1), int(end)],[mean(mean(log(intensity(intensity > min_intensity)))),mean(mean(log(intensity(intensity > min_intensity))))],'k--')
plot(int(1:end),my_mean(1:end)+ my_std(1:end),'r--')
plot(int(1:end),my_mean(1:end)- my_std(1:end),'r--')
set(gca,'Xlim',[int(1), int(end)]);

% plot(int(2:end),my_upperq(2:end),'r--')
% plot(int(2:end),my_lowerq(2:end),'r--')